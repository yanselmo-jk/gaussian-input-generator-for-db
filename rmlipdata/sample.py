"""Structure samplers: monomers, pairs, clusters, bond scans, pseudo-TS reactions, cation hops, nanoreactor.

`rotors`, `bends`, `min_far_distance`, `sample`, `_forcefield` copied from
7_apply_distortion_to_molecule/src/moldistort/geometry.py.
"""
import warnings
import zlib
from itertools import combinations

import numpy as np
from ase import Atoms
from ase.data import atomic_masses, atomic_numbers, covalent_radii
from ase.io import read as ase_read
from rdkit import Chem
from rdkit.Chem import rdMolTransforms as T

from . import build, species, xtb
from .g16 import n_electrons

BOHR_A = 0.529177210903
CM_PER_SQRT_EH_BOHR2_AMU = 5140.49
DONORS = {"O", "N", "F"}
ELECTROPHILES = {"P", "S", "B"}
LIMIT = 10 ** 9
HOP_MAX_DEV_A = 1.0  # hop frames whose cation is farther than this from the intended path point are dropped


# ---- distortion helpers (moldistort) --------------------------------------------------------------
# from 7_apply_distortion_to_molecule/src/moldistort/geometry.py
def rotors(mol):
    """(i, j, k, l) for every non-ring single bond j-k with both ends non-terminal."""
    out = []
    for b in mol.GetBonds():
        j, k = b.GetBeginAtom(), b.GetEndAtom()
        if b.GetBondType() != Chem.BondType.SINGLE or b.IsInRing() or j.GetDegree() < 2 or k.GetDegree() < 2:
            continue
        i = next(n.GetIdx() for n in j.GetNeighbors() if n.GetIdx() != k.GetIdx())
        l = next(n.GetIdx() for n in k.GetNeighbors() if n.GetIdx() != j.GetIdx())
        out.append((i, j.GetIdx(), k.GetIdx(), l))
    return out


def bends(mol):
    """(i, j, k) angles that RDKit can set: the moved side j-k must not be a ring bond."""
    out = []
    for a in mol.GetAtoms():
        j = a.GetIdx()
        for i, k in combinations([n.GetIdx() for n in a.GetNeighbors()], 2):
            if not mol.GetBondBetweenAtoms(j, k).IsInRing():
                out.append((i, j, k))
            elif not mol.GetBondBetweenAtoms(j, i).IsInRing():
                out.append((k, j, i))
    return out


def min_far_distance(mol, xyz, topo=None):
    """Shortest distance between atoms more than 2 bonds apart (clash check)."""
    far = (Chem.GetDistanceMatrix(mol) if topo is None else topo) > 2
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    return d[far].min(initial=np.inf)


def sample(mol, n, torsion=True, angle_delta=10.0, clash_min=1.0, rng=None):
    """Frame 0 = input geometry, then up to n random torsion/angle distortions without clashes."""
    rng = rng or np.random.default_rng(0)
    m = Chem.Mol(mol)
    conf = m.GetConformer()
    pos0 = conf.GetPositions()
    tors = rotors(m) if torsion else []
    angs = [(a, T.GetAngleDeg(conf, *a)) for a in bends(m)] if angle_delta else []
    angs = [(a, eq) for a, eq in angs if eq < 175.0]
    frames = [pos0]
    if not tors and not angs:
        warnings.warn("no rotatable torsion or settable angle; only the input geometry is returned")
        return frames
    topo = Chem.GetDistanceMatrix(m)
    for _ in range(50 * n):
        if len(frames) > n:
            break
        conf.SetPositions(pos0)
        for t in tors:
            T.SetDihedralDeg(conf, *t, rng.uniform(-180.0, 180.0))
        for a, eq in angs:
            T.SetAngleDeg(conf, *a, eq + rng.uniform(-angle_delta, angle_delta))
        xyz = conf.GetPositions()
        if min_far_distance(m, xyz, topo) >= clash_min:
            frames.append(xyz)
    return frames


def _forcefield(mol):
    from rdkit.Chem import AllChem
    if AllChem.MMFFHasAllMoleculeParams(mol):
        return AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    return AllChem.UFFGetMoleculeForceField(mol)


# ---- common --------------------------------------------------------------------------------------
def _limit(scfg):
    return int(scfg.get("max_structures", LIMIT))


def _uhf(symbols, charge):
    return n_electrons(symbols, int(charge)) % 2


def _timeout(ctx):
    return int(ctx.cfg.get("xtb", {}).get("timeout_s", 1800))


def _formal(ctx, names):
    return sum(ctx.cat[n].charge for n in names)


def _cation_only(ctx, names):
    """True if every member is a cation (a bare Li+/Na+ cluster)."""
    return all(ctx.cat[n].kind == "cation" for n in names)


def _multi_radical(ctx, names):
    """Two or more radical members (catalog radicals: EC_ro, MeOCO2, CH3, CH3O): the even-electron singlet is written
    as a UKS guess=mix (BS) singlet. RKS pairs the two electrons, which for separated radicals gives a state several
    eV above the open-shell singlet that the same fragments get elsewhere (bond_scan BS)."""
    return sum(any(a.GetNumRadicalElectrons() for a in m.GetAtoms())
               for m in (Chem.MolFromSmiles(ctx.cat[n].smiles or "") for n in names) if m is not None) >= 2


def _bs(ctx, names, delta_q):
    """BS singlet: two or more radical members, or delta_q <= -2 (_q_ok: one Li/Na per added electron, i.e. two
    reduced metal centres)."""
    return int(delta_q) <= -2 or _multi_radical(ctx, names)


def _rec(ctx, symbols, positions, charge, sampler, names, job, bs, params):
    return dict(symbols=[str(s) for s in symbols], positions=np.asarray(positions, float), charge=int(charge),
                sampler=sampler, composition=build.composition(list(names)), job=job, bs=bool(bs),
                params=params)


def _with_positions(mol, positions):
    """Copy of the topology mol carrying the given conformer."""
    m = Chem.Mol(mol)
    m.GetConformer().SetPositions(np.asarray(positions, float))
    return m


def _q_ok(ctx, q, names=None, delta_q=0):
    """Total charge inside cfg charge.range (default [-2, 2]); beyond it gas-phase clusters lose the extra
    electron (anions) or Coulomb-explode (cations). With names, a reduction (delta_q < 0) is only a real
    bound state when the cluster holds Li/Na (down to q = -1) or a neutral radical (whose extra electron
    closes a shell)."""
    lo, hi = ctx.cfg.get("charge", {}).get("range", (-2, 2))
    if not lo <= q <= hi:
        return False
    if names is not None and int(delta_q) <= -2:
        n_li = sum(1 for n in names if ctx.cat[n].kind == "cation"
                   or "Li" in ctx.cat[n].smiles or "Na" in ctx.cat[n].smiles)
        if n_li < -int(delta_q):
            return False  # a multi-electron reduction needs one Li/Na per added electron
    if names is None or int(delta_q) >= 0:
        return True
    if any(ctx.cat[n].kind == "cation" or "Li" in ctx.cat[n].smiles or "Na" in ctx.cat[n].smiles
           for n in names):
        return q >= -1  # Li + 2 PF6 at q=-2 is Li0 + two PF6-: an unbound gas-phase dianion
    if q >= -1:
        for n in names:
            sp = ctx.cat[n]
            if sp.charge != 0:
                continue
            mol = Chem.MolFromSmiles(sp.smiles)
            if mol is not None and any(a.GetNumRadicalElectrons() > 0 for a in mol.GetAtoms()):
                return True  # e.g. CH3O + e- -> CH3O-
    return False


def _too_close(symbols, x, f, mask=None):
    """True if any atom pair is closer than max(f x covalent-radius sum, 0.6 A) (an absolute clash floor)."""
    x = np.asarray(x, float)
    r = np.array([_cov(symbols, i) for i in range(len(symbols))])
    d = np.linalg.norm(x[:, None] - x[None], axis=-1) + np.eye(len(x)) * 1e9
    hit = d < np.maximum(f * (r[:, None] + r[None]), 0.6)
    return bool((hit if mask is None else hit & mask).any())


def _cluster(ctx, combo, delta_q):
    q = _formal(ctx, combo) + int(delta_q)
    if _cation_only(ctx, combo) and int(delta_q) > 0:
        return q, None  # no extra oxidation of a bare cation cluster
    key = (tuple(combo), q)
    if not _q_ok(ctx, q, combo, delta_q):
        return q, None
    if key not in ctx.clusters:
        rng = np.random.default_rng(int(ctx.cfg.get("seed", 0)) + zlib.crc32(f"{key}".encode()))
        cl = build.assemble(combo, ctx.cat, ctx.workdir, q, rng)
        if _fragmented(cl):  # e.g. a +3 gas-phase cluster that Coulomb-exploded during the xTB opt
            warnings.warn(f"cluster {list(combo)} q={q} fell apart in the xtb opt; skipped")
            cl = None
        ctx.clusters[key] = cl
    return q, ctx.clusters[key]


def _fragmented(cl, cutoff=4.0):
    """True if some member has no atom within cutoff A of any other member."""
    return build.fragmented(cl["positions"], [idx for _, idx in cl["members"]], cutoff)


def _etemp(*results):
    """Provenance param when an xTB step behind a record needed the --etemp retry (SCC not converged)."""
    return {"xtb_etemp_K": xtb.ETEMP_RETRY_K} if any(r and r.get("etemp_K") for r in results) else {}


def _stream(ctx, sampler, names, q):
    """Per-(sampler, combo, q) RNG: an xtb failure on one combo must not shift later combos' draws."""
    return np.random.default_rng(int(ctx.cfg.get("seed", 0)) + zlib.crc32(f"{sampler}{list(names)}{q}".encode()))


def _combos(ctx):
    cfg = ctx.cfg
    if not hasattr(ctx, "combos") or ctx.combos is None:
        rng = np.random.default_rng(int(cfg.get("seed", 0)))
        ctx.combos = build.combos(ctx.cat, cfg["species"], cfg.get("cluster", {}),
                                  int(cfg.get("limits", {}).get("max_atoms", 10 ** 9)), rng)
    return ctx.combos


def _mols(ctx, combo):
    return [(n, species.geometry(ctx.cat[n], ctx.workdir)[2]) for n in combo]


def _cov(symbols, i):
    return float(covalent_radii[atomic_numbers[symbols[i]]])


def _bonds(symbols, positions, atoms=None, scale=1.2):
    """Covalent bonds (i < j) shorter than scale x covalent-radius sum; Li/Na contacts are ionic, so ignored."""
    pos = np.asarray(positions, float)
    idx = list(range(len(symbols))) if atoms is None else list(atoms)
    out = set()
    for a, i in enumerate(idx):
        for j in idx[a + 1:]:
            if symbols[i] in ("Li", "Na") or symbols[j] in ("Li", "Na"):
                continue
            if np.linalg.norm(pos[i] - pos[j]) < scale * (_cov(symbols, i) + _cov(symbols, j)):
                out.add((min(i, j), max(i, j)))
    return out


def _components(symbols, positions):
    """Covalently connected atom groups (per _bonds); Li/Na are single-atom groups."""
    parent = list(range(len(symbols)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in _bonds(symbols, positions):
        parent[root(i)] = root(j)
    groups = {}
    for i in range(len(symbols)):
        groups.setdefault(root(i), []).append(i)
    return list(groups.values())


def _com(pos, idx):
    return pos[list(idx)].mean(axis=0)


def _min_between(pos, ia, ib):
    return float(np.linalg.norm(pos[list(ia)][:, None] - pos[list(ib)][None], axis=-1).min())


def _shift_to(pos, ia, ib, u, target, span):
    """Copy of pos with member B moved along u so that the shortest A-B distance is target (bisection)."""
    lo, hi = -span, 3.0 * span
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        moved = pos.copy()
        moved[ib] += mid * u
        lo, hi = (mid, hi) if _min_between(moved, ia, ib) < target else (lo, mid)
    moved = pos.copy()
    moved[ib] += hi * u
    return moved


def _until_jump(E, tol=5.0):
    """Number of leading scan points before the first consecutive xTB energy jump > tol eV (SCC state change)."""
    for k in range(1, len(E)):
        if abs(E[k] - E[k - 1]) > tol:
            warnings.warn(f"xtb scan energy jump {E[k] - E[k - 1]:+.1f} eV at point {k}; later points dropped")
            return k
    return len(E)


def _until_apart(cl, frames):
    """Number of leading scan frames before a member drifts > 4 A from the rest (relaxed scans let e.g. Li+
    leave a stretched EC radical cation in a q=+2 cluster)."""
    for k, f in enumerate(frames):
        if _fragmented({"members": cl["members"], "positions": np.asarray(f, float)}):
            return k
    return len(frames)


def _kabsch(x, ref, idx):
    """x rigidly moved so that atoms idx best overlay ref."""
    a, b = x[idx] - x[idx].mean(0), ref[idx] - ref[idx].mean(0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    return (x - x[idx].mean(0)) @ rot + ref[idx].mean(0)


# ---- samplers ------------------------------------------------------------------------------------
def monomer(ctx, scfg, delta_q):
    """One sp for cations; opt + classical normal-mode + moldistort frames for molecules and anions."""
    n = 0
    for kind in ("cations", "anions", "solvents"):
        for name in ctx.cfg["species"].get(kind, []):
            getattr(ctx, "progress", lambda *a: None)(name)
            if n >= _limit(scfg):
                return
            sp = ctx.cat[name]
            q = sp.charge + int(delta_q)
            if not _q_ok(ctx, q, [name], delta_q):
                continue
            if sp.kind == "cation":
                if q > 1:  # bare Li2+/Na2+ is not an electrolyte state
                    continue
                yield dict(_rec(ctx, [name], np.zeros((1, 3)), q, "monomer", [name], "sp", False,
                                {"kind": "cation", "index": 0}),
                           mults=[_hund_mult(n_electrons([name], q))])  # Li+ triplet is a 1s2s core excitation
                n += 1
                continue
            symbols, positions, mol = species.geometry(sp, ctx.workdir, charge=q)
            if q != sp.charge:
                s0, p0, _ = species.geometry(sp, ctx.workdir)
                comps = _components(symbols, positions)  # Li/Na are single-atom groups
                if len(comps) != len(_components(s0, p0)) or build.fragmented(positions, comps):
                    warnings.warn(f"monomer {name} q={q}: xtb opt dissociated the molecule; skipped")
                    continue
            yield _rec(ctx, symbols, positions, q, "monomer", [name], "opt", False, {"kind": "opt", "index": 0})
            n += 1
            if n >= _limit(scfg):
                return

            seed = int(ctx.cfg.get("seed", 0))
            rng = np.random.default_rng(seed + zlib.crc32(f"monomer{name}{q}".encode()))
            rng_d = np.random.default_rng(seed + zlib.crc32(f"monomer_distort{name}{q}".encode()))  # a failed hess must not shift these
            r = xtb.run_xtb(symbols, positions, q, _uhf(symbols, q), ctx.workdir, mode="hess",
                            timeout=_timeout(ctx))
            for k in range(int(scfg.get("n_normal_mode", 0))):
                if n >= _limit(scfg):
                    return
                pos = _normal_mode(symbols, positions, r.get("hessian"), float(scfg.get("T_K", 500)), rng) \
                    if r.get("ok") else None
                if pos is None:
                    break
                yield _rec(ctx, symbols, pos, q, "monomer", [name], "sp", False,
                           {"kind": "normal_mode", "T_K": float(scfg.get("T_K", 500)), "index": k} | _etemp(r))
                n += 1

            frames = sample(_with_positions(mol, positions), int(scfg.get("n_distort", 0)), rng=rng_d)
            for k, pos in enumerate(frames[1:int(scfg.get("n_distort", 0)) + 1]):
                if n >= _limit(scfg):
                    return
                yield _rec(ctx, symbols, pos, q, "monomer", [name], "sp", False,
                           {"kind": "distort", "sample": k})
                n += 1


MADELUNG_L = (0, 0, 1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1)  # 1s 2s 2p 3s 3p 4s 3d 4p ... 7p


def _hund_mult(n):
    """Ground-state multiplicity of an atom or atomic ion with n electrons: Hund's rule on the last open subshell
    in aufbau order. Right for main-group atoms (H-Ar, K, Ca, Ga-Kr, ...); None for an open d or f shell, where
    aufbau exceptions (Cr, Cu) and the 4s/3d order of ions make the rule unreliable."""
    for l in MADELUNG_L:
        cap = 2 * (2 * l + 1)
        if n <= cap:
            return None if l >= 2 and 0 < n < cap else min(n, cap - n) + 1
        n -= cap
    return None


def atoms(ctx, scfg, delta_q):
    """Isolated-atom reference decks (delta_q 0 only): every element in the configured species, q = -1..+1."""
    if int(delta_q) != 0:
        return
    els = set()
    for kind in ("cations", "anions", "solvents"):
        for name in ctx.cfg["species"].get(kind, []):
            mol = Chem.MolFromSmiles(ctx.cat[name].smiles)
            if mol is not None:
                els.update(a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms())
    lo, hi = ctx.cfg.get("charge", {}).get("range", (-2, 2))
    n = 0
    for el in sorted(els):
        for c in (-1, 0, 1):
            if not lo <= c <= hi or n_electrons([el], c) == 0:
                continue  # H+ has no electrons
            if n >= _limit(scfg):
                return
            mult = _hund_mult(n_electrons([el], c))
            if mult is None:
                warnings.warn(f"atoms: {el} q={c} has an open d/f shell; ground-state multiplicity not set, skipped")
                continue
            yield dict(symbols=[el], positions=np.zeros((1, 3)), charge=c, sampler="atoms",
                       composition=el, job="sp", bs=False, params={"kind": "atom", "element": el}, mults=[mult])
            n += 1


def _normal_mode(symbols, positions, hessian, T_K, rng, max_move=0.5, tries=100):
    """Classical harmonic sample: q_k ~ N(0, kT/lambda_k), displacement M^-1/2 sum q_k L_k (bohr -> A)."""
    if hessian is None:
        return None
    n = len(symbols)
    m3 = np.repeat([atomic_masses[atomic_numbers[s]] for s in symbols], 3)
    hm = np.asarray(hessian, float) / np.sqrt(np.outer(m3, m3))
    lam, L = np.linalg.eigh(hm)
    pos = np.asarray(positions, float)
    centred = pos - pos.mean(axis=0)
    linear = np.linalg.matrix_rank(centred, tol=1e-3) <= 1
    keep = slice(5 if linear else 6, None)
    lam, L = lam[keep], L[:, keep]
    freq = CM_PER_SQRT_EH_BOHR2_AMU * np.sqrt(np.clip(lam, 0.0, None))
    sel = (lam > 0) & (freq >= 50.0)
    lam, L = lam[sel], L[:, sel]
    if not len(lam):
        return None
    kT = 3.166811563e-6 * float(T_K)
    scale = np.sqrt(kT / lam)
    for _ in range(tries):
        q = rng.normal(size=len(lam)) * scale
        disp = (L @ q / np.sqrt(m3)).reshape(n, 3) * BOHR_A
        if np.linalg.norm(disp, axis=1).max() <= max_move:
            return pos + disp
    return None


def pair_scan(ctx, scfg, delta_q):
    """Contact pair assembly, then rigid translation of member B along the COM axis at each scale."""
    names = [n for k in ("cations", "anions", "solvents") for n in ctx.cfg["species"].get(k, [])]
    skip = ctx.cfg.get("_skip_pairs", ())  # cluster.skip_done
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i:] if (a, b) not in skip]
    max_pairs = int(scfg.get("max_pairs", 10 ** 9))
    if len(pairs) > max_pairs:
        pairs = [pairs[i] for i in ctx.rng.permutation(len(pairs))[:max_pairs]]
    n = 0
    for pair in pairs:
        getattr(ctx, "progress", lambda *a: None)("+".join(pair))
        if _cation_only(ctx, pair) and int(delta_q) > 0:
            continue  # no extra oxidation of bare cation pairs
        q = _formal(ctx, pair) + int(delta_q)
        rng = _stream(ctx, "pair_scan", pair, q)
        for orient in range(1 if _cation_only(ctx, pair) else int(scfg.get("n_orient", 1))):  # bare ions: no orientation
            if n >= _limit(scfg):
                return
            if not _q_ok(ctx, q, pair, delta_q):
                break
            extra = {}
            if _cation_only(ctx, pair):
                sa = species.geometry(ctx.cat[pair[0]], ctx.workdir)[0]
                sb = species.geometry(ctx.cat[pair[1]], ctx.workdir)[0]
                dc = _cov(sa, 0) + _cov(sb, 0)  # cation radii sum is the contact assumption (Li-Li 2.56 A)
                pos = np.zeros((len(sa) + len(sb), 3))
                pos[len(sa):, 0] = dc  # B on +x, A at the origin
                cl = {"symbols": sa + sb, "positions": pos,
                      "members": [(pair[0], list(range(len(sa)))), (pair[1], list(range(len(sa), len(sa) + len(sb))))]}
                (_, ia), (_, ib) = cl["members"]
                d0 = dc
                u = np.array([1.0, 0.0, 0.0])
                r = None
                extra["contact"] = "rcov_sum"
            else:
                cl = build.assemble(list(pair), ctx.cat, ctx.workdir, q, rng)
                (_, ia), (_, ib) = cl["members"]
                pos = cl["positions"]
                d0 = float(np.linalg.norm(_com(pos, ib) - _com(pos, ia)))
                if d0 < 1e-6:
                    break
                u = (_com(pos, ib) - _com(pos, ia)) / d0
                r = None
                if _min_between(pos, ia, ib) > 3.0:  # gas-phase xTB opt left the pair apart: bring into contact
                    pos = _shift_to(pos, ia, ib, u, 2.5, d0)
                    r = xtb.run_xtb(cl["symbols"], pos, q, _uhf(cl["symbols"], q), ctx.workdir, mode="opt",
                                    timeout=_timeout(ctx))
                    new = np.asarray(r["positions"], float) if r.get("ok") and "positions" in r else None
                    if new is not None and _min_between(new, ia, ib) <= 3.0:
                        pos = new  # an opt that pushes the pair apart again (e.g. two anions) keeps the 2.5 A contact
                    else:
                        r = None
                    d0 = float(np.linalg.norm(_com(pos, ib) - _com(pos, ia)))
                    u = (_com(pos, ib) - _com(pos, ia)) / d0
            dc = _min_between(pos, ia, ib)  # contact = shortest inter-member distance
            if not extra:  # an xtb-relaxed adduct: a rigid distance scan would stretch its new bond
                sa, sb = set(ia), set(ib)
                if any((i in sa and j in sb) or (j in sa and i in sb) for i, j in _bonds(cl["symbols"], pos)):
                    warnings.warn(f"pair_scan {list(pair)} q={q} orient {orient}: xtb opt bonded the members "
                                  f"(d_min {dc:.2f} A); orientation skipped")
                    continue
            for s in scfg.get("scales", [1.0]):
                if n >= _limit(scfg):
                    return
                moved = _shift_to(pos, ia, ib, u, s * dc, d0)
                yield _rec(ctx, cl["symbols"], moved, q, "pair_scan", pair, "sp", _bs(ctx, pair, delta_q),
                           {"pair": list(pair), "orient": orient, "scale": float(s), "contact_A": dc,
                            "d_min_A": _min_between(moved, ia, ib)} | extra | _etemp(cl, r))
                n += 1


def cluster(ctx, scfg, delta_q):
    """Assembled cluster opt + Gaussian-rattled sp frames + rigid member dilations."""
    n = 0
    for combo in _combos(ctx):
        if len(combo) < 2 or n >= _limit(scfg):
            continue
        getattr(ctx, "progress", lambda *a: None)(build.composition(list(combo)))  # plan counts >= 2 members
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
        rng = _stream(ctx, "cluster", combo, q)
        yield _rec(ctx, cl["symbols"], cl["positions"], q, "cluster", combo, "opt", _bs(ctx, combo, delta_q),
                   {"combo": list(combo), "kind": "opt"} | _etemp(cl))
        n += 1
        if n >= _limit(scfg):
            return
        rattle = scfg.get("rattle_A", 0.05)
        sigmas = rattle if isinstance(rattle, (list, tuple)) else [rattle]
        k = 0
        for sigma in sigmas:
            for _ in range(int(scfg.get("n_distort", 0))):
                if n >= _limit(scfg):
                    return
                pos = _rattle(cl["positions"], float(sigma), rng)
                if pos is None:
                    continue
                yield _rec(ctx, cl["symbols"], pos, q, "cluster", combo, "sp", _bs(ctx, combo, delta_q),
                           {"combo": list(combo), "kind": "distort", "index": k,
                            "sigma_A": float(sigma)} | _etemp(cl))
                n += 1
                k += 1
        if len(cl["members"]) < 2:
            continue
        member_of = np.empty(len(cl["symbols"]), int)
        for m, (_, idxs) in enumerate(cl["members"]):
            member_of[idxs] = m
        other = member_of[:, None] != member_of[None]
        com = cl["positions"].mean(axis=0)
        for f in scfg.get("dilate", [0.9, 1.1, 1.25]):
            if n >= _limit(scfg):
                return
            moved = cl["positions"].copy()
            for _, idxs in cl["members"]:
                moved[idxs] += (float(f) - 1.0) * (_com(cl["positions"], idxs) - com)
            if _too_close(cl["symbols"], moved, 0.7, other):
                continue  # rigid dilation pushed two members into each other
            yield _rec(ctx, cl["symbols"], moved, q, "cluster", combo, "sp", _bs(ctx, combo, delta_q),
                       {"combo": list(combo), "kind": "dilate", "factor": float(f)} | _etemp(cl))
            n += 1
        cats = [(name, idxs) for name, idxs in cl["members"] if ctx.cat[name].kind == "cation"]
        if len(cats) != 1:
            continue
        cat_pos = _com(cl["positions"], cats[0][1])
        for f in scfg.get("m_stretch", [1.05, 1.10]):
            if n >= _limit(scfg):
                return
            moved = cl["positions"].copy()
            for name, idxs in cl["members"]:
                if ctx.cat[name].kind == "cation":
                    continue
                u = _com(cl["positions"], idxs) - cat_pos
                d = float(np.linalg.norm(u))
                if d < 1e-8:
                    continue
                # GFN2 puts M-O ~0.1 A (Li) to ~0.3 A (Na) shorter than DFT: stretch the ionic contact
                moved[idxs] += (float(f) - 1.0) * _min_between(cl["positions"], cats[0][1], idxs) * u / d
            if _too_close(cl["symbols"], moved, 0.7):
                continue
            yield _rec(ctx, cl["symbols"], moved, q, "cluster", combo, "sp", _bs(ctx, combo, delta_q),
                       {"combo": list(combo), "kind": "m_stretch", "factor": float(f)} | _etemp(cl))
            n += 1


def _rattle(positions, sigma, rng, min_dist=0.7, tries=100):
    n = len(positions)
    for _ in range(tries):
        pos = positions + rng.normal(0.0, sigma, size=positions.shape)
        d = np.linalg.norm(pos[:, None] - pos[None], axis=-1)
        d[np.diag_indices(n)] = np.inf
        if d.min() >= min_dist:
            return pos
    return None


def _rattled(x, sigma, rng, ok=None):
    """Rattled copy of a path frame when it passes ok, else None (the caller keeps the exact frame)."""
    if sigma <= 0:
        return None
    pos = _rattle(np.asarray(x, float), sigma, rng)
    return pos if pos is not None and (ok is None or ok(pos)) else None


def _bond_candidates(ctx, cl, combo, max_bonds, rng):
    """One intramolecular bond per symmetry class per non-cation member, at most max_bonds (random subset).
    Returns [(i, j, member atom indices)] in cluster indices."""
    seen, cand = set(), []
    for (name, idxs), (_, mol) in zip(cl["members"], _mols(ctx, combo)):
        if ctx.cat[name].kind == "cation":
            continue
        rank = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
        for bd in mol.GetBonds():
            i, j = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()
            key = (name, *sorted((rank[i], rank[j])))
            if key not in seen:
                seen.add(key)
                cand.append((idxs[i], idxs[j], mol, idxs))
    if len(cand) > max_bonds:
        cand = [cand[i] for i in rng.permutation(len(cand))[:max_bonds]]
    return cand


def bond_scan1(ctx, scfg, delta_q):
    """Relaxed xtb scans of one intramolecular bond per symmetry class per member species."""
    bs = bool(scfg.get("bs", False))
    npts = int(scfg.get("n_points", 12))
    max_bonds = int(scfg.get("max_bonds_per_cluster", 10 ** 9))
    sigma = float(scfg.get("path_rattle_A", 0.07))
    n = 0
    for combo in _combos(ctx):
        getattr(ctx, "progress", lambda *a: None)(build.composition(list(combo)))
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
        rng = _stream(ctx, "bond_scan1", combo, q)
        symbols, pos = cl["symbols"], cl["positions"]
        for i, j, _, _ in _bond_candidates(ctx, cl, combo, max_bonds, rng):
            prng = np.random.default_rng(int(rng.integers(2 ** 63)))  # drawn before xtb: a failure must not shift later noise
            r0 = float(np.linalg.norm(pos[i] - pos[j]))
            r1 = 1.05 * (_cov(symbols, i) + _cov(symbols, j)) + 1.5
            r = xtb.run_xtb(symbols, pos, q, _uhf(symbols, q), ctx.workdir, mode="sp",
                            constraints=[(i, j, r0)], scan=dict(n=npts, ranges=[(0.85 * r0, r1)]),
                            timeout=_timeout(ctx))
            if not r.get("ok"):
                continue
            kmax = min(_until_jump(r["E_eV"]), _until_apart(cl, r["frames"]))
            for k, (frame, energy) in enumerate(zip(r["frames"][:kmax], r["E_eV"][:kmax])):
                if n >= _limit(scfg):
                    return
                frame = np.asarray(frame, float)
                rr = _rattled(frame, sigma, prng) if k % 2 else None
                yield _rec(ctx, symbols, frame if rr is None else rr, q, "bond_scan1", combo, "sp", bs,
                           {"bond": [i, j], "elements": [symbols[i], symbols[j]],
                            "r_A": float(np.linalg.norm(frame[i] - frame[j])), "index": k,
                            "E_xtb_eV": float(energy)}
                           | ({"rattle_A": sigma} if rr is not None else {}) | _etemp(cl, r))
                n += 1


def _idpp(images):
    """IDPP interpolation (ase.mep) of images in place; ASE's NEB default-method notice is not ours to report."""
    from ase.mep import NEB
    from ase.mep.neb import idpp_interpolate
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="The default method has changed")
        NEB(images, method="improvedtangent").interpolate()
        idpp_interpolate(images, traj=None, log=None)


def _side(mol, i, j):
    """Local atom indices on j's side of bond i-j; just [j] when the bond is in a ring."""
    side, todo = {j}, [j]
    while todo:
        a = todo.pop()
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            k = nb.GetIdx()
            if (a, k) != (j, i) and k not in side:
                side.add(k)
                todo.append(k)
    return [j] if i in side else sorted(side)


def bond_scan2(ctx, scfg, delta_q):
    """Unrelaxed bond-breaking paths: IDPP interpolation (ase.mep) from the cluster to one xtb endpoint.

    Endpoint: j's side of the bond (only atom j for a ring bond) is moved along the bond to
    r1 = 1.05 x (sum of covalent radii) + 1.5 A, then optimised by xtb with that distance constrained. The
    n_points images between the cluster and the endpoint get no energy model, as in the OMol25 geodesic /
    Popcornn paths, so they carry no xtb scan artefacts (state switches, rearrangements along the scan)."""
    bs = bool(scfg.get("bs", False))
    npts = int(scfg.get("n_points", 12))
    max_bonds = int(scfg.get("max_bonds_per_cluster", 10 ** 9))
    sigma = float(scfg.get("path_rattle_A", 0.07))
    n = 0
    for combo in _combos(ctx):
        getattr(ctx, "progress", lambda *a: None)(build.composition(list(combo)))
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
        rng = _stream(ctx, "bond_scan2", combo, q)
        symbols, pos = cl["symbols"], np.asarray(cl["positions"], float)
        rcov = np.array([_cov(symbols, k) for k in range(len(symbols))])
        bonds = _bonds(symbols, pos)
        for i, j, mol, idxs in _bond_candidates(ctx, cl, combo, max_bonds, rng):
            prng = np.random.default_rng(int(rng.integers(2 ** 63)))  # drawn before xtb: a failure must not shift later noise
            li, lj = idxs.index(i), idxs.index(j)
            key = (min(i, j), max(i, j))
            side = [idxs[k] for k in _side(mol, li, lj)]
            r0 = float(np.linalg.norm(pos[i] - pos[j]))
            r1 = 1.05 * (rcov[i] + rcov[j]) + 1.5
            start = pos.copy()
            start[side] += (r1 - r0) * (pos[j] - pos[i]) / r0
            r = xtb.run_xtb(symbols, start, q, _uhf(symbols, q), ctx.workdir, mode="opt",
                            constraints=[(i, j, r1)], timeout=_timeout(ctx))
            if not (r.get("ok") and "positions" in r):
                continue
            end, endpoint = np.asarray(r["positions"], float), "xtb"
            be = _bonds(symbols, end)
            unexpected = (be - bonds) | (bonds - be - {key})
            if unexpected:  # the endpoint rearranged beyond breaking i-j
                if mol.GetBondBetweenAtoms(li, lj).IsInRing():
                    warnings.warn(f"bond_scan2 {list(combo)} q={q} bond {symbols[i]}{i}-{symbols[j]}{j}: "
                                  f"xtb endpoint rearranged; skipped")
                    continue
                if _bonds(symbols, start) != bonds - {key}:  # moving the fragment hit another atom
                    warnings.warn(f"bond_scan2 {list(combo)} q={q} bond {symbols[i]}{i}-{symbols[j]}{j}: "
                                  f"xtb endpoint rearranged and the rigid one makes new bonds; skipped")
                    continue
                end, endpoint = start, "rigid"  # keep the unrelaxed displaced geometry as the endpoint
            if build.fragmented(end, [ix for _, ix in cl["members"]]):
                continue  # endpoint opt pushed whole members apart: not a bond-breaking path
            images = [Atoms(symbols, positions=pos)] + [Atoms(symbols, positions=pos) for _ in range(npts - 2)] \
                + [Atoms(symbols, positions=end)]
            _idpp(images)

            def ok_at(y, key=key):
                if _too_close(symbols, y, 0.7):
                    return False  # two atoms pushed together
                return not any(np.linalg.norm(y[a] - y[b]) > 1.3 * np.linalg.norm(pos[a] - pos[b])
                               for a, b in bonds if (a, b) != key)  # another bond torn

            def ok_far(y):
                return ok_at(y) and not (_bonds(symbols, y) - bonds)  # no new bond to a spectator

            e = 0  # emission index within this path: odd frames carry thermal noise
            for k, im in enumerate(images[1:], 1):  # image 0 is the cluster itself (cluster opt deck)
                if n >= _limit(scfg):
                    return
                x = im.get_positions()
                if not ok_at(x):
                    continue
                rr = _rattled(x, sigma, prng, ok_at) if e % 2 else None
                yield _rec(ctx, symbols, x if rr is None else rr, q, "bond_scan2", combo, "sp", bs,
                           {"bond": [i, j], "elements": [symbols[i], symbols[j]], "method": "idpp",
                            "r_A": float(np.linalg.norm(x[i] - x[j])), "index": k,
                            "ring": bool(mol.GetBondBetweenAtoms(li, lj).IsInRing()),
                            "endpoint": endpoint}
                           | ({"rattle_A": sigma} if rr is not None else {})
                           | ({"E_xtb_end_eV": float(r["E_eV"])} if endpoint == "xtb" else {}) | _etemp(cl, r))
                e += 1
                n += 1
            if not mol.GetBondBetweenAtoms(li, lj).IsInRing():  # far points off a ring bond are meaningless
                vec = end[j] - end[i]
                d_end = float(np.linalg.norm(vec))
                if d_end < 1e-8:
                    continue
                for k, dist in enumerate([float(v) for v in scfg.get("far_A", [4.0, 5.0]) if float(v) > d_end],
                                         len(images)):
                    if n >= _limit(scfg):
                        return
                    x = end.copy()
                    x[side] += (dist - d_end) * vec / d_end  # j-side fragment moved rigidly along i->j
                    if not ok_far(x):
                        continue
                    rr = _rattled(x, sigma, prng, ok_far) if e % 2 else None
                    yield _rec(ctx, symbols, x if rr is None else rr, q, "bond_scan2", combo, "sp", bs,
                               {"bond": [i, j], "elements": [symbols[i], symbols[j]], "method": "far",
                                "r_A": float(np.linalg.norm(x[i] - x[j])), "index": k,
                                "ring": False, "endpoint": endpoint, "far": True}
                               | ({"rattle_A": sigma} if rr is not None else {})
                               | ({"E_xtb_end_eV": float(r["E_eV"])} if endpoint == "xtb" else {}) | _etemp(cl, r))
                    e += 1
                    n += 1


def _electrophiles(mol):
    """Carbonyl C (C=O) and P/S/B atoms."""
    out = []
    for a in mol.GetAtoms():
        i = a.GetIdx()
        if a.GetSymbol() in ELECTROPHILES:
            out.append(i)
        elif a.GetSymbol() == "C" and any(
                b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetSymbol() == "O"
                for b in a.GetBonds()):
            out.append(i)
    return out


def _carbonyl_oxygens(mol, c):
    return {b.GetOtherAtom(mol.GetAtomWithIdx(c)).GetIdx() for b in mol.GetAtomWithIdx(c).GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(mol.GetAtomWithIdx(c)).GetSymbol() == "O"}


def _reaction_candidates(ctx, cl, combo, extra_radicals=(), oxidised=False):
    """[(template, pairs, ranges)] over member pairs within 6 A; pairs[k] carries ranges[k] = (start, target).
    extra_radicals: cluster atoms to treat as radical centres (charge-derived, see _charge_radicals)."""
    infos, cat_infos = [], []
    xr = set(extra_radicals)
    for (name, idxs), (_, mol) in zip(cl["members"], _mols(ctx, combo)):
        info = dict(name=name, idxs=idxs, mol=mol,
                    rad=list(dict.fromkeys([idxs[i] for i in _radical_atoms(mol)]
                                           + [k for k in idxs if k in xr])))
        (cat_infos if ctx.cat[name].kind == "cation" else infos).append(info)
    pos = cl["positions"]
    rcov = np.array([_cov(cl["symbols"], k) for k in range(len(cl["symbols"]))])
    cand = []
    # ring_open_M first (cation-assisted): opening of a ring C-O single bond
    for ci in cat_infos:
        for a in infos:
            if np.linalg.norm(_com(pos, ci["idxs"]) - _com(pos, a["idxs"])) > 6.0:
                continue
            rings = []
            for bd in a["mol"].GetBonds():
                if not bd.IsInRing():
                    continue
                c_atom, o_atom = bd.GetBeginAtom(), bd.GetEndAtom()
                if o_atom.GetSymbol() == "C":
                    c_atom, o_atom = o_atom, c_atom
                if c_atom.GetSymbol() != "C" or o_atom.GetSymbol() != "O":
                    continue
                if any(b.GetBondType() == Chem.BondType.DOUBLE for b in o_atom.GetBonds()):
                    continue
                rings.append((o_atom.GetIdx(), c_atom.GetIdx()))
            for o_local, c_local in sorted(rings):
                c, o = a["idxs"][c_local], a["idxs"][o_local]
                for m in ci["idxs"]:
                    cand.append(("ring_open_M", [(c, o), (m, o)],
                                 [(float(np.linalg.norm(pos[c] - pos[o])),
                                   1.05 * (rcov[c] + rcov[o]) + 1.5),
                                  (float(np.linalg.norm(pos[m] - pos[o])),
                                   build.D_M.get(ci["name"], 2.0))]))
    for ia, a in enumerate(infos):
        donors = [a["idxs"][i] for i in range(a["mol"].GetNumAtoms())
                  if a["mol"].GetAtomWithIdx(i).GetSymbol() in DONORS]
        for b in infos:
            if b is a:
                continue
            if np.linalg.norm(_com(pos, a["idxs"]) - _com(pos, b["idxs"])) > 6.0:
                continue
            for e_local in _electrophiles(b["mol"]):
                e = b["idxs"][e_local]
                for nu in donors:
                    d0 = float(np.linalg.norm(pos[nu] - pos[e]))
                    tgt = 1.05 * (rcov[nu] + rcov[e])
                    cand.append(("add", [(nu, e)], [(d0, tgt)]))
                    carbonyl_o = _carbonyl_oxygens(b["mol"], e_local)
                    for x in b["mol"].GetAtomWithIdx(e_local).GetNeighbors():
                        if x.GetSymbol() not in ("O", "F") or x.GetIdx() in carbonyl_o:
                            continue
                        l = b["idxs"][x.GetIdx()]
                        r0 = float(np.linalg.norm(pos[e] - pos[l]))
                        tgt2 = 1.05 * (rcov[e] + rcov[l]) + 1.5
                        cand.append(("substitute", [(nu, e), (e, l)], [(d0, tgt), (r0, tgt2)]))
    # transfer: X (F/H) bonded to Y in a non-cation member A -> cation, or an electrophile of another member
    cations = [k for ci in cat_infos for k in ci["idxs"]]
    for a in infos:
        for x_local in range(a["mol"].GetNumAtoms()):
            if a["mol"].GetAtomWithIdx(x_local).GetSymbol() not in ("F", "H"):
                continue
            for y_atom in a["mol"].GetAtomWithIdx(x_local).GetNeighbors():
                x, y = a["idxs"][x_local], a["idxs"][y_atom.GetIdx()]
                rxy = float(np.linalg.norm(pos[y] - pos[x]))
                tgt = 1.05 * (rcov[y] + rcov[x]) + 1.2
                for z in cations:
                    if np.linalg.norm(_com(pos, a["idxs"]) - pos[z]) > 6.0:
                        continue
                    dz = float(np.linalg.norm(pos[z] - pos[x]))
                    ztgt = 0.9 * (rcov[z] + rcov[x])
                    cand.append(("transfer", [(y, x), (z, x)], [(rxy, tgt), (dz, ztgt)]))
                for b in infos:
                    if b is a or np.linalg.norm(_com(pos, a["idxs"]) - _com(pos, b["idxs"])) > 6.0:
                        continue
                    for e_local in _electrophiles(b["mol"]):
                        z = b["idxs"][e_local]
                        dz = float(np.linalg.norm(pos[z] - pos[x]))
                        ztgt = 1.0 * (rcov[z] + rcov[x])
                        cand.append(("transfer", [(y, x), (z, x)], [(rxy, tgt), (dz, ztgt)]))
    # h_transfer_ox (oxidised clusters only): an ethylene H of A moves to a carbonyl O of a neutral member
    # or an F/O/N of an anion member. Borodin & Xing, PCCP 14, 12838 (2012): in EC radical cation dimers the
    # ethylene H is abstracted by the carbonyl O of a neighbouring EC (exothermic). Xing et al., J. Phys.
    # Chem. A 115, 13896 (2011): oxidised PC-PF6 forms HF + PF5 in the first step.
    if oxidised:
        for a in infos:
            hs = [(a["idxs"][h.GetIdx()], a["idxs"][c.GetIdx()])
                  for c in a["mol"].GetAtoms() if c.GetSymbol() == "C"
                  for h in c.GetNeighbors() if h.GetSymbol() == "H"]
            for h, c in hs:
                rch = float(np.linalg.norm(pos[c] - pos[h]))
                for b in infos:
                    if b is a or np.linalg.norm(_com(pos, a["idxs"]) - _com(pos, b["idxs"])) > 6.0:
                        continue
                    sp_b = ctx.cat[b["name"]]
                    if sp_b.kind == "anion":
                        accept = [k for k, at in enumerate(b["mol"].GetAtoms())
                                  if at.GetSymbol() in ("F", "O", "N")]
                    elif int(sp_b.charge) == 0:
                        accept = [k for k, at in enumerate(b["mol"].GetAtoms())
                                  if at.GetSymbol() == "O" and any(
                                      bd.GetBondType() == Chem.BondType.DOUBLE for bd in at.GetBonds())]
                    else:
                        continue
                    for k in accept:
                        y = b["idxs"][k]
                        dyh = float(np.linalg.norm(pos[y] - pos[h]))
                        cand.append(("h_transfer_ox", [(c, h), (y, h)],
                                     [(rch, 1.05 * (rcov[c] + rcov[h]) + 1.2),
                                      (dyh, 1.05 * (rcov[y] + rcov[h]))]))
    # radical_couple: two radical centres of different members
    for ai, a in enumerate(infos):
        ra = a["rad"]
        if not ra:
            continue
        for b in infos[ai + 1:]:
            if np.linalg.norm(_com(pos, a["idxs"]) - _com(pos, b["idxs"])) > 6.0:
                continue
            for i in ra:
                for j in b["rad"]:
                    d0 = float(np.linalg.norm(pos[i] - pos[j]))
                    cand.append(("radical_couple", [(i, j)], [(d0, 1.05 * (rcov[i] + rcov[j]))]))
    # h_abstract: a radical centre of A takes an H bonded to Y of a different member B
    for a in infos:
        ra = a["rad"]
        if not ra:
            continue
        for b in infos:
            if b is a or np.linalg.norm(_com(pos, a["idxs"]) - _com(pos, b["idxs"])) > 6.0:
                continue
            for h_local in range(b["mol"].GetNumAtoms()):
                if b["mol"].GetAtomWithIdx(h_local).GetSymbol() != "H":
                    continue
                h = b["idxs"][h_local]
                for y_atom in b["mol"].GetAtomWithIdx(h_local).GetNeighbors():
                    y = b["idxs"][y_atom.GetIdx()]
                    r0 = float(np.linalg.norm(pos[y] - pos[h]))
                    tgt0 = 1.05 * (rcov[y] + rcov[h]) + 1.2
                    for r in ra:
                        d0 = float(np.linalg.norm(pos[r] - pos[h]))
                        cand.append(("h_abstract", [(y, h), (r, h)],
                                     [(r0, tgt0), (d0, 1.05 * (rcov[r] + rcov[h]))]))
    # intra_shift: H/F of one member moves to an atom 2-3 bonds away in the same member
    for a in infos:
        dm = Chem.GetDistanceMatrix(a["mol"])
        for x_local in range(a["mol"].GetNumAtoms()):
            xa = a["mol"].GetAtomWithIdx(x_local)
            if xa.GetSymbol() not in ("H", "F"):
                continue
            for y_atom in xa.GetNeighbors():
                y_local = y_atom.GetIdx()
                x, y = a["idxs"][x_local], a["idxs"][y_local]
                rxy = float(np.linalg.norm(pos[y] - pos[x]))
                tgt0 = 1.05 * (rcov[y] + rcov[x]) + 1.2
                for z_local in range(a["mol"].GetNumAtoms()):
                    if z_local == y_local or a["mol"].GetAtomWithIdx(z_local).GetSymbol() == "H":
                        continue
                    if int(dm[x_local, z_local]) not in (2, 3):
                        continue
                    z = a["idxs"][z_local]
                    dz = float(np.linalg.norm(pos[z] - pos[x]))
                    cand.append(("intra_shift", [(y, x), (z, x)],
                                 [(rxy, tgt0), (dz, 1.05 * (rcov[z] + rcov[x]))]))
    return cand


def _pick(cand, k, rng):
    """At most k candidates, taken round-robin over templates (each template shuffled), so a template with many
    candidates (intra_shift: every H to every atom 2-3 bonds away) does not crowd out the others."""
    groups = {}
    for c in cand:
        groups.setdefault(c[0], []).append(c)
    groups = [[g[i] for i in rng.permutation(len(g))] for _, g in sorted(groups.items())]
    groups = [groups[i] for i in rng.permutation(len(groups))]  # no template always first
    out = []
    while len(out) < k and any(groups):
        for g in groups:
            if g and len(out) < k:
                out.append(g.pop(0))
    return out


def _charge_radicals(ctx, cl, q, delta_q):
    """Cluster atoms whose xtb charge shifted most between q and q - delta_q, when the q state is a radical.

    Reduction (delta_q < 0) picks atoms that gained electron density (delta <= -0.15), oxidation (delta_q > 0)
    atoms that lost it (delta >= +0.15); heavy atoms of non-cation members only, up to 2, most |delta| first.
    The SMILES of a delta_q != 0 cluster has no radical marked, so these are its radical sites."""
    dq = int(delta_q)
    if dq == 0 or n_electrons(cl["symbols"], q) % 2 == 0:
        return []
    q0 = int(q) - dq
    r1 = xtb.run_xtb(cl["symbols"], cl["positions"], q, _uhf(cl["symbols"], q), ctx.workdir, mode="sp",
                     timeout=_timeout(ctx))
    r0 = xtb.run_xtb(cl["symbols"], cl["positions"], q0, _uhf(cl["symbols"], q0), ctx.workdir, mode="sp",
                     timeout=_timeout(ctx))
    if "charges" not in r1 or "charges" not in r0:
        return []  # older cached run or failed sp: radical sites unknown
    d = np.asarray(r1["charges"], float) - np.asarray(r0["charges"], float)
    cand = [i for name, idxs in cl["members"] if ctx.cat[name].kind != "cation"
            for i in idxs if cl["symbols"][i] not in ("H", "Li", "Na")  # Li/Na: ionic, not a radical centre
            and (d[i] <= -0.15 if dq < 0 else d[i] >= 0.15)]
    return sorted(cand, key=lambda i: -abs(d[i]))[:2]


def reaction(ctx, scfg, delta_q):
    """Pseudo-TS reaction paths (add / substitute / transfer): scan, keep the energy peak +/- k_peak."""
    bs = bool(scfg.get("bs", False))
    templates = set(scfg.get("templates", []))
    npts = int(scfg.get("n_points", 14))
    k_peak = int(scfg.get("k_peak", 2))
    max_per = int(scfg.get("max_per_cluster", 10 ** 9))
    sigma = float(scfg.get("path_rattle_A", 0.07))
    peak_n = int(scfg.get("peak_rattle_n", 2))
    n = 0
    for combo in _combos(ctx):
        getattr(ctx, "progress", lambda *a: None)(build.composition(list(combo)))
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
        rng = _stream(ctx, "reaction", combo, q)
        extra = _charge_radicals(ctx, cl, q, delta_q)
        cand = [c for c in _reaction_candidates(ctx, cl, combo, extra, int(delta_q) > 0) if c[0] in templates]
        cand = _pick(cand, max_per, rng)
        for template, pairs, ranges in cand:
            prng = np.random.default_rng(int(rng.integers(2 ** 63)))  # drawn before xtb: a failure must not shift later noise
            if n >= _limit(scfg):
                return
            atoms = list(dict.fromkeys(i for p in pairs for i in p))
            cons = [(a, b, r0) for (a, b), (r0, _) in zip(pairs, ranges)]
            r = xtb.run_xtb(cl["symbols"], cl["positions"], q, _uhf(cl["symbols"], q), ctx.workdir, mode="sp",
                            constraints=cons, scan=dict(n=npts, ranges=ranges), timeout=_timeout(ctx))
            if not r.get("ok"):
                continue
            E = r["E_eV"][:min(_until_jump(r["E_eV"]), _until_apart(cl, r["frames"]))]
            if len(E) < 3:
                continue
            p = 1 + int(np.argmax(E[1:-1]))
            keep = sorted({len(E) - 1} | set(range(max(1, p - k_peak), min(p + k_peak, len(E) - 2) + 1)))
            e = 0  # emission index within this path: odd frames carry thermal noise
            for k in keep:
                if n >= _limit(scfg):
                    return
                x = np.asarray(r["frames"][k], float)
                rr = _rattled(x, sigma, prng) if e % 2 and k not in (p, len(E) - 1) else None  # exact peak and end
                yield _rec(ctx, cl["symbols"], x if rr is None else rr, q, "reaction", combo, "sp", bs,
                           {"template": template, "atoms": atoms, "pairs": [list(pr) for pr in pairs],
                            "ranges": [list(rg) for rg in ranges],
                            "index": int(k), "peak_index": int(p), "is_peak": bool(k == p),
                            "charge_radicals": list(extra), "E_xtb_eV": float(E[k])}
                           | ({"rattle_A": sigma} if rr is not None else {}) | _etemp(cl, r))
                e += 1
                n += 1
            peak = np.asarray(r["frames"][p], float)
            for j in range(peak_n if sigma > 0 else 0):
                if n >= _limit(scfg):
                    return
                rr = _rattle(peak, sigma, prng)
                if rr is None:
                    continue
                yield _rec(ctx, cl["symbols"], rr, q, "reaction", combo, "sp", bs,
                           {"template": template, "atoms": atoms, "pairs": [list(pr) for pr in pairs],
                            "ranges": [list(rg) for rg in ranges],
                            "index": int(p), "peak_index": int(p), "is_peak": False, "peak_rattle": j,
                            "rattle_A": sigma, "charge_radicals": list(extra), "E_xtb_eV": float(E[p])}
                           | _etemp(cl, r))
                n += 1


def _radical_atoms(mol):
    """Local atom indices with unpaired electrons (catalog radicals: EC_ro, MeOCO2, CH3, CH3O)."""
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetNumRadicalElectrons() > 0]


def _endpoint_check(symbols, b0, end, pairs, ranges, tol=0.3):
    """('', formed, broken) if the endpoint made exactly the template's bond changes, else (reason, ...).
    A pair whose target is below 1.3 x r_cov sum (no Li/Na) must form, one stretched by >= 1.2 A must break, every
    constrained distance (ionic M-X ones included) must be within tol A of its target, and the cluster must not
    already be in the product state."""
    rc = [_cov(symbols, k) for k in range(len(symbols))]
    formed, broken = set(), set()
    for (a, b), (_, t) in zip(pairs, ranges):
        pair = (min(a, b), max(a, b))
        if t < 1.3 * (rc[a] + rc[b]) and not {"Li", "Na"} & {symbols[a], symbols[b]}:
            formed.add(pair)
        if t >= 1.05 * (rc[a] + rc[b]) + 1.2:
            broken.add(pair)
    b1 = _bonds(symbols, end)
    if formed & b0 or not broken <= b0:
        return "already in the product state", formed, broken
    if any(abs(float(np.linalg.norm(end[a] - end[b])) - t) > tol for (a, b), (_, t) in zip(pairs, ranges)):
        return "missed a target distance", formed, broken
    if not formed <= b1 or broken & b1:
        return "did not make the template bonds", formed, broken
    if (b0 ^ b1) - formed - broken:
        return "rearranged", formed, broken
    return "", formed, broken


def reaction2(ctx, scfg, delta_q):
    """Pseudo-TS reaction paths between the cluster and an xTB-optimised product endpoint (IDPP interpolation)."""
    bs = bool(scfg.get("bs", False))
    templates = set(scfg.get("templates", []))
    npts = int(scfg.get("n_points", 14))
    max_per = int(scfg.get("max_per_cluster", 10 ** 9))
    free_opt = bool(scfg.get("free_opt", True))
    sigma = float(scfg.get("path_rattle_A", 0.07))
    n = 0
    for combo in _combos(ctx):
        getattr(ctx, "progress", lambda *a: None)(build.composition(list(combo)))
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
        rng = _stream(ctx, "reaction2", combo, q)
        symbols, pos = cl["symbols"], np.asarray(cl["positions"], float)
        rcov = np.array([_cov(symbols, k) for k in range(len(symbols))])
        members = [ix for _, ix in cl["members"]]
        b0 = _bonds(symbols, pos)
        extra = _charge_radicals(ctx, cl, q, delta_q)
        cand = [c for c in _reaction_candidates(ctx, cl, combo, extra, int(delta_q) > 0) if c[0] in templates]
        cand = _pick(cand, max_per, rng)
        for template, pairs, ranges in cand:
            prng = np.random.default_rng(int(rng.integers(2 ** 63)))  # drawn before xtb: a failure must not shift later noise
            if n >= _limit(scfg):
                return
            atoms = list(dict.fromkeys(i for p in pairs for i in p))
            cons = [(a, b, t) for (a, b), (_, t) in zip(pairs, ranges)]
            r = xtb.run_xtb(symbols, pos, q, _uhf(symbols, q), ctx.workdir, mode="opt",
                            constraints=cons, timeout=_timeout(ctx))
            if not (r.get("ok") and "positions" in r):
                continue
            end = np.asarray(r["positions"], float)
            if build.fragmented(end, members):
                continue
            why, formed, broken = _endpoint_check(symbols, b0, end, pairs, ranges)
            if why:
                warnings.warn(f"reaction2 {template} {list(combo)} q={q}: endpoint {why}; skipped")
                continue
            b1 = _bonds(symbols, end)
            endpoint = "constrained"
            if free_opt:
                rf = xtb.run_xtb(symbols, end, q, _uhf(symbols, q), ctx.workdir, mode="opt",
                                 timeout=_timeout(ctx))
                e2 = np.asarray(rf.get("positions"), float) if rf.get("ok") and "positions" in rf else None
                if e2 is not None and not build.fragmented(e2, members) and _bonds(symbols, e2) == b1:
                    end, endpoint, r = e2, "free", rf
            images = [Atoms(symbols, positions=pos)] + [Atoms(symbols, positions=pos) for _ in range(npts - 2)] \
                + [Atoms(symbols, positions=end)]
            _idpp(images)

            def ok_at(y):
                if _too_close(symbols, y, 0.7):
                    return False  # two atoms pushed together
                return not any(np.linalg.norm(y[i] - y[j]) > 1.3 * np.linalg.norm(pos[i] - pos[j])
                               for i, j in b0 - broken)  # a surviving bond stretched far beyond its cluster length

            e = 0  # emission index within this path: odd frames carry thermal noise
            for k, im in enumerate(images[1:], 1):  # image 0 is the cluster itself (cluster opt deck)
                if n >= _limit(scfg):
                    return
                x = im.get_positions()
                if not ok_at(x):
                    continue
                rr = _rattled(x, sigma, prng, ok_at) if e % 2 else None
                yield _rec(ctx, symbols, x if rr is None else rr, q, "reaction2", combo, "sp", bs,
                           {"template": template, "atoms": atoms, "pairs": [list(p) for p in pairs],
                            "index": k, "endpoint": endpoint, "charge_radicals": list(extra),
                            "E_xtb_end_eV": float(r["E_eV"])}
                           | ({"rattle_A": sigma} if rr is not None else {}) | _etemp(cl, r))
                e += 1
                n += 1


def hop(ctx, scfg, delta_q):
    """Cation hopping path P1 -> P2 (in front of the target donor). At each point the cation's distances to the
    start and target donors (and their mutual distance) are constrained and the rest is xTB-relaxed. xtb $fix was
    not exact for these clusters (fixed atoms moved up to 2 A), so distance constraints are used. Frames where the
    cation is closer than 0.75 x covalent-radius sum to any atom are dropped."""
    npts = int(scfg.get("n_points", 10))
    max_per = int(scfg.get("max_per_cluster", 10 ** 9))
    sigma = float(scfg.get("path_rattle_A", 0.07))
    n = 0
    for combo in _combos(ctx):
        getattr(ctx, "progress", lambda *a: None)(build.composition(list(combo)))
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
        rng = _stream(ctx, "hop", combo, q)
        symbols, pos = cl["symbols"], cl["positions"]
        mols = _mols(ctx, combo)
        for cation, cidx in cl["members"]:
            if ctx.cat[cation].kind != "cation":
                continue
            m = cidx[0]
            donors = [(k, name, idxs, i) for k, ((name, idxs), (_, mol))
                      in enumerate(zip(cl["members"], mols))
                      if ctx.cat[name].kind != "cation"
                      for i in idxs if symbols[i] in DONORS]
            if not donors:
                continue
            d1 = int(min(donors, key=lambda z: np.linalg.norm(pos[m] - pos[z[3]]))[3])
            # the start donor itself is a target when the cation sits > 2.6 A from it (e.g. Na0): exclude it
            targets = [z for z in donors if z[3] != d1 and 2.6 < np.linalg.norm(pos[m] - pos[z[3]]) < 6.0]
            if len(targets) > max_per:
                targets = [targets[i] for i in rng.permutation(len(targets))[:max_per]]
            for mem, name, idxs, d2 in targets:
                pseed = int(rng.integers(2 ** 63))  # drawn before xtb; noise per path point: failures shift nothing
                mol = mols[mem][1]
                local = idxs.index(d2)
                nb = [idxs[a.GetIdx()] for a in mol.GetAtomWithIdx(local).GetNeighbors()]
                u = pos[d2] - (pos[nb].mean(axis=0) if nb else _com(pos, idxs))
                if np.linalg.norm(u) < 1e-8:
                    continue
                p2 = pos[d2] + build.D_M.get(cation, 2.0) * u / np.linalg.norm(u)
                p1 = pos[m].copy()
                rest = [i for i in range(len(symbols)) if i != m]
                rcov = np.array([_cov(symbols, i) for i in range(len(symbols))])
                r3, best = None, -1.0  # third reference atom: pins the cation to the d1-d2 axis
                for k, s in enumerate(symbols):
                    if s == "H" or s in ("Li", "Na") or k in (d1, int(d2)):
                        continue
                    if np.linalg.norm(pos[k] - 0.5 * (p1 + p2)) > 6.0:
                        continue
                    v = pos[k] - pos[d1]
                    ab = pos[int(d2)] - pos[d1]
                    perp = np.linalg.norm(v - (v @ ab) / (ab @ ab) * ab)
                    if perp > best:  # ties: lowest index wins
                        r3, best = k, perp
                if best < 0.5:
                    r3 = None  # (nearly) on the d1-d2 axis: it would not pin the cation
                cur = pos.copy()

                def ok_at(y, m=m, rest=rest, rcov=rcov):
                    return not (np.linalg.norm(y[rest] - y[m], axis=1) < 0.75 * (rcov[m] + rcov[rest])).any()

                for kt, t in enumerate(np.linspace(0.0, 1.0, npts)[1:], 1):  # t=0 is the relaxed cluster (its opt deck)
                    if n >= _limit(scfg):
                        return
                    pt = (1 - t) * p1 + t * p2
                    cur = _kabsch(cur, pos, rest)  # back into the frame where p1, p2 are defined
                    cur[m] = pt
                    cons = [(m, d1, float(np.linalg.norm(pt - pos[d1]))), (m, int(d2), float(np.linalg.norm(pt - pos[d2]))),
                            (d1, int(d2), float(np.linalg.norm(pos[d1] - pos[d2])))]
                    if r3 is not None:
                        cons += [(m, r3, float(np.linalg.norm(pt - pos[r3]))),
                                 (d1, r3, float(np.linalg.norm(pos[d1] - pos[r3]))),
                                 (int(d2), r3, float(np.linalg.norm(pos[int(d2)] - pos[r3])))]
                    if any(r0 < 0.75 * (rcov[m] + rcov[d]) for _, d, r0 in cons[:2]):
                        continue  # path runs through a donor; the guard below would drop this frame anyway
                    r = xtb.run_xtb(symbols, cur, q, _uhf(symbols, q), ctx.workdir, mode="opt",
                                    constraints=cons, timeout=_timeout(ctx))
                    if not (r.get("ok") and "positions" in r):
                        continue
                    cur = np.asarray(r["positions"], float)
                    dev = float(np.linalg.norm(_kabsch(cur, pos, rest)[m] - pt))
                    if dev > HOP_MAX_DEV_A:
                        continue  # e.g. mirror image across the reference plane: not on the intended path
                    if not ok_at(cur):
                        continue  # cation runs into an atom: no useful PES point
                    rr = _rattled(cur, sigma, np.random.default_rng([pseed, kt]), ok_at) if kt % 2 == 0 else None
                    yield _rec(ctx, symbols, cur if rr is None else rr, q, "hop", combo, "sp", _bs(ctx, combo, delta_q),
                               {"kind": "hop", "cation": cation, "cation_index": int(m), "from_donor": d1,
                                "to_donor": int(d2), "t": float(t),
                                "d_from_A": float(np.linalg.norm(cur[m] - cur[d1])),
                                "d_to_A": float(np.linalg.norm(cur[m] - cur[d2])),
                                "path_dev_A": dev, "ref3": r3}
                               | ({"rattle_A": sigma} if rr is not None else {}) | _etemp(cl, r))
                    n += 1
        if bool(scfg.get("exchange", True)):
            emax = int(scfg.get("exchange_max_per_cluster", 2))
            ext = float(scfg.get("exchange_extend_A", 3.0))
            pairs = [(cidx[0], i) for c, cidx in cl["members"] if ctx.cat[c].kind == "cation"
                     for name, idxs in cl["members"] if ctx.cat[name].kind != "cation"
                     for i in idxs if symbols[i] in DONORS and np.linalg.norm(pos[cidx[0]] - pos[i]) <= 2.6]
            if len(pairs) > emax:
                pairs = [pairs[i] for i in rng.permutation(len(pairs))[:emax]]
            for m, d in pairs:
                prng = np.random.default_rng(int(rng.integers(2 ** 63)))  # drawn before xtb: a failure must not shift later noise
                if n >= _limit(scfg):
                    return
                r0 = float(np.linalg.norm(pos[m] - pos[d]))
                r = xtb.run_xtb(symbols, pos, q, _uhf(symbols, q), ctx.workdir, mode="sp",
                                constraints=[(m, d, r0)], scan=dict(n=npts, ranges=[(r0, r0 + ext)]),
                                timeout=_timeout(ctx))
                if not r.get("ok"):
                    continue
                kmax = min(_until_jump(r["E_eV"]), _until_apart(cl, r["frames"]))
                others = [k2 for k2 in range(len(symbols)) if k2 != m]
                rcov_others = 0.75 * (_cov(symbols, m) + np.array([_cov(symbols, k2) for k2 in others]))

                def ok_at(y, m=m, others=others, rcov_others=rcov_others):
                    return not (np.linalg.norm(y[others] - y[m], axis=1) < rcov_others).any()

                e = 0  # emission index within this path: odd frames carry thermal noise
                for k, (frame, energy) in enumerate(zip(r["frames"][:kmax], r["E_eV"][:kmax])):
                    if n >= _limit(scfg):
                        return
                    if k == 0:
                        continue  # k=0 is the relaxed cluster (its opt deck)
                    frame = np.asarray(frame, float)
                    if not ok_at(frame):
                        continue  # cation runs into an atom, as in the hop guard
                    rr = _rattled(frame, sigma, prng, ok_at) if e % 2 else None
                    yield _rec(ctx, symbols, frame if rr is None else rr, q, "hop", combo, "sp", _bs(ctx, combo, delta_q),
                               {"kind": "exchange_out", "cation_index": int(m), "donor": int(d),
                                "r_A": float(np.linalg.norm(frame[m] - frame[d])), "index": k,
                                "E_xtb_eV": float(energy)}
                               | ({"rattle_A": sigma} if rr is not None else {}) | _etemp(cl, r))
                    e += 1
                    n += 1


def _rmsd(x, ref):
    """All-atom RMSD (A) of x after Kabsch alignment onto ref."""
    idx = list(range(len(ref)))
    moved = _kabsch(np.asarray(x, float), np.asarray(ref, float), idx)
    return float(np.sqrt(((moved - np.asarray(ref, float)) ** 2).sum(axis=1).mean()))


def _maxmin(frames, k):
    """Indices of up to k frames picked by greedy max-min RMSD: start from frame 0, then repeatedly add the frame
    whose minimum RMSD to the already picked frames is largest."""
    picked = [0] if frames and k > 0 else []
    while len(picked) < min(int(k), len(frames)):
        best, bi = -1.0, None
        for j in range(len(frames)):
            if j in picked:
                continue
            d = min(_rmsd(frames[j], frames[i]) for i in picked)
            if d > best:
                best, bi = d, j
        picked.append(bi)
    return picked


def md_compositions(cfg, cat, scfg):
    """cluster_md compositions: cation + n solvents (+ one anion replacing a solvent), within max_atoms.
    n_solvent is a list (every cation) or a dict keyed by cation name; scfg max_atoms overrides limits.max_atoms."""
    sp = cfg["species"]
    max_atoms = int(scfg.get("max_atoms", cfg.get("limits", {}).get("max_atoms", 10 ** 9)))
    ns = scfg.get("n_solvent", [4])
    names_all = []
    for c in sp.get("cations", []):
        for s in sp.get("solvents", []):
            for n in (ns.get(c, []) if isinstance(ns, dict) else ns):
                cand = [[c] + [s] * int(n)]
                if scfg.get("with_anion", True):
                    cand += [[c] + [s] * (int(n) - 1) + [a] for a in sp.get("anions", [])]
                for nm in cand:
                    if sum(build.n_atoms(cat[x]) for x in nm) <= max_atoms:
                        names_all.append(nm)
    return names_all


def cluster_md(ctx, scfg, delta_q):
    """First-solvation-shell clusters: xTB MD at several temperatures, frames picked by greedy max-min RMSD."""
    names_all = md_compositions(ctx.cfg, ctx.cat, scfg)
    kmax = int(scfg.get("max_clusters", 20))
    if len(names_all) > kmax:
        names_all = [names_all[i] for i in ctx.rng.permutation(len(names_all))[:kmax]]
    T_K = [int(T) for T in scfg.get("T_K", [300, 450, 600])]
    time_ps = float(scfg.get("time_ps", 10))
    dump_fs = float(scfg.get("dump_fs", 50))
    equil_ps = float(scfg.get("equil_ps", 2))
    n_frames = int(scfg.get("n_frames", 10))
    n = 0
    for names in names_all:
        getattr(ctx, "progress", lambda *a: None)(build.composition(names))
        q = _formal(ctx, names) + int(delta_q)
        if not _q_ok(ctx, q, names, delta_q):
            continue
        cl = build.assemble(names, ctx.cat, ctx.workdir, q,
                            np.random.default_rng(int(ctx.cfg.get("seed", 0)) + zlib.crc32(f"md{names}{q}".encode())))
        members = [idx for _, idx in cl["members"]]
        if build.fragmented(cl["positions"], members):
            warnings.warn(f"cluster_md {names} q={q} fell apart in the xtb opt; skipped")
            continue
        symbols, pos = cl["symbols"], np.asarray(cl["positions"], float)
        pos = pos - pos.mean(axis=0)  # the xtb wall sphere is centred at the origin
        R = float(np.linalg.norm(pos, axis=1).max()) + 2.0
        skip = int(np.ceil(equil_ps / (dump_fs / 1000.0)))
        for T in T_K:
            r = xtb.run_xtb(symbols, pos, q, _uhf(symbols, q), ctx.workdir, timeout=_timeout(ctx),
                            md=dict(T_K=T, time_ps=time_ps, dump_fs=dump_fs, wall_R_A=R))
            if not (r.get("ok") and r.get("frames")):
                err = (str(r.get("err") or "no frames").strip().splitlines() or ["?"])[-1]
                warnings.warn(f"cluster_md {list(names)} q={q} T={T}: xtb md failed ({err}); skipped")
                continue
            raw = [np.asarray(f, float) for f in r["frames"]]
            idx = [k for k in range(skip, len(raw))
                   if not build.fragmented(raw[k], members)
                   and not build.fragmented(raw[k], _components(symbols, raw[k]))
                   and not _too_close(symbols, raw[k], 0.6)]
            cand = [raw[k] for k in idx]
            picked = _maxmin(cand, n_frames)
            for j, p in enumerate(picked):
                if n >= _limit(scfg):
                    return
                n += 1
                rmsd = None if j == 0 else min(_rmsd(cand[p], cand[h]) for h in picked[:j])
                yield _rec(ctx, symbols, cand[p], q, "cluster_md", names, "sp", _bs(ctx, names, delta_q),
                           {"T_K": float(T), "frame": int(idx[p]), "rmsd_min_A": rmsd, "names": list(names)}
                           | _etemp(r))


def nanoreactor(ctx, scfg, delta_q):
    """xTB nanoreactor MD frames (default off)."""
    T_K = int(scfg.get("T_K", 1500))
    ncyc = int(scfg.get("ncyc", 2))
    every = int(scfg.get("every_n_frames", 50))
    n = 0
    for combo in _combos(ctx):
        getattr(ctx, "progress", lambda *a: None)(build.composition(list(combo)))
        q = _formal(ctx, combo) + int(delta_q)
        if not _q_ok(ctx, q, combo, delta_q):
            continue
        mols = [(list(s), p) for s, p in (species.geometry(ctx.cat[nm], ctx.workdir)[:2] for nm in combo)]
        label = "_".join(combo)
        rundir = ctx.workdir / "nanoreactor" / f"{label}_q{q}_T{T_K}_c{ncyc}_s{int(ctx.cfg.get('seed', 0))}"
        xtb.build(rundir, mols, label, q, _uhf([s for s, _ in mols for s in s], q), T_K,
                  int(ctx.cfg.get("seed", 0)), ncyc)
        traj = xtb.run(rundir)
        frames = ase_read(str(traj), ":")
        frames = frames if isinstance(frames, list) else [frames]
        for k in range(0, len(frames), every):
            if n >= _limit(scfg):
                return
            a = frames[k]
            yield _rec(ctx, a.get_chemical_symbols(), a.get_positions(), q, "nanoreactor", combo, "sp", False,
                       {"combo": list(combo), "frame": k})
            n += 1


SAMPLERS = {"monomer": monomer, "atoms": atoms, "pair_scan": pair_scan, "cluster": cluster, "cluster_md": cluster_md, "bond_scan1": bond_scan1, "bond_scan2": bond_scan2,
            "reaction": reaction, "reaction2": reaction2, "hop": hop, "nanoreactor": nanoreactor}
