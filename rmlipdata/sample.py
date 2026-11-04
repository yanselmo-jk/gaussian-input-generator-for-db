"""Structure samplers: monomers, pairs, clusters, bond scans, pseudo-TS reactions, cation hops, nanoreactor.

`rotors`, `bends`, `min_far_distance`, `sample`, `_forcefield` copied from
7_apply_distortion_to_molecule/src/moldistort/geometry.py.
"""
import warnings
import zlib
from itertools import combinations

import numpy as np
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


def _rec(ctx, symbols, positions, charge, sampler, names, job, bs, params):
    return dict(symbols=[str(s) for s in symbols], positions=np.asarray(positions, float), charge=int(charge),
                sampler=sampler, composition=build.composition(list(names)), job=job, bs=bool(bs),
                params=params)


def _with_positions(mol, positions):
    """Copy of the topology mol carrying the given conformer."""
    m = Chem.Mol(mol)
    m.GetConformer().SetPositions(np.asarray(positions, float))
    return m


def _q_ok(ctx, q):
    """Total charge inside cfg charge.range (default [-2, 2]); beyond it gas-phase clusters lose the extra
    electron (anions) or Coulomb-explode (cations)."""
    lo, hi = ctx.cfg.get("charge", {}).get("range", (-2, 2))
    return lo <= q <= hi


def _cluster(ctx, combo, delta_q):
    q = _formal(ctx, combo) + int(delta_q)
    key = (tuple(combo), q)
    if not _q_ok(ctx, q):
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
    mem = [idx for _, idx in cl["members"]]
    if len(mem) < 2:
        return False
    return any(min(_min_between(cl["positions"], a, b) for b in mem if b is not a) > cutoff for a in mem)


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
            if n >= _limit(scfg):
                return
            sp = ctx.cat[name]
            q = sp.charge + int(delta_q)
            if not _q_ok(ctx, q):
                continue
            if sp.kind == "cation":
                if q > 1:  # bare Li2+/Na2+ is not an electrolyte state
                    continue
                yield _rec(ctx, [name], np.zeros((1, 3)), q, "monomer", [name], "sp", False,
                           {"kind": "cation", "index": 0})
                n += 1
                continue
            symbols, positions, mol = species.geometry(sp, ctx.workdir, charge=q)
            yield _rec(ctx, symbols, positions, q, "monomer", [name], "opt", False, {"kind": "opt", "index": 0})
            n += 1
            if n >= _limit(scfg):
                return

            r = xtb.run_xtb(symbols, positions, q, _uhf(symbols, q), ctx.workdir, mode="hess",
                            timeout=_timeout(ctx))
            for k in range(int(scfg.get("n_normal_mode", 0))):
                if n >= _limit(scfg):
                    return
                pos = _normal_mode(symbols, positions, r.get("hessian"), float(scfg.get("T_K", 500)), ctx.rng) \
                    if r.get("ok") else None
                if pos is None:
                    break
                yield _rec(ctx, symbols, pos, q, "monomer", [name], "sp", False,
                           {"kind": "normal_mode", "T_K": float(scfg.get("T_K", 500)), "index": k})
                n += 1

            frames = sample(_with_positions(mol, positions), int(scfg.get("n_distort", 0)), rng=ctx.rng)
            for k, pos in enumerate(frames[1:int(scfg.get("n_distort", 0)) + 1]):
                if n >= _limit(scfg):
                    return
                yield _rec(ctx, symbols, pos, q, "monomer", [name], "sp", False,
                           {"kind": "distort", "sample": k})
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
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i:]]
    pairs = [p for p in pairs if not (ctx.cat[p[0]].kind == "cation" and ctx.cat[p[1]].kind == "cation")]
    max_pairs = int(scfg.get("max_pairs", 10 ** 9))
    if len(pairs) > max_pairs:
        pairs = [pairs[i] for i in ctx.rng.permutation(len(pairs))[:max_pairs]]
    n = 0
    for pair in pairs:
        for orient in range(int(scfg.get("n_orient", 1))):
            if n >= _limit(scfg):
                return
            q = _formal(ctx, pair) + int(delta_q)
            if not _q_ok(ctx, q):
                break
            cl = build.assemble(list(pair), ctx.cat, ctx.workdir, q, ctx.rng)
            (_, ia), (_, ib) = cl["members"]
            pos = cl["positions"]
            d0 = float(np.linalg.norm(_com(pos, ib) - _com(pos, ia)))
            if d0 < 1e-6:
                break
            u = (_com(pos, ib) - _com(pos, ia)) / d0
            if _min_between(pos, ia, ib) > 3.0:  # gas-phase xTB opt left a neutral pair apart: bring into contact
                pos = _shift_to(pos, ia, ib, u, 2.5, d0)
                r = xtb.run_xtb(cl["symbols"], pos, q, _uhf(cl["symbols"], q), ctx.workdir, mode="opt",
                                timeout=_timeout(ctx))
                if r.get("ok") and "positions" in r:
                    pos = np.asarray(r["positions"], float)
                d0 = float(np.linalg.norm(_com(pos, ib) - _com(pos, ia)))
                u = (_com(pos, ib) - _com(pos, ia)) / d0
            dc = _min_between(pos, ia, ib)  # contact = shortest inter-member distance
            for s in scfg.get("scales", [1.0]):
                if n >= _limit(scfg):
                    return
                moved = _shift_to(pos, ia, ib, u, s * dc, d0)
                yield _rec(ctx, cl["symbols"], moved, q, "pair_scan", pair, "sp", False,
                           {"pair": list(pair), "orient": orient, "scale": float(s), "contact_A": dc,
                            "d_min_A": _min_between(moved, ia, ib)})
                n += 1


def cluster(ctx, scfg, delta_q):
    """Assembled cluster opt + Gaussian-rattled sp frames."""
    n = 0
    for combo in _combos(ctx):
        if len(combo) < 2 or n >= _limit(scfg):
            continue
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
        yield _rec(ctx, cl["symbols"], cl["positions"], q, "cluster", combo, "opt", False,
                   {"combo": list(combo), "kind": "opt"})
        n += 1
        for k in range(int(scfg.get("n_distort", 0))):
            if n >= _limit(scfg):
                return
            pos = _rattle(cl["positions"], float(scfg.get("rattle_A", 0.1)), ctx.rng)
            if pos is None:
                continue
            yield _rec(ctx, cl["symbols"], pos, q, "cluster", combo, "sp", False,
                       {"combo": list(combo), "kind": "distort", "index": k})
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


def bond_scan(ctx, scfg, delta_q):
    """Relaxed xtb scans of one intramolecular bond per symmetry class per member species."""
    bs = bool(scfg.get("bs", False))
    npts = int(scfg.get("n_points", 12))
    max_bonds = int(scfg.get("max_bonds_per_cluster", 10 ** 9))
    n = 0
    for combo in _combos(ctx):
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
        symbols, pos = cl["symbols"], cl["positions"]
        seen, cand = set(), []
        for (name, idxs), (_, mol) in zip(cl["members"], _mols(ctx, combo)):
            if ctx.cat[name].kind == "cation":
                continue
            rank = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
            for b in mol.GetBonds():
                i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                key = (name, *sorted((rank[i], rank[j])))
                if key not in seen:
                    seen.add(key)
                    cand.append((idxs[i], idxs[j]))
        if len(cand) > max_bonds:
            cand = [cand[i] for i in ctx.rng.permutation(len(cand))[:max_bonds]]
        for i, j in cand:
            r0 = float(np.linalg.norm(pos[i] - pos[j]))
            r1 = 1.05 * (_cov(symbols, i) + _cov(symbols, j)) + 1.5
            r = xtb.run_xtb(symbols, pos, q, _uhf(symbols, q), ctx.workdir, mode="sp",
                            constraints=[(i, j, r0)], scan=dict(n=npts, ranges=[(0.85 * r0, r1)]),
                            timeout=_timeout(ctx))
            if not r.get("ok"):
                continue
            kmax = _until_jump(r["E_eV"])
            for k, (frame, energy) in enumerate(zip(r["frames"][:kmax], r["E_eV"][:kmax])):
                if n >= _limit(scfg):
                    return
                frame = np.asarray(frame, float)
                yield _rec(ctx, symbols, frame, q, "bond_scan", combo, "sp", bs,
                           {"bond": [i, j], "elements": [symbols[i], symbols[j]],
                            "r_A": float(np.linalg.norm(frame[i] - frame[j])), "index": k,
                            "E_xtb_eV": float(energy)})
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


def _reaction_candidates(ctx, cl, combo):
    """[(template, atoms, ranges)] over member pairs within 6 A."""
    infos = []
    for (name, idxs), (_, mol) in zip(cl["members"], _mols(ctx, combo)):
        if ctx.cat[name].kind == "cation":
            continue
        infos.append(dict(name=name, idxs=idxs, mol=mol))
    pos = cl["positions"]
    cand = []
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
                    tgt = 1.05 * (_cov(cl["symbols"], nu) + _cov(cl["symbols"], e))
                    cand.append(("add", [nu, e], [(d0, tgt)]))
                    carbonyl_o = _carbonyl_oxygens(b["mol"], e_local)
                    for x in b["mol"].GetAtomWithIdx(e_local).GetNeighbors():
                        if x.GetSymbol() not in ("O", "F") or x.GetIdx() in carbonyl_o:
                            continue
                        l = b["idxs"][x.GetIdx()]
                        r0 = float(np.linalg.norm(pos[e] - pos[l]))
                        tgt2 = 1.05 * (_cov(cl["symbols"], e) + _cov(cl["symbols"], l)) + 1.5
                        cand.append(("substitute", [nu, e, l], [(d0, tgt), (r0, tgt2)]))
    # transfer: X (F/H) bonded to Y in a non-cation member A -> cation, or an electrophile of another member
    cations = [k for (name, idxs) in cl["members"] if ctx.cat[name].kind == "cation" for k in idxs]
    for a in infos:
        for x_local in range(a["mol"].GetNumAtoms()):
            if a["mol"].GetAtomWithIdx(x_local).GetSymbol() not in ("F", "H"):
                continue
            for y_atom in a["mol"].GetAtomWithIdx(x_local).GetNeighbors():
                x, y = a["idxs"][x_local], a["idxs"][y_atom.GetIdx()]
                rxy = float(np.linalg.norm(pos[y] - pos[x]))
                tgt = 1.05 * (_cov(cl["symbols"], y) + _cov(cl["symbols"], x)) + 1.2
                for z in cations:
                    if np.linalg.norm(_com(pos, a["idxs"]) - pos[z]) > 6.0:
                        continue
                    dz = float(np.linalg.norm(pos[z] - pos[x]))
                    ztgt = 0.9 * (_cov(cl["symbols"], z) + _cov(cl["symbols"], x))
                    cand.append(("transfer", [y, x, z], [(rxy, tgt), (dz, ztgt)]))
                for b in infos:
                    if b is a or np.linalg.norm(_com(pos, a["idxs"]) - _com(pos, b["idxs"])) > 6.0:
                        continue
                    for e_local in _electrophiles(b["mol"]):
                        z = b["idxs"][e_local]
                        dz = float(np.linalg.norm(pos[z] - pos[x]))
                        ztgt = 1.0 * (_cov(cl["symbols"], z) + _cov(cl["symbols"], x))
                        cand.append(("transfer", [y, x, z], [(rxy, tgt), (dz, ztgt)]))
    return cand


def reaction(ctx, scfg, delta_q):
    """Pseudo-TS reaction paths (add / substitute / transfer): scan, keep the energy peak +/- k_peak."""
    bs = bool(scfg.get("bs", False))
    templates = set(scfg.get("templates", []))
    npts = int(scfg.get("n_points", 14))
    k_peak = int(scfg.get("k_peak", 2))
    max_per = int(scfg.get("max_per_cluster", 10 ** 9))
    n = 0
    for combo in _combos(ctx):
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
        cand = [c for c in _reaction_candidates(ctx, cl, combo) if c[0] in templates]
        if len(cand) > max_per:
            cand = [cand[i] for i in ctx.rng.permutation(len(cand))[:max_per]]
        for template, atoms, ranges in cand:
            if n >= _limit(scfg):
                return
            if template == "transfer":
                cons = [(atoms[0], atoms[1], ranges[0][0]), (atoms[2], atoms[1], ranges[1][0])]
            else:
                cons = [(atoms[k], atoms[k + 1], ranges[k][0]) for k in range(len(ranges))]
            r = xtb.run_xtb(cl["symbols"], cl["positions"], q, _uhf(cl["symbols"], q), ctx.workdir, mode="sp",
                            constraints=cons, scan=dict(n=npts, ranges=ranges), timeout=_timeout(ctx))
            if not r.get("ok"):
                continue
            E = r["E_eV"][:_until_jump(r["E_eV"])]
            if len(E) < 3:
                continue
            p = 1 + int(np.argmax(E[1:-1]))
            keep = sorted({0, len(E) - 1} | set(range(max(1, p - k_peak), min(p + k_peak, len(E) - 2) + 1)))
            for k in keep:
                if n >= _limit(scfg):
                    return
                yield _rec(ctx, cl["symbols"], r["frames"][k], q, "reaction", combo, "sp", bs,
                           {"template": template, "atoms": list(atoms), "ranges": [list(x) for x in ranges],
                            "index": int(k), "peak_index": int(p), "is_peak": bool(k == p),
                            "E_xtb_eV": float(E[k])})
                n += 1


def hop(ctx, scfg, delta_q):
    """Cation hopping path P1 -> P2 (in front of the target donor). At each point the cation's distances to the
    start and target donors (and their mutual distance) are constrained and the rest is xTB-relaxed. xtb $fix was
    not exact for these clusters (fixed atoms moved up to 2 A), so distance constraints are used. Frames where the
    cation is closer than 0.75 x covalent-radius sum to any atom are dropped."""
    npts = int(scfg.get("n_points", 10))
    max_per = int(scfg.get("max_per_cluster", 10 ** 9))
    n = 0
    for combo in _combos(ctx):
        q, cl = _cluster(ctx, combo, delta_q)
        if cl is None:
            continue
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
            near = min(donors, key=lambda z: np.linalg.norm(pos[m] - pos[z[3]]))
            targets = [z for z in donors if 2.6 < np.linalg.norm(pos[m] - pos[z[3]]) < 6.0]
            if len(targets) > max_per:
                targets = [targets[i] for i in ctx.rng.permutation(len(targets))[:max_per]]
            for mem, name, idxs, d2 in targets:
                mol = mols[mem][1]
                local = idxs.index(d2)
                nb = [idxs[a.GetIdx()] for a in mol.GetAtomWithIdx(local).GetNeighbors()]
                u = pos[d2] - (pos[nb].mean(axis=0) if nb else _com(pos, idxs))
                if np.linalg.norm(u) < 1e-8:
                    continue
                p2 = pos[d2] + build.D_M.get(cation, 2.0) * u / np.linalg.norm(u)
                d1 = int(near[3])
                p1 = pos[m].copy()
                rest = [i for i in range(len(symbols)) if i != m]
                rcov = np.array([_cov(symbols, i) for i in range(len(symbols))])
                cur = pos.copy()
                for t in np.linspace(0.0, 1.0, npts):
                    if n >= _limit(scfg):
                        return
                    pt = (1 - t) * p1 + t * p2
                    cur = _kabsch(cur, pos, rest)  # back into the frame where p1, p2 are defined
                    cur[m] = pt
                    cons = [(m, d1, float(np.linalg.norm(pt - pos[d1]))), (m, int(d2), float(np.linalg.norm(pt - pos[d2]))),
                            (d1, int(d2), float(np.linalg.norm(pos[d1] - pos[d2])))]
                    r = xtb.run_xtb(symbols, cur, q, _uhf(symbols, q), ctx.workdir, mode="opt",
                                    constraints=cons, timeout=_timeout(ctx))
                    if not (r.get("ok") and "positions" in r):
                        continue
                    cur = np.asarray(r["positions"], float)
                    dm = np.linalg.norm(cur[rest] - cur[m], axis=1)
                    if (dm < 0.75 * (rcov[m] + rcov[rest])).any():
                        continue  # cation runs into an atom: no useful PES point
                    yield _rec(ctx, symbols, cur, q, "hop", combo, "sp", False,
                               {"cation": cation, "cation_index": int(m), "from_donor": d1, "to_donor": int(d2), "t": float(t),
                                "d_from_A": float(np.linalg.norm(cur[m] - cur[d1])),
                                "d_to_A": float(np.linalg.norm(cur[m] - cur[d2]))})
                    n += 1


def nanoreactor(ctx, scfg, delta_q):
    """xTB nanoreactor MD frames (default off)."""
    T_K = int(scfg.get("T_K", 1500))
    ncyc = int(scfg.get("ncyc", 2))
    every = int(scfg.get("every_n_frames", 50))
    n = 0
    for combo in _combos(ctx):
        q = _formal(ctx, combo) + int(delta_q)
        if not _q_ok(ctx, q):
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


SAMPLERS = {"monomer": monomer, "pair_scan": pair_scan, "cluster": cluster, "bond_scan": bond_scan,
            "reaction": reaction, "hop": hop, "nanoreactor": nanoreactor}
