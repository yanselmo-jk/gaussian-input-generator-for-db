"""Cluster enumeration, composition strings and donor-directed assembly of Li+/Na+ electrolyte clusters.

`_rodrigues` copied from 9_build_automatic_reaction_search_agent/rxnagent/core/solvation.py.
"""
import warnings
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement

import numpy as np
from ase import Atoms
from rdkit import Chem

from . import species, xtb

DONORS = {"O", "N", "F"}
D_M = {"Li": 2.0, "Na": 2.35}
_TETRA = np.array([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], float)
TETRA = _TETRA / np.linalg.norm(_TETRA, axis=1)[:, None]


# from 9_build_automatic_reaction_search_agent/rxnagent/core/solvation.py
def _rodrigues(a, b):
    """Rotation matrix taking vector a onto vector b."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    v, c = np.cross(a, b), float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-9:
        if c > 0:
            return np.eye(3)
        axis = np.cross(a, [1.0, 0.0, 0.0])
        if np.linalg.norm(axis) < 1e-9:
            axis = np.cross(a, [0.0, 1.0, 0.0])
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.eye(3) + vx + vx @ vx / (1.0 + c)


@lru_cache(maxsize=256)
def _n_atoms_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.AddHs(mol).GetNumAtoms() if mol is not None else 0


def n_atoms(sp):
    """Atom count including H, from the SMILES (no embedding or xtb)."""
    return _n_atoms_smiles(sp.smiles)


def composition(names):
    """Counts per name, cations then anions then solvents, alphabetical: ["Li","EC","EC","PF6"] -> 'Li1 PF6_1 EC2'."""
    cat = species.catalog()
    order = {"cation": 0, "anion": 1, "solvent": 2}
    counts = Counter(names)
    keys = sorted(counts, key=lambda n: (order.get(cat[n].kind, 3), n))
    return " ".join(f"{n}{counts[n]}" if not n[-1].isdigit() else f"{n}_{counts[n]}" for n in keys)


def combos(cat, cfg_species, cfg_cluster, max_atoms, rng):
    """All cation/component multisets within max_atoms, deterministically ordered then subsampled by rng."""
    cations = list(cfg_species.get("cations", []))
    comps = list(cfg_species.get("anions", [])) + list(cfg_species.get("solvents", []))
    min_c = int(cfg_cluster.get("min_components", 1))
    max_c = int(cfg_cluster.get("max_components", 3))
    max_cat = int(cfg_cluster.get("max_cations", 2))
    max_combos = int(cfg_cluster.get("max_combos", 10 ** 9))
    out = []
    for nc in range(max_cat + 1):
        for cats in combinations_with_replacement(cations, nc):
            for ncomp in range(min_c, max_c + 1):
                for cs in combinations_with_replacement(comps, ncomp):
                    names = tuple(cats) + tuple(cs)
                    if sum(n_atoms(cat[n]) for n in names) <= max_atoms:
                        out.append(names)
    order = rng.permutation(len(out))
    return [out[i] for i in order[:max_combos]]


def _place(pos, center, used, placed, distance, donors, rng):
    """Orient a block so COM->donor points at the cation; donor at `distance` along a free tetrahedral direction."""
    com = pos.mean(axis=0)
    for _ in range(200):
        donor = int(donors[rng.integers(len(donors))])
        free = [i for i in range(len(TETRA)) if i not in used] or list(range(len(TETRA)))
        di = free[int(rng.integers(len(free)))]
        u = TETRA[di]
        a = com - pos[donor]
        if np.linalg.norm(a) < 1e-8:
            a = np.array([0.0, 0.0, 1.0])
        new = (pos - pos[donor]) @ _rodrigues(a, u).T + center + distance * u
        if all(np.min(np.linalg.norm(new[:, None] - q[None], axis=2)) >= 1.6 for q in placed.values()):
            used.append(di)
            return new
    return None


def _blocks(names, workdir):
    blocks = []
    for n in names:
        symbols, positions, _ = species.geometry(species.catalog()[n], workdir)
        blocks.append((n, list(symbols), np.asarray(positions, float)))
    return blocks


def _fallback_pack(blocks, rng):
    mols = [Atoms(symbols=s, positions=p) for _, s, p in blocks]
    whole = Atoms(sum((s for _, s, _ in blocks), []))
    at = xtb.pack(mols, rng, xtb.r_liq(whole) + 1.0)
    sizes = [len(s) for _, s, _ in blocks]
    out, k = [], 0
    for n in sizes:
        out.append(np.asarray(at.positions[k:k + n], float))
        k += n
    return _contact(out)


def _contact(blocks, target=2.5):
    """Move each block (rigidly, toward the blocks before it) until its shortest distance to them is target A.
    Gas-phase xTB opt does not pull neutral molecules together from a few A apart."""
    out = [blocks[0]]
    for b in blocks[1:]:
        ref = np.vstack(out)
        v = ref.mean(axis=0) - b.mean(axis=0)
        span = float(np.linalg.norm(v))
        if span < 1e-8:
            out.append(b)
            continue
        u = v / span
        dmin = lambda t: np.linalg.norm((b + t * u)[:, None] - ref[None], axis=-1).min()
        lo, hi = -span, span  # largest shift toward ref that keeps the contact >= target (bisection)
        if dmin(0.0) <= target:
            lo, hi = -span, 0.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if dmin(mid) > target else (lo, mid)
        out.append(b + lo * u)
    return out


def assemble(names, cat, workdir, charge, rng):
    """Place members around cations (or pack), then xTB-optimise. Returns symbols/positions/members/composition."""
    names = list(names)
    blocks = _blocks(names, workdir)
    is_cation = [cat[n].kind == "cation" for n in names]
    cation_idx = [k for k, f in enumerate(is_cation) if f]
    positions = [None] * len(names)
    if cation_idx:
        used = {}
        for k, i in enumerate(cation_idx):
            positions[i] = np.zeros(3) if k == 0 else np.array([4.5 * k, 0.0, 0.0])
            used[i] = []
        for j, k in enumerate([k for k, f in enumerate(is_cation) if not f]):
            target = cation_idx[j % len(cation_idx)]
            donors = [a for a, s in enumerate(blocks[k][1]) if s in DONORS]
            placed = {i: p for i, p in enumerate(positions) if p is not None}
            got = _place(blocks[k][2], positions[target], used[target], placed,
                         D_M.get(names[target], 2.0), donors, rng) if donors else None
            if got is None:
                positions = _fallback_pack(blocks, rng)
                break
            positions[k] = got
    else:
        positions = _fallback_pack(blocks, rng)

    symbols, members, k = [], [], 0
    for (n, sym, _), p in zip(blocks, positions):
        symbols += sym
        members.append((n, list(range(k, k + len(sym)))))
        k += len(sym)
    coords = np.vstack(positions)

    r = xtb.run_xtb(symbols, coords, int(charge), _uhf(symbols, charge), workdir, mode="opt")
    out = {"symbols": symbols, "positions": coords, "members": members, "composition": composition(names)}
    if not (r["ok"] and "positions" in r):
        warnings.warn(f"assemble {list(names)}: xtb opt failed, keeping pre-opt geometry ({r.get('err', '')[:120]})")
    elif r.get("etemp_K") and fragmented(r["positions"], [idx for _, idx in members]):
        # smeared opt converged but the members flew apart (no bound gas-phase minimum, e.g. two anions):
        # the contact geometry is the useful PES point
        warnings.warn(f"assemble {list(names)} q={charge}: --etemp opt separated the members; keeping pre-opt geometry")
    else:
        out["positions"] = np.asarray(r["positions"], float)
        if r.get("etemp_K"):
            out["etemp_K"] = r["etemp_K"]
    return out


def fragmented(positions, members, cutoff=4.0):
    """True if some member (list of atom-index lists) has no atom within cutoff A of any other member."""
    if len(members) < 2:
        return False
    pos = np.asarray(positions, float)
    return any(min(float(np.linalg.norm(pos[a][:, None] - pos[b][None], axis=-1).min())
                   for b in members if b is not a) > cutoff for a in members)


def _uhf(symbols, charge):
    from .g16 import n_electrons
    return n_electrons(symbols, int(charge)) % 2
