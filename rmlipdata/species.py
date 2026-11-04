"""Species catalog (anions, solvents, Li+/Na+) and their embedded + xTB-optimised geometries.

`embed`, `_forcefield` copied from 7_apply_distortion_to_molecule/src/moldistort/geometry.py.
"""
import csv
import re
import warnings
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from . import xtb
from .g16 import n_electrons

RDLogger.DisableLog("rdApp.*")
DATA = Path(__file__).parent / "data"
CATIONS = {"Li": "[Li+]", "Na": "[Na+]"}
SUFFIX = re.compile(r"\(([^()]*)\)\s*$")


@dataclass
class Species:
    name: str
    smiles: str
    charge: int
    kind: str  # cation | anion | solvent


# from 7_apply_distortion_to_molecule/src/moldistort/geometry.py
def embed(smiles, seed=0):
    """SMILES -> 3D molecule with H (ETKDGv3 + MMFF, UFF if MMFF lacks parameters)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"cannot parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        # distance-geometry start fails on some cages/macrocycles; random start handles them
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError(f"3D embedding failed: {smiles}")
    _forcefield(mol).Minimize(maxIts=2000)
    return mol


def _forcefield(mol):
    if AllChem.MMFFHasAllMoleculeParams(mol):
        return AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    return AllChem.UFFGetMoleculeForceField(mol)  # Li+ complexes land here


def _tsv(path):
    lines = [l for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def _abbrev(name):
    m = SUFFIX.search(name)
    return m[1] if m else re.sub(r"[^A-Za-z0-9]", "_", name)


def _checked(rec, name, smiles, charge, kind):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"{name}: cannot parse SMILES {smiles!r}")
    got = Chem.GetFormalCharge(mol)
    if got != charge:
        raise ValueError(f"{name}: SMILES formal charge {got} != tabulated {charge}")
    return Species(name, smiles, charge, kind)


@lru_cache(maxsize=1)
def catalog() -> dict:
    """name -> Species; anions keyed by abbrev, solvents by the parenthesised abbreviation."""
    out = {}
    for r in _tsv(DATA / "anions.tsv"):
        if r["note"].strip().upper().startswith("EXCLUDED"):
            continue
        name = r["abbrev"].strip()
        out[name] = _checked(r, name, r["smiles"].strip(), int(r["charge"]), "anion")
    for r in _tsv(DATA / "solvents.tsv"):
        name = _abbrev(r["name"].strip())
        out[name] = _checked(r, name, r["smiles"].strip(), 0, "solvent")
    for name, smiles in CATIONS.items():
        out[name] = _checked(None, name, smiles, 1, "cation")
    return out


def _read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    rows = [l.split() for l in lines[2:2 + n]]
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])


def geometry(sp: Species, workdir, charge=None):
    """(symbols, positions, mol): RDKit embed, then xTB opt. Cached at workdir/species/<name>_q<charge>.xyz.

    mol keeps H in the same atom order as symbols (topology only). Single atoms sit at the origin.
    """
    charge = sp.charge if charge is None else int(charge)
    mol = embed(sp.smiles)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    if mol.GetNumAtoms() == 1:
        return symbols, np.zeros((1, 3)), mol
    cache = Path(workdir) / "species" / f"{sp.name}_q{charge}.xyz"
    if cache.exists():
        symbols, positions = _read_xyz(cache)
        return symbols, positions, mol
    positions = mol.GetConformer().GetPositions()
    r = xtb.run_xtb(symbols, positions, charge, n_electrons(symbols, charge) % 2, workdir, mode="opt")
    if r["ok"] and "positions" in r:
        positions = np.asarray(r["positions"], float)
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(f"{len(symbols)}\n{sp.name} q={charge}\n"
                         + "".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n"
                                   for s, (x, y, z) in zip(symbols, positions)))
        symbols, positions = _read_xyz(cache)  # same rounded floats as later cache hits, so reruns reproduce ids
    else:
        warnings.warn(f"{sp.name}: xtb opt failed, keeping force-field geometry ({r.get('err', '')[:120]})")
    return symbols, positions, mol
