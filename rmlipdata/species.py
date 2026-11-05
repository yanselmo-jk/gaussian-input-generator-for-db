"""Species catalog (anions, solvents, Li+/Na+) and their embedded + xTB-optimised geometries.

`embed`, `_forcefield` copied from 7_apply_distortion_to_molecule/src/moldistort/geometry.py.
"""
import csv
import hashlib
import os
import re
import warnings
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDetermineBonds

from . import xtb
from .g16 import BASIS, n_electrons

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
    group: str = ""
    molfile: str | None = None


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
    start = mol.GetConformer().GetPositions()
    _forcefield(mol).Minimize(maxIts=2000)
    if _min_heavy(mol) < 1.0:  # UFF types hypervalent P/S as sp3 and collapses e.g. PF5 (F-F 0.87 A): keep ETKDG
        mol.GetConformer().SetPositions(start)
    return mol


def _min_heavy(mol):
    """Shortest distance between two non-H atoms (shortest real bond, C#O, is 1.13 A)."""
    x = mol.GetConformer().GetPositions()[[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]]
    d = np.linalg.norm(x[:, None] - x[None], axis=-1)
    return d[~np.eye(len(x), dtype=bool)].min(initial=np.inf)


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


def _checked(rec, name, smiles, charge, kind, group=""):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"{name}: cannot parse SMILES {smiles!r}")
    got = Chem.GetFormalCharge(mol)
    if got != charge:
        raise ValueError(f"{name}: SMILES formal charge {got} != tabulated {charge}")
    return Species(name, smiles, charge, kind, group)


def _builtin():
    out = {}
    for r in _tsv(DATA / "anions.tsv"):
        if r["note"].strip().upper().startswith("EXCLUDED"):
            continue
        name = r["abbrev"].strip()
        out[name] = _checked(r, name, r["smiles"].strip(), int(r["charge"]), "anion", "anion")
    for r in _tsv(DATA / "solvents.tsv"):
        name = _abbrev(r["name"].strip())
        out[name] = _checked(r, name, r["smiles"].strip(), 0, "solvent", r["family"].strip())
    for r in _tsv(DATA / "extra.tsv"):
        name = r["name"].strip()
        if name in out:
            raise ValueError(f"extra.tsv: {name} already in the catalog")
        out[name] = _checked(r, name, r["smiles"].strip(), int(r["charge"]), r["kind"].strip(),
                             r["group"].strip())
    for name, smiles in CATIONS.items():
        out[name] = _checked(None, name, smiles, 1, "cation", "metal_cation")
    return out


def home() -> Path:
    """$RMLIPDATA_HOME, else ${XDG_DATA_HOME:-~/.local/share}/rmlipdata."""
    if os.environ.get("RMLIPDATA_HOME"):
        return Path(os.environ["RMLIPDATA_HOME"])
    return Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share") / "rmlipdata"


def _user_rows():
    p = home() / "user_species.tsv"
    return p, (_tsv(p) if p.exists() else [])


@lru_cache(maxsize=1)
def catalog() -> dict:
    """name -> Species; anions keyed by abbrev, solvents by the parenthesised abbreviation."""
    out = _builtin()
    _, rows = _user_rows()
    for r in rows:
        name = r["name"].strip()
        if name in out:
            raise ValueError(f"user species {name!r} collides with a built-in species")
        out[name] = Species(name, r["smiles"].strip(), int(r["charge"]), r["kind"].strip(), "user",
                            str(home() / "species" / f"{name}.mol"))
    return out


NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_+-]*$")


def _mol_from_file(path):
    """Read a .mol/.sdf/.mol2/.xyz/.pdb file into an RDKit mol with bonds where possible."""
    suffix = path.suffix.lower()
    if suffix in (".mol", ".sdf"):
        mol = Chem.MolFromMolFile(str(path), removeHs=False)
        if mol is None:
            raise ValueError(f"cannot read {path}")
        return mol, Chem.GetFormalCharge(mol)
    if suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(path), removeHs=False)
        if mol is None:
            raise ValueError(f"cannot read {path}")
        return mol, Chem.GetFormalCharge(mol)
    if suffix not in (".xyz", ".pdb"):
        raise ValueError(f"unsupported file type {path.suffix!r}: use .mol, .mol2, .xyz or .pdb")
    if suffix == ".pdb":
        pdb = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False)
        if pdb is None or pdb.GetNumConformers() == 0:
            raise ValueError(f"cannot read atoms and coordinates from {path}")
        conf = pdb.GetConformer()
        symbols = [a.GetSymbol() for a in pdb.GetAtoms()]
        xyz = [(s, conf.GetAtomPosition(i)) for i, s in enumerate(symbols)]
    else:
        symbols, positions = _read_xyz(path)
        xyz = list(zip(symbols, positions))
    block = f"{len(xyz)}\n{path.name}\n" + "".join(
        f"{s} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n" for s, p in xyz)
    mol = Chem.MolFromXYZBlock(block)
    if mol is None:
        raise ValueError(f"cannot read atoms and coordinates from {path}")
    return mol, None


def _save_user_row(name, smiles, charge, kind, source):
    p, rows = _user_rows()
    p.parent.mkdir(parents=True, exist_ok=True)
    rows = [r for r in rows if r["name"].strip() != name] + [
        dict(name=name, smiles=smiles, charge=str(charge), kind=kind, source=str(source),
             added=datetime.now().isoformat(timespec="seconds"))]
    tmp = p.with_name(f"{p.name}.{os.getpid()}.tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "smiles", "charge", "kind", "source", "added"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, p)  # atomic: a crash never truncates user_species.tsv


def add_user_species(path, name, charge=None, kind=None, force=False) -> Species:
    """Validate a mol/mol2/xyz/pdb file, store it as $RMLIPDATA_HOME/species/<name>.mol, refresh catalog."""
    path = Path(path)
    if not NAME_RE.match(name):
        raise ValueError(f"invalid species name {name!r}: use letters, digits, _ + -")
    if name in _builtin():
        raise ValueError(f"{name!r} is a built-in species name")
    _, rows = _user_rows()
    if any(r["name"].strip() == name for r in rows) and not force:
        raise ValueError(f"user species {name!r} already exists; pass force=True to overwrite")
    if not path.is_file():
        raise ValueError(f"no such file: {path}")

    mol, file_charge = _mol_from_file(path)
    if file_charge is None:  # xyz/pdb: perceive bonds and formal charges
        charge = 0 if charge is None else int(charge)
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=charge)
        except Exception as exc:
            raise ValueError(f"{name}: cannot perceive bonds from {path}; "
                             f"supply .mol or .mol2 with bonds ({exc})") from exc
    else:
        if charge is not None and int(charge) != file_charge:
            raise ValueError(f"{name}: charge {int(charge)} != formal charge {file_charge} in {path}")
        charge = file_charge

    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    missing = sorted(set(symbols) - set(BASIS))
    if missing:
        raise ValueError(f"{name}: no def2-TZVPPD basis for {missing}")
    if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
        raise ValueError(f"{name}: {path} has no 3D coordinates")
    if sum(a.GetNumImplicitHs() for a in mol.GetAtoms()):
        warnings.warn(f"{name}: implicit hydrogens in {path}; adding them with coordinates")
        mol = Chem.AddHs(mol, addCoords=True)
    if len(Chem.GetMolFrags(mol)) > 1:
        warnings.warn(f"{name}: {path} contains more than one fragment")
    if kind is None:
        kind = "anion" if charge < 0 else "solvent"
    if kind not in ("anion", "solvent"):
        raise ValueError(f"{name}: kind must be 'anion' or 'solvent', got {kind!r}")

    dest = home() / "species" / f"{name}.mol"
    dest.parent.mkdir(parents=True, exist_ok=True)
    Chem.MolToMolFile(mol, str(dest))
    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
    _save_user_row(name, smiles, charge, kind, path.absolute())
    catalog.cache_clear()
    return Species(name, smiles, charge, kind, "user", str(dest))


def _read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    rows = [l.split() for l in lines[2:2 + n]]
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])


def geometry(sp: Species, workdir, charge=None):
    """(symbols, positions, mol): RDKit embed, then xTB opt. Cached at workdir/species/<name>_q<charge>.xyz.

    User species start from their stored .mol coordinates, and their cache name carries a hash of the
    molblock so re-adding the same name never reuses a stale geometry. mol keeps H in the same atom
    order as symbols (topology only). Single atoms sit at the origin.
    """
    charge = sp.charge if charge is None else int(charge)
    if sp.molfile:
        mol = Chem.MolFromMolFile(sp.molfile, removeHs=False)
        if mol is None:
            raise ValueError(f"{sp.name}: cannot read {sp.molfile}")
        tag = "_" + hashlib.sha1(Chem.MolToMolBlock(mol).encode()).hexdigest()[:8]
    else:
        mol = embed(sp.smiles)
        tag = ""
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    if mol.GetNumAtoms() == 1:
        return symbols, np.zeros((1, 3)), mol
    cache = Path(workdir) / "species" / f"{sp.name}_q{charge}{tag}.xyz"
    if cache.exists():
        symbols, positions = _read_xyz(cache)
        return symbols, positions, mol
    positions = mol.GetConformer().GetPositions()
    r = xtb.run_xtb(symbols, positions, charge, n_electrons(symbols, charge) % 2, workdir, mode="opt")
    if r["ok"] and "positions" in r:
        positions = np.asarray(r["positions"], float)
        cache.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache.with_name(f"{cache.name}.{os.getpid()}.tmp")
        tmp.write_text(f"{len(symbols)}\n{sp.name} q={charge}\n"
                       + "".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n"
                                 for s, (x, y, z) in zip(symbols, positions)))
        os.replace(tmp, cache)  # atomic: parallel generate stages may build the same species
        symbols, positions = _read_xyz(cache)  # same rounded floats as later cache hits, so reruns reproduce ids
    else:
        warnings.warn(f"{sp.name}: xtb opt failed, keeping force-field geometry ({r.get('err', '')[:120]})")
    return symbols, positions, mol
