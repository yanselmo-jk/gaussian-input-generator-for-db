import warnings

import numpy as np
import pytest
import yaml
from ase.io import read
from rdkit import Chem

from rmlipdata import build, cli, species

METHANOL_MOL2 = """@<TRIPOS>MOLECULE
methanol
 6 5 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 MET       0.0000
      2 O1          1.4300    0.0000    0.0000 O.3       1 MET       0.0000
      3 H1         -0.3600    1.0300    0.0000 H         1 MET       0.0000
      4 H2         -0.3600   -0.5100    0.8900 H         1 MET       0.0000
      5 H3         -0.3600   -0.5100   -0.8900 H         1 MET       0.0000
      6 H4          1.8700    0.8900    0.0000 H         1 MET       0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     1     3    1
     3     1     4    1
     4     1     5    1
     5     2     6    1
"""


@pytest.fixture(autouse=True)
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("RMLIPDATA_HOME", str(tmp_path / "home"))
    species.catalog.cache_clear()
    yield
    species.catalog.cache_clear()


def _ec():
    return species.embed("O=C1OCCO1")


def _write_xyz(path, mol):
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), "test"] + [
        f"{a.GetSymbol()} {conf.GetAtomPosition(a.GetIdx()).x} "
        f"{conf.GetAtomPosition(a.GetIdx()).y} {conf.GetAtomPosition(a.GetIdx()).z}"
        for a in mol.GetAtoms()]
    path.write_text("\n".join(lines) + "\n")


def _ec_mol(tmp_path):
    p = tmp_path / "ec.mol"
    Chem.MolToMolFile(_ec(), str(p))
    return p


def _check(name, charge, kind, n_atoms):
    sp = species.catalog()[name]
    assert (sp.charge, sp.kind, build.n_atoms(sp)) == (charge, kind, n_atoms)
    assert sp.group == "user" and sp.molfile and sp.molfile.endswith(f"{name}.mol")


def test_add_ec_mol(tmp_path):
    species.add_user_species(_ec_mol(tmp_path), "ECmol")
    _check("ECmol", 0, "solvent", 10)


def test_add_ec_xyz(tmp_path):
    p = tmp_path / "ec.xyz"
    _write_xyz(p, _ec())
    species.add_user_species(p, "ECxyz")
    _check("ECxyz", 0, "solvent", 10)


def test_add_ec_pdb(tmp_path):
    p = tmp_path / "ec.pdb"
    Chem.MolToPDBFile(_ec(), str(p))
    species.add_user_species(p, "ECpdb")
    _check("ECpdb", 0, "solvent", 10)


def test_add_methanol_mol2(tmp_path):
    p = tmp_path / "methanol.mol2"
    p.write_text(METHANOL_MOL2)
    species.add_user_species(p, "MeOH")
    _check("MeOH", 0, "solvent", 6)


def test_add_formate_xyz_charge(tmp_path):
    p = tmp_path / "formate.xyz"
    p.write_text("4\nformate\nC 0 0 0\nO 0 0 1.2\nO 0 0 -1.2\nH 1.0 0 0\n")
    species.add_user_species(p, "HCOO", charge=-1)
    _check("HCOO", -1, "anion", 4)


def test_collision_and_force(tmp_path):
    mol = _ec_mol(tmp_path)
    with pytest.raises(ValueError):
        species.add_user_species(mol, "EC")
    species.add_user_species(mol, "Mine")
    with pytest.raises(ValueError):
        species.add_user_species(mol, "Mine")
    species.add_user_species(mol, "Mine", force=True)
    assert sum(1 for sp in species.catalog().values() if sp.name == "Mine") == 1


def test_element_without_basis(tmp_path):
    p = tmp_path / "hbr.xyz"
    p.write_text("2\nhbr\nBr 0 0 0\nH 0 0 1.5\n")
    with pytest.raises(ValueError, match="Br"):
        species.add_user_species(p, "HBr")


def test_geometry_keeps_file_atom_order(tmp_path, monkeypatch):
    p = _ec_mol(tmp_path)
    species.add_user_species(p, "Order")
    monkeypatch.setattr(species.xtb, "run_xtb", lambda *a, **k: {"ok": False, "err": "skipped"})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        symbols, positions, mol = species.geometry(species.catalog()["Order"], tmp_path)
    assert symbols == [a.GetSymbol() for a in Chem.MolFromMolFile(str(p), removeHs=False).GetAtoms()]
    assert positions.shape == (10, 3) and len(symbols) == 10


def test_cli_species_lists_user(tmp_path, capsys):
    p = _ec_mol(tmp_path)
    species.add_user_species(p, "Mine")
    assert cli.main(["species"]) == 0
    out = capsys.readouterr().out
    assert "PF6" in out and "EC" in out and "Mine" in out


def test_cli_view_save(tmp_path):
    species.add_user_species(_ec_mol(tmp_path), "Mine")
    out = tmp_path / "view.xyz"
    assert cli.main(["view", "Mine", "--save", str(out)]) == 0
    atoms = read(out)
    assert len(atoms) == 10 and atoms.get_chemical_formula() == "C3H4O3"


def test_generate_monomer_with_user_species(tmp_path):
    species.add_user_species(_ec_mol(tmp_path), "Mine")
    cfg = {
        "out": str(tmp_path / "run"),
        "seed": 0,
        "species": {"cations": ["Li"], "anions": [], "solvents": ["Mine"]},
        "charge": {"delta_q": [0], "mult_offsets": [0, 2]},
        "g16": {"nproc": 1, "mem": "1GB"},
        "limits": {"max_atoms": 60},
        "xtb": {"threads": 2, "timeout_s": 600},
        "cluster": {"min_components": 1, "max_components": 1, "max_cations": 1, "max_combos": 5},
        "samplers": {"monomer": {"on": True, "n_normal_mode": 1, "n_distort": 1,
                                 "T_K": 500, "max_structures": 4}},
    }
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0
    assert list((tmp_path / "run" / "gjf" / "monomer").glob("*.gjf"))
