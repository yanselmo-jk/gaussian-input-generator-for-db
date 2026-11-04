import csv

import numpy as np
import yaml

from rmlipdata import build, cli, registry, species

TINY = {
    "out": "run",
    "seed": 0,
    "species": {"cations": ["Li"], "anions": [], "solvents": ["EC"]},
    "charge": {"delta_q": [0], "mult_offsets": [0, 2]},
    "g16": {"nproc": 1, "mem": "1GB"},
    "limits": {"max_atoms": 60},
    "xtb": {"threads": 2, "timeout_s": 600},
    "cluster": {"min_components": 1, "max_components": 1, "max_cations": 1, "max_combos": 5},
    "samplers": {
        "monomer": {"on": True, "n_normal_mode": 1, "n_distort": 1, "T_K": 500, "max_structures": 4},
        "pair_scan": {"on": True, "max_pairs": 1, "n_orient": 1, "scales": [0.9, 1.5], "max_structures": 2},
        "cluster": {"on": True, "n_distort": 1, "rattle_A": 0.1},
        "bond_scan": {"on": True, "n_points": 3, "max_bonds_per_cluster": 1, "bs": True, "max_structures": 3},
        "reaction": {"on": True, "templates": ["add", "substitute", "transfer"], "n_points": 3, "k_peak": 1,
                     "max_per_cluster": 1, "bs": True},
        "hop": {"on": True, "n_points": 3, "max_per_cluster": 1},
        "nanoreactor": {"on": False, "T_K": 1500, "ncyc": 2, "every_n_frames": 50},
    },
}


def test_catalog():
    cat = species.catalog()
    for name in ("EC", "PF6", "FSI", "Li", "Na"):
        assert name in cat
    assert "AsF6" not in cat
    assert cat["EC"].charge == 0
    assert cat["PF6"].charge == -1 and cat["FSI"].charge == -1
    assert cat["Li"].charge == 1 and cat["Na"].charge == 1


def test_combos():
    cat = species.catalog()
    sp = {"cations": ["Li"], "anions": ["PF6"], "solvents": ["EC"]}
    cluster = {"min_components": 1, "max_components": 2, "max_cations": 1, "max_combos": 1000}
    out = build.combos(cat, sp, cluster, 60, np.random.default_rng(0))
    assert len(out) == 10
    small = build.combos(cat, sp, cluster, 8, np.random.default_rng(0))
    assert 0 < len(small) < len(out)


def test_assemble_li_ec(tmp_path):
    cat = species.catalog()
    cl = build.assemble(["Li", "EC"], cat, tmp_path, 0, np.random.default_rng(0))
    assert cl["composition"] == "Li1 EC1"
    li = cl["symbols"].index("Li")
    d = [np.linalg.norm(cl["positions"][li] - cl["positions"][i])
         for i, s in enumerate(cl["symbols"]) if s == "O"]
    assert min(d) < 2.3


def _manifest_rows(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {(r["id"], r["gjf_path"], r["sampler"], str(r["n_atoms"]), str(r["charge"]), str(r["mult"]),
             str(r["bs"]), r["job"]) for r in rows}


def test_generate_end_to_end(tmp_path):
    cfg = dict(TINY, out=str(tmp_path / "run"))
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0

    db = tmp_path / "run" / "registry.sqlite"
    con = registry.connect(db)
    samplers = {r[0] for r in con.execute("SELECT DISTINCT sampler FROM structures")}
    assert {"monomer", "pair_scan", "cluster", "bond_scan", "reaction", "hop"} <= samplers
    rows = con.execute("SELECT id, gjf_path, sampler, n_atoms, charge, mult, bs, job FROM structures "
                       "WHERE status='generated'").fetchall()
    assert rows and all(r["gjf_path"] and __import__("pathlib").Path(r["gjf_path"]).exists() for r in rows)
    manifest = _manifest_rows(tmp_path / "run" / "manifest.csv")
    reg = {(r["id"], r["gjf_path"], r["sampler"], str(r["n_atoms"]), str(r["charge"]), str(r["mult"]),
            str(r["bs"]), r["job"]) for r in rows}
    assert manifest == reg
    n_before = con.execute("SELECT COUNT(*) FROM structures").fetchone()[0]
    con.close()

    assert cli.main(["generate", str(path)]) == 0
    con = registry.connect(db)
    assert con.execute("SELECT COUNT(*) FROM structures").fetchone()[0] == n_before
    con.close()

    assert cli.main(["resources", "--registry", str(db), "--nproc", "4", "--mem", "8GB",
                     "--sampler", "monomer"]) == 0
    with open(tmp_path / "run" / "manifest.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert all((r["nproc"], r["mem"]) == (("4", "8GB") if r["sampler"] == "monomer" else ("1", "1GB")) for r in rows)
    mono = next(r for r in rows if r["sampler"] == "monomer")
    assert open(mono["gjf_path"]).read().startswith("%nprocshared=4\n%mem=8GB\n")


def test_sampler_guards():
    import types
    import warnings
    from rmlipdata import sample
    ctx = types.SimpleNamespace(cfg={"charge": {"range": [-2, 2]}})
    assert sample._q_ok(ctx, -2) and sample._q_ok(ctx, 2) and not sample._q_ok(ctx, -3) and not sample._q_ok(ctx, 3)
    two = {"members": [("A", [0]), ("B", [1])], "positions": np.array([[0, 0, 0], [0, 0, 3.0]])}
    assert not sample._fragmented(two)
    two["positions"][1, 2] = 40.0                      # Coulomb-exploded pair
    assert sample._fragmented(two)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert sample._until_jump([0.0, 1.0, 2.0, -64.0, -63.0]) == 3   # drop points after a > 5 eV jump
    assert sample._until_jump([0.0, 1.0, 2.0]) == 3
    pos = np.array([[0, 0, 0], [0, 0, 5.0]])
    moved = sample._shift_to(pos, [0], [1], np.array([0, 0, 1.0]), 1.5, 5.0)
    assert np.isclose(moved[1, 2], 1.5, atol=1e-6) and (moved[0] == 0).all()
