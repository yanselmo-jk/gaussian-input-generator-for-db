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
        "bond_scan1": {"on": True, "n_points": 3, "max_bonds_per_cluster": 1, "bs": True, "max_structures": 3},
        "bond_scan2": {"on": True, "n_points": 4, "max_bonds_per_cluster": 1, "bs": True, "max_structures": 3},
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


def test_generate_end_to_end(tmp_path, capsys):
    cfg = dict(TINY, out=str(tmp_path / "run"))
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0

    db = tmp_path / "run" / "registry.sqlite"
    con = registry.connect(db)
    samplers = {r[0] for r in con.execute("SELECT DISTINCT sampler FROM structures")}
    assert {"monomer", "pair_scan", "cluster", "bond_scan1", "bond_scan2", "reaction", "hop"} <= samplers
    rows = con.execute("SELECT id, gjf_path, sampler, n_atoms, charge, mult, bs, job FROM structures "
                       "WHERE status='generated'").fetchall()
    assert rows and all(r["gjf_path"] and __import__("pathlib").Path(r["gjf_path"]).exists() for r in rows)
    manifest = _manifest_rows(tmp_path / "run" / "manifest.csv")
    reg = {(r["id"], r["gjf_path"], r["sampler"], str(r["n_atoms"]), str(r["charge"]), str(r["mult"]),
            str(r["bs"]), r["job"]) for r in rows}
    assert manifest == reg
    n_before = con.execute("SELECT COUNT(*) FROM structures").fetchone()[0]
    con.close()

    capsys.readouterr()
    assert cli.main(["generate", str(path)]) == 0           # rerun: xtb results come from <out>/xtb_cache
    err = capsys.readouterr().err
    assert sum(int(n) for n in __import__("re").findall(r"\(cached (\d+)\), ", err)) > 0
    con = registry.connect(db)
    assert con.execute("SELECT COUNT(*) FROM structures").fetchone()[0] == n_before
    con.close()

    cfg1 = dict(cfg, out=str(tmp_path / "run1"), xtb=dict(TINY["xtb"], n_max_process=1))
    path1 = tmp_path / "cfg1.yaml"
    path1.write_text(yaml.safe_dump(cfg1))
    assert cli.main(["generate", str(path1)]) == 0           # one stage at a time: same ids as 2 at a time
    con = registry.connect(tmp_path / "run1" / "registry.sqlite")
    ids1 = {r[0] for r in con.execute("SELECT id FROM structures")}
    con.close()
    con = registry.connect(db)
    assert ids1 == {r[0] for r in con.execute("SELECT id FROM structures")}
    con.close()

    mono_dir = tmp_path / "run" / "gjf" / "monomer"
    n_mono = len(list(mono_dir.glob("*.gjf")))
    assert cli.main(["submit", str(mono_dir), "--registry", str(db)]) == 0
    with open(tmp_path / "run" / "manifest.csv", newline="") as f:
        assert not any(r["sampler"] == "monomer" for r in csv.DictReader(f))   # submitted: out of the manifest
    capsys.readouterr()
    assert cli.main(["list", "--registry", str(db)]) == 0
    assert len(capsys.readouterr().out.splitlines()) == n_mono + 1
    assert cli.main(["generate", str(path)]) == 0           # a rerun does not put them back
    con = registry.connect(db)
    assert con.execute("SELECT COUNT(*) FROM structures WHERE status='submitted'").fetchone()[0] == n_mono
    con.close()
    assert cli.main(["submit", str(mono_dir), "--registry", str(db), "--undo"]) == 0
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
    frames = [[[0, 0, 0], [0, 0, z]] for z in (2.0, 3.5, 4.5, 3.0)]
    assert sample._until_apart(two, frames) == 2       # frames after the first > 4 A separation are dropped
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert sample._until_jump([0.0, 1.0, 2.0, -64.0, -63.0]) == 3   # drop points after a > 5 eV jump
    assert sample._until_jump([0.0, 1.0, 2.0]) == 3
    pos = np.array([[0, 0, 0], [0, 0, 5.0]])
    moved = sample._shift_to(pos, [0], [1], np.array([0, 0, 1.0]), 1.5, 5.0)
    assert np.isclose(moved[1, 2], 1.5, atol=1e-6) and (moved[0] == 0).all()


def test_xtb_scc_retry(tmp_path, monkeypatch):
    """SCC non-convergence -> one rerun with --etemp 1000; result carries etemp_K. CO2 + CO2 at q=-1 (from the
    stress matrix) fails at 300 K; whether it does is CPU-dependent, so the retry path is also forced below."""
    import subprocess
    import warnings
    from rmlipdata import xtb
    sym = ["O", "C", "O", "O", "C", "O"]
    pos = [[-1.7135, 0.3328, 1.2102], [-1.1690, -0.6660, 1.3293], [-0.6245, -1.6647, 1.4483],
           [-1.4815, -0.0896, -1.2429], [-1.5066, 0.4431, -2.2547], [-1.5316, 0.9758, -3.2665]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = xtb.run_xtb(sym, pos, -1, 1, tmp_path, mode="opt")
    assert r["ok"] and r.get("etemp_K", 1000) == 1000, r.get("err")
    calls, real = [], subprocess.run

    def fake(args, **kw):
        calls.append(list(args))
        if "--etemp" not in args:
            return subprocess.CompletedProcess(args, 1, "-1- scf: " + xtb.SCC_FAIL + "\n", "")
        return real(args, **kw)
    monkeypatch.setattr(xtb.subprocess, "run", fake)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = xtb.run_xtb(sym, pos, -1, 1, tmp_path, mode="sp")
    assert r["ok"] and r["etemp_K"] == 1000 and len(calls) == 2 and calls[1][-2:] == ["--etemp", "1000"]
    two = [[0, 1, 2], [3, 4, 5]]
    assert not build.fragmented(pos, two) and build.fragmented(np.array(pos) + np.array([[0, 0, 0]] * 3 + [[0, 0, -9.0]] * 3), two)


def test_view_gjf(tmp_path):
    from ase.io import read
    from rmlipdata import g16
    (tmp_path / "a.gjf").write_text(g16.deck(
        ["O", "H", "H"], [[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]], 0, 1, False, "sp", "run1-000001", 1, "1GB"))
    out = tmp_path / "out.xyz"
    assert cli.main(["view", str(tmp_path), "--save", str(out)]) == 0
    atoms = read(out, format="extxyz")
    assert len(atoms) == 3
    assert (atoms.info["charge"], atoms.info["spin"], atoms.info["id"]) == (0, 1, "run1-000001")
    assert cli.main(["view", "EC", "--save", str(tmp_path / "ec.xyz")]) == 0


def test_generate_warnings_log(tmp_path, monkeypatch, capsys):
    import warnings
    from rmlipdata import sample
    cfg = dict(TINY, out=str(tmp_path / "run"))
    cfg["samplers"] = {name: dict(sc, on=name == "monomer") for name, sc in TINY["samplers"].items()}
    real = sample.SAMPLERS["monomer"]

    def noisy(ctx, scfg, dq):
        warnings.warn("xtb scan energy jump -6.6 eV at point 10; later points dropped")
        warnings.warn("xtb scan energy jump +5.1 eV at point 3; later points dropped")
        yield from real(ctx, scfg, dq)

    monkeypatch.setitem(sample.SAMPLERS, "monomer", noisy)
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0
    log = tmp_path / "run" / "warnings.log"
    assert len(log.read_text().splitlines()) >= 2
    assert "2  xtb scan energy jump # eV at point #; later points dropped" in capsys.readouterr().err


def test_generate_progress(tmp_path, capsys):
    cfg = dict(TINY, out=str(tmp_path / "run"))
    cfg["samplers"] = {name: dict(sc, on=name in ("monomer", "cluster"))
                       for name, sc in TINY["samplers"].items()}
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0
    err = capsys.readouterr().err
    assert "stages (delta_q x sampler)" in err
    assert "[1/" in err
    assert " done: " in err
    assert cli.main(["generate", str(path), "--dry"]) == 0
    assert "monomer: " in capsys.readouterr().out
