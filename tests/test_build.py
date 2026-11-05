import csv
import json
import re

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
        "cluster": {"on": True, "n_distort": 1, "rattle_A": [0.05, 0.15]},
        "cluster_md": {"on": True, "n_solvent": [2], "with_anion": False, "T_K": [400], "time_ps": 0.6,
                       "equil_ps": 0.1, "dump_fs": 50, "n_frames": 2, "max_clusters": 1},
        "bond_scan1": {"on": True, "n_points": 3, "max_bonds_per_cluster": 1, "bs": True, "max_structures": 3},
        "bond_scan2": {"on": True, "n_points": 3, "max_bonds_per_cluster": 1, "bs": True, "max_structures": 3},
        "reaction": {"on": True, "templates": ["add", "substitute", "transfer"], "n_points": 3, "k_peak": 1,
                     "max_per_cluster": 1, "bs": True},
        "reaction2": {"on": True, "templates": ["ring_open_M"],
                      "n_points": 4, "max_per_cluster": 4, "bs": True},
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


def test_reaction_candidates_li_ec(tmp_path):
    """Li+EC: ring_open_M and intra_shift candidates; every pair has exactly one range."""
    import types
    from rmlipdata import sample
    cat = species.catalog()
    combo = ["Li", "EC"]
    cl = build.assemble(combo, cat, tmp_path, 1, np.random.default_rng(0))
    ctx = types.SimpleNamespace(cat=cat, workdir=tmp_path, cfg={}, rng=np.random.default_rng(0))
    cand = sample._reaction_candidates(ctx, cl, combo)
    assert {"ring_open_M", "intra_shift"} <= {c[0] for c in cand}
    assert all(len(pairs) == len(ranges) for _, pairs, ranges in cand)


def test_reaction_candidates_radicals(tmp_path):
    """Two radical members give radical_couple; a radical taking a neighbour's H gives h_abstract."""
    import types
    from rmlipdata import sample
    cat = species.catalog()
    combo = ["CH3", "CH3O"]
    cl = build.assemble(combo, cat, tmp_path, 0, np.random.default_rng(0))
    ctx = types.SimpleNamespace(cat=cat, workdir=tmp_path, cfg={}, rng=np.random.default_rng(0))
    assert {"radical_couple", "h_abstract"} <= {c[0] for c in sample._reaction_candidates(ctx, cl, combo)}


def test_charge_radicals(tmp_path):
    """delta_q != 0: the xtb charge difference marks EC atoms as radical sites; delta_q 0 adds none, and an
    extra radical on EC with an H on another member gives an h_abstract candidate."""
    import types
    from rmlipdata import sample
    cat = species.catalog()
    combo = ["Li", "EC"]
    cl = build.assemble(combo, cat, tmp_path, 0, np.random.default_rng(0))  # Li+EC at q=0 = delta_q -1
    ctx = types.SimpleNamespace(cat=cat, workdir=tmp_path, cfg={}, rng=np.random.default_rng(0))
    extra = sample._charge_radicals(ctx, cl, 0, -1)
    ec = next(idxs for name, idxs in cl["members"] if name == "EC")
    assert extra and all(i in ec and cl["symbols"][i] != "H" for i in extra)
    assert sample._charge_radicals(ctx, cl, 0, 0) == []

    combo3 = ["Li", "EC", "H2O"]
    cl3 = build.assemble(combo3, cat, tmp_path, 0, np.random.default_rng(0))
    extra3 = sample._charge_radicals(ctx, cl3, 0, -1)
    assert extra3
    assert "h_abstract" in {c[0] for c in sample._reaction_candidates(ctx, cl3, combo3, extra3)}


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
    assert {"monomer", "pair_scan", "cluster", "cluster_md", "bond_scan1", "bond_scan2", "reaction", "reaction2",
            "hop"} <= samplers
    rows = con.execute("SELECT id, gjf_path, sampler, n_atoms, charge, mult, bs, job FROM structures "
                       "WHERE status='generated'").fetchall()
    assert rows and all(r["gjf_path"] and __import__("pathlib").Path(r["gjf_path"]).exists() for r in rows)
    for r in rows:  # deck folders by job type: gjf/<job>/<sampler>/<id>.gjf
        p = __import__("pathlib").Path(r["gjf_path"])
        assert (p.parent.name, p.parent.parent.name, p.parent.parent.parent.name) == \
            (r["sampler"], r["job"], "gjf")
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
    for sampler, key, allowed in (("hop", "path_dev_A", None), ("hop", "ref3", None),
                                  ("bond_scan2", "endpoint", {"xtb", "rigid"})):
        rows = con.execute("SELECT params FROM structures WHERE sampler=?", (sampler,)).fetchall()
        params = [json.loads(r["params"]) for r in rows]
        if sampler == "hop":
            params = [p for p in params if p.get("kind") == "hop"]
        assert params and all(key in p for p in params), (sampler, key)
        if allowed:
            assert all(p[key] in allowed for p in params)
    for sampler, kind in (("cluster", "dilate"), ("hop", "exchange_out")):
        rows = con.execute("SELECT params FROM structures WHERE sampler=?", (sampler,)).fetchall()
        assert any(json.loads(r["params"]).get("kind") == kind for r in rows), (sampler, kind)
    rows = con.execute("SELECT params FROM structures WHERE sampler='bond_scan2'").fetchall()
    assert any(json.loads(r["params"]).get("far") for r in rows)
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

    mono_dirs = sorted((tmp_path / "run" / "gjf").glob("*/monomer"))   # gjf/sp/monomer and gjf/opt/monomer
    n_mono = sum(len(list(d.glob("*.gjf"))) for d in mono_dirs)
    assert cli.main(["submit", *map(str, mono_dirs), "--registry", str(db)]) == 0
    with open(tmp_path / "run" / "manifest.csv", newline="") as f:
        assert not any(r["sampler"] == "monomer" for r in csv.DictReader(f))   # submitted: out of the manifest
    capsys.readouterr()
    assert cli.main(["list", "--registry", str(db)]) == 0
    assert len(capsys.readouterr().out.splitlines()) == n_mono + 1
    assert cli.main(["generate", str(path)]) == 0           # a rerun does not put them back
    con = registry.connect(db)
    assert con.execute("SELECT COUNT(*) FROM structures WHERE status='submitted'").fetchone()[0] == n_mono
    con.close()
    assert cli.main(["submit", *map(str, mono_dirs), "--registry", str(db), "--undo"]) == 0
    assert cli.main(["resources", "--registry", str(db), "--nproc", "4", "--mem", "8GB",
                     "--sampler", "monomer"]) == 0
    with open(tmp_path / "run" / "manifest.csv", newline="") as f:
        rows = list(csv.DictReader(f))
    assert all((r["nproc"], r["mem"]) == (("4", "8GB") if r["sampler"] == "monomer" else ("1", "1GB")) for r in rows)
    mono = next(r for r in rows if r["sampler"] == "monomer")
    assert open(mono["gjf_path"]).read().startswith("%nprocshared=4\n%mem=8GB\n")


def test_fraction_controls(tmp_path, monkeypatch):
    """redox_fraction 0 drops every delta_q != 0 record; high_spin_fraction 0 drops non-BS decks above the
    lowest multiplicity (the broken-symmetry singlet stays)."""
    from ase.io import read
    from rmlipdata import g16, sample
    syms = ["O", "H", "H"]
    base = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])

    def fake(ctx, scfg, dq):
        yield dict(symbols=syms, positions=base, charge=int(dq), sampler="fake", composition="H2O1",
                   job="sp", bs=False, params={})
        if int(dq) == 0:  # same formal charge, broken-symmetry record: its singlet deck must survive
            yield dict(symbols=syms, positions=base + [[0.3, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                       charge=0, sampler="fake",
                       composition="H2O1", job="sp", bs=True, params={})

    monkeypatch.setitem(sample.SAMPLERS, "fake", fake)
    cfg = dict(TINY, out=str(tmp_path / "run"),
               charge=dict(TINY["charge"], delta_q=[0, -1], redox_fraction=0.0, high_spin_fraction=0.0),
               samplers={"fake": {"on": True}})
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0

    con = registry.connect(tmp_path / "run" / "registry.sqlite")
    rows = con.execute("SELECT composition, charge, mult, bs, gjf_path FROM structures").fetchall()
    con.close()
    assert rows and all("/gjf/sp/fake/" in r["gjf_path"] for r in rows)
    cat = species.catalog()
    for r in rows:  # only delta_q 0 records: charge equals the formal charge of the composition
        formal = sum(cat[m[1]].charge * int(m[2]) for m in
                     (re.fullmatch(r"(.+?)_?(\d+)", tok) for tok in r["composition"].split()))
        assert r["charge"] == formal, r
    for r in rows:  # no non-BS deck above the lowest multiplicity
        if r["bs"]:
            continue
        atoms = read(r["gjf_path"], format="gaussian-in")
        mult = 1 if g16.n_electrons(atoms.get_chemical_symbols(), r["charge"]) % 2 == 0 else 2
        assert r["mult"] == mult, r
    assert any(r["bs"] and r["mult"] == 1 for r in rows)


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


def test_cluster_dilation(monkeypatch):
    """Dilation scales each member COM from the cluster COM and keeps the member geometry rigid."""
    from types import SimpleNamespace
    from rmlipdata import sample
    pos = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0],
                    [-3.0, 0.8, 0.0], [-2.1, 0.8, 0.0], [-3.8, 0.8, 0.0]])
    cl = {"symbols": ["C", "O", "C", "C", "O"], "positions": pos,
          "members": [("EC", [0, 1]), ("CO2", [2, 3, 4])], "composition": "EC1 CO2_1"}
    monkeypatch.setattr(sample, "_combos", lambda ctx: [("EC", "CO2")])
    monkeypatch.setattr(sample, "_cluster", lambda ctx, combo, dq: (0, cl))
    ctx = SimpleNamespace(cfg={}, cat=species.catalog(), workdir=None, rng=np.random.default_rng(0))
    recs = list(sample.cluster(ctx, {"n_distort": 0, "dilate": [0.9, 1.5]}, 0))
    assert [r["params"]["kind"] for r in recs] == ["opt", "dilate", "dilate"]
    com = pos.mean(axis=0)
    for rec in recs[1:]:
        f, moved = rec["params"]["factor"], rec["positions"]
        for _, idxs in cl["members"]:
            assert np.isclose(np.linalg.norm(sample._com(moved, idxs) - com),
                              f * np.linalg.norm(sample._com(pos, idxs) - com))
            assert np.allclose(np.linalg.norm(moved[idxs][:, None] - moved[idxs][None], axis=-1),
                               np.linalg.norm(pos[idxs][:, None] - pos[idxs][None], axis=-1))


def test_maxmin_picks_diverse_frames():
    """Greedy max-min RMSD: frame 0, then the frame farthest from it; rigid translations/rotations count as equal."""
    from rmlipdata import sample
    base = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.5, 0.0]])
    frames = [base + [0.0, 0.0, float(i)] for i in range(4)] + [base * 2.0]   # 0-3 identical, 4 expanded
    assert sample._maxmin(frames, 2) == [0, 4]
    assert sample._maxmin(frames, 1) == [0]
    assert sorted(sample._maxmin(frames, 10)) == list(range(5))
    assert sample._maxmin([], 3) == [] and sample._maxmin(frames, 0) == []
    scaled = [base * s for s in (1.0, 2.0, 1.9, 1.5)]   # 1.9 is close to 2.0, so 1.5 is the third pick
    assert sample._maxmin(scaled, 3) == [0, 1, 3]


def test_bonds_ignore_li_na():
    from rmlipdata import sample
    pos = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [0.0, 0.0, 1.9], [0.0, 0.0, 10.0]])
    assert sample._bonds(["C", "O", "Li", "Na"], pos) == {(0, 1)}


def test_pair_scan_bonded_members(tmp_path, monkeypatch):
    """A1: an assembled pair whose atoms ended up covalently bonded must not be distance-scanned."""
    import warnings
    from types import SimpleNamespace
    from rmlipdata import sample
    cl = {"symbols": ["C", "O"], "positions": np.array([[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]]),
          "members": [("EC", [0]), ("CO2", [1])], "composition": "EC1 CO2_1"}
    monkeypatch.setattr(build, "assemble", lambda *a, **k: cl)
    cfg = {"charge": {"range": [-2, 2]},
           "species": {"cations": [], "anions": [], "solvents": ["EC", "CO2"]}}
    ctx = SimpleNamespace(cfg=cfg, cat=species.catalog(), workdir=tmp_path, rng=np.random.default_rng(0))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert list(sample.pair_scan(ctx, {"n_orient": 1, "scales": [1.0]}, 0)) == []
    assert any("bonded the members" in str(w.message) for w in caught)


def test_pair_scan_cations(tmp_path, monkeypatch):
    """A2/A3: cation pairs are placed by covalent-radius sum at delta_q 0 and skipped at delta_q +1."""
    from types import SimpleNamespace
    from rmlipdata import sample
    cat = species.catalog()
    cfg = {"charge": {"range": [-2, 2]},
           "species": {"cations": ["Li", "Na"], "anions": [], "solvents": []}}

    def boom(*a, **k):
        raise AssertionError("a cation-only pair must not call assemble/xtb")

    monkeypatch.setattr(build, "assemble", boom)
    monkeypatch.setattr(sample.xtb, "run_xtb", boom)
    scfg = {"scales": [0.75, 1.0], "n_orient": 1}
    ctx = SimpleNamespace(cfg=cfg, cat=cat, workdir=tmp_path, rng=np.random.default_rng(0))
    recs = list(sample.pair_scan(ctx, scfg, 0))
    assert len(recs) == 6
    first = {}
    for r in recs:
        p = r["params"]
        assert p["contact"] == "rcov_sum"
        first.setdefault(tuple(p["pair"]), p["d_min_A"])
    assert set(first) == {("Li", "Li"), ("Li", "Na"), ("Na", "Na")}
    for (a, b), d in first.items():
        assert np.isclose(d, 0.75 * (sample._cov([a, b], 0) + sample._cov([a, b], 1)), atol=1e-6)
    ctx = SimpleNamespace(cfg=cfg, cat=cat, workdir=tmp_path, rng=np.random.default_rng(0))
    assert list(sample.pair_scan(ctx, scfg, 1)) == []
    assert sample._cluster(SimpleNamespace(cfg=cfg, cat=cat), ("Li", "Na"), 1) == (3, None)


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


def test_pick_balances_templates():
    """Round-robin over templates in a random template order: with k=1 every template can be picked."""
    from rmlipdata import sample
    cand = [("add", i) for i in range(10)] + [("intra_shift", i) for i in range(50)] + [("transfer", 0)]
    first = {sample._pick(cand, 1, np.random.default_rng(s))[0][0] for s in range(30)}
    assert first == {"add", "intra_shift", "transfer"}
    three = sample._pick(cand, 3, np.random.default_rng(0))
    assert sorted(c[0] for c in three) == ["add", "intra_shift", "transfer"]


def test_reaction2_endpoint_check():
    """reaction2 keeps an endpoint only if it made exactly the template's bond changes."""
    from rmlipdata import sample
    sym = ["C", "C", "H", "O"]
    start = np.array([[0, 0, 0], [3.0, 0, 0], [-1.09, 0, 0], [6.0, 0, 0]], float)
    b0 = sample._bonds(sym, start)
    pairs, ranges = [(0, 1)], [(3.0, 1.6)]                      # radical_couple C0-C1 at 1.6 A
    assert sample._endpoint_check(sym, b0, start, pairs, ranges)[0] == "missed a target distance"
    coupled = start.copy()
    coupled[1, 0] = 1.6
    assert sample._endpoint_check(sym, b0, coupled, pairs, ranges)[0] == ""
    moved = coupled.copy()
    moved[3] = [2.9, 0, 0]                                       # O now bonded to C1 too: rearranged
    assert sample._endpoint_check(sym, b0, moved, pairs, ranges)[0] == "rearranged"
    assert sample._endpoint_check(sym, sample._bonds(sym, coupled), coupled, pairs, ranges)[0] == \
        "already in the product state"


def test_components_flag_detached_atom():
    from rmlipdata import build, sample
    sym = ["Li", "O", "H", "H", "H"]
    pos = [[0, 0, 0], [1.9, 0, 0], [2.5, 0.75, 0], [2.5, -0.75, 0], [9.0, 0, 0]]
    comps = sample._components(sym, pos)
    assert sorted(map(sorted, comps)) == [[0], [1, 2, 3], [4]]
    assert build.fragmented(pos, comps) and not build.fragmented(pos[:4], sample._components(sym[:4], pos[:4]))


def _fake_cfg(tmp_path, **charge):
    return dict(TINY, out=str(tmp_path / "run"), samplers={"fake": {"on": True}},
                charge=dict(TINY["charge"], **charge))


def test_vertical_charge_pairs(tmp_path, monkeypatch, capsys):
    """vertical_fraction 1: a Li sp record adds q-1 and q+1 decks with the base deck as parent_id; a record
    without Li/Na gets no q-1 deck and an opt record gets no vertical decks."""
    from ase.io import read
    from rmlipdata import sample
    li = np.array([[0.0, 0.0, 0.0], [1.9, 0.0, 0.0], [2.5, 0.7, 0.0], [2.5, -0.7, 0.0]])
    h2o = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])

    def fake(ctx, scfg, dq):
        yield dict(symbols=["Li", "O", "H", "H"], positions=li, charge=0, sampler="fake",
                   composition="Li1 H2O1", job="sp", bs=False, params={"kind": "sp"})
        yield dict(symbols=["O", "H", "H"], positions=h2o, charge=0, sampler="fake",
                   composition="H2O1", job="sp", bs=False, params={})
        yield dict(symbols=["Li", "O", "H", "H"], positions=li + [0.0, 0.0, 1.0], charge=0, sampler="fake",
                   composition="Li1 H2O1", job="opt", bs=False, params={})

    monkeypatch.setitem(sample.SAMPLERS, "fake", fake)
    cfg = dict(_fake_cfg(tmp_path, vertical_fraction=1.0, high_spin_fraction=0.0),
               limits=dict(TINY["limits"], dedup_A=0.0))
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0

    con = registry.connect(tmp_path / "run" / "registry.sqlite")
    rows = con.execute("SELECT id, parent_id, charge, mult, job, params, composition, gjf_path FROM structures"
                       ).fetchall()
    con.close()
    sp_li = [r for r in rows if r["job"] == "sp" and r["composition"] == "Li1 H2O1"]
    assert {r["charge"] for r in sp_li} == {-1, 0, 1}
    vert = [(r, json.loads(r["params"])["vertical_dq"]) for r in rows
            if json.loads(r["params"]).get("vertical_dq")]
    assert sorted(v for _, v in vert) == [-1, 1, 1]          # Li: +-1; H2O: +1 only
    base = {r["composition"]: r for r in rows if r["job"] == "sp" and r["charge"] == 0}
    assert all(r["parent_id"] == base[r["composition"]]["id"] for r, _ in vert)
    assert all("Li" in r["composition"] or "Na" in r["composition"] for r, v in vert if v == -1)
    assert not [(r, v) for r, v in vert if r["composition"] == "H2O1" and v == -1]   # no cation
    assert not [(r, v) for r, v in vert if r["job"] == "opt"]                # vertical pairs only for sp
    x = read(base["Li1 H2O1"]["gjf_path"], format="gaussian-in").get_positions()
    assert all(np.allclose(read(r["gjf_path"], format="gaussian-in").get_positions(), x)
               for r, _ in vert if r["composition"] == "Li1 H2O1")


def test_dedup_near_duplicates(tmp_path, monkeypatch, capsys):
    """A rotated copy is the same geometry and is dropped; a 0.3 A displacement is a different point."""
    from rmlipdata import sample
    base = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])
    rot = base @ np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    moved = base.copy()
    moved[1, 0] += 0.3

    def fake(ctx, scfg, dq):
        for pos in (base, rot, moved):
            yield dict(symbols=["O", "H", "H"], positions=pos, charge=0, sampler="fake", composition="H2O1",
                       job="sp", bs=False, params={})

    monkeypatch.setitem(sample.SAMPLERS, "fake", fake)
    cfg = dict(_fake_cfg(tmp_path, vertical_fraction=0.0, high_spin_fraction=0.0),
               limits=dict(TINY["limits"], dedup_A=0.1))
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0
    assert "1 near-duplicates dropped" in capsys.readouterr().out

    con = registry.connect(tmp_path / "run" / "registry.sqlite")
    n = con.execute("SELECT COUNT(*) FROM structures").fetchone()[0]
    con.close()
    assert n == 2


def test_per_sampler_high_spin_fraction(tmp_path, monkeypatch):
    """samplers.<name>.high_spin_fraction overrides the global value for that sampler only."""
    from rmlipdata import sample
    pos = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])

    def make(name):
        def gen(ctx, scfg, dq):
            yield dict(symbols=["O", "H", "H"], positions=pos, charge=0, sampler=name, composition="H2O1",
                       job="sp", bs=False, params={})
        return gen

    monkeypatch.setitem(sample.SAMPLERS, "fakelo", make("fakelo"))
    monkeypatch.setitem(sample.SAMPLERS, "fakehi", make("fakehi"))
    cfg = dict(TINY, out=str(tmp_path / "run"),
               charge=dict(TINY["charge"], vertical_fraction=0.0, high_spin_fraction=1.0),
               limits=dict(TINY["limits"], dedup_A=0.0),
               samplers={"fakelo": {"on": True, "high_spin_fraction": 0.0}, "fakehi": {"on": True}})
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0

    con = registry.connect(tmp_path / "run" / "registry.sqlite")
    got = {name: sorted(r[0] for r in con.execute("SELECT mult FROM structures WHERE sampler=?", (name,)))
           for name in ("fakelo", "fakehi")}
    con.close()
    assert got == {"fakelo": [1], "fakehi": [1, 3]}


def test_q_ok_unbound_electron():
    """A reduction needs Li/Na or a neutral radical; otherwise the extra electron is unbound (silent skip)."""
    import types
    from rmlipdata import sample
    cat = species.catalog()
    ctx = types.SimpleNamespace(cat=cat, cfg={"charge": {"range": [-2, 2]}})

    def ok(names):
        return sample._q_ok(ctx, sample._formal(ctx, names) - 1, names, -1)

    assert not ok(["EC"])
    assert ok(["Li", "EC"])
    assert ok(["CH3O"])
    assert not ok(["PF6"])
    assert ok(["LiF"])


def test_too_close_absolute_floor():
    from rmlipdata import sample
    assert sample._too_close(["H", "H"], np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.5]]), 0.7)
    assert not sample._too_close(["H", "H"], np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.7]]), 0.7)


def test_reaction_hop_skip_cluster_minimum(tmp_path):
    """The relaxed cluster is already the cluster opt deck: no reaction index 0, no hop t=0, no exchange k=0."""
    cfg = dict(TINY, out=str(tmp_path / "run"))
    cfg["samplers"] = {name: dict(sc, on=name in ("reaction", "hop")) for name, sc in TINY["samplers"].items()}
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0
    con = registry.connect(tmp_path / "run" / "registry.sqlite")
    rows = con.execute("SELECT sampler, params FROM structures WHERE sampler IN ('reaction', 'hop')").fetchall()
    con.close()
    assert {s for s, _ in rows} == {"reaction", "hop"}
    for sampler, raw in rows:
        p = json.loads(raw)
        if sampler == "reaction":
            assert p["index"] != 0
        elif p["kind"] == "hop":
            assert p["t"] > 0
        else:
            assert p["index"] != 0


def test_md_compositions_dict_and_max_atoms():
    from rmlipdata import sample
    cat = species.catalog()
    cfg = {"species": {"cations": ["Li", "Na"], "anions": [], "solvents": ["EC"]},
           "limits": {"max_atoms": 60}}
    scfg = {"n_solvent": {"Li": [2, 3], "Na": [4, 5]}, "with_anion": False}
    out = sample.md_compositions(cfg, cat, scfg)
    assert len([n for n in out if n[0] == "Li"]) == 2
    assert len([n for n in out if n[0] == "Na"]) == 2
    limited = sample.md_compositions(cfg, cat, dict(scfg, max_atoms=25))
    assert 0 < len(limited) < len(out)
    assert all(sum(build.n_atoms(cat[x]) for x in n) <= 25 for n in limited)


def test_monomer_rng_per_species(tmp_path, monkeypatch):
    """A failed hess on the first species must not shift the next species' normal-mode/distortion draws."""
    from types import SimpleNamespace
    from rmlipdata import sample
    cfg = {"seed": 0, "charge": {"range": [-2, 2]},
           "species": {"cations": [], "anions": [], "solvents": ["EC", "H2O"]}}
    real = sample.xtb.run_xtb

    def fake(*a, mode="sp", **kw):
        if mode == "hess" and "C" in a[0]:
            return {"ok": False, "err": "forced hess failure"}
        return real(*a, mode=mode, **kw)

    monkeypatch.setattr(sample.xtb, "run_xtb", fake)
    scfg = {"n_normal_mode": 1, "n_distort": 2, "T_K": 500}

    def run(names):
        c = dict(cfg, species=dict(cfg["species"], solvents=names))
        ctx = SimpleNamespace(cfg=c, cat=species.catalog(), workdir=tmp_path, rng=np.random.default_rng(0))
        return [r for r in sample.monomer(ctx, scfg, 0) if r["composition"] == build.composition(["H2O"])]

    both, alone = run(["EC", "H2O"]), run(["H2O"])
    assert both and len(both) == len(alone)
    for a, b in zip(both, alone):
        assert np.allclose(a["positions"], b["positions"])


def test_per_combo_rng_isolated(tmp_path, monkeypatch):
    """A failed xtb opt on the first combo must not shift a later combo's rattle draws."""
    from types import SimpleNamespace
    from rmlipdata import sample
    cat = species.catalog()
    combo1, combo2 = ("Li", "PF6"), ("Li", "EC")

    def cl_for(combo, bound):
        if not bound:  # a Coulomb-exploded cluster: _cluster drops it
            return {"symbols": ["Li", "F"], "positions": np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 40.0]]),
                    "members": [("Li", [0]), ("F", [1])]}
        if combo == combo1:
            pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.9], [0.0, 1.2, 2.3]])
            return {"symbols": ["Li", "P", "F"], "positions": pos,
                    "members": [("Li", [0]), ("PF6", [1, 2])]}
        pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.5], [0.0, 1.3, 1.1]])
        return {"symbols": ["Li", "C", "O"], "positions": pos,
                "members": [("Li", [0]), ("EC", [1, 2])]}

    def fake_assemble(combo, cat_, workdir, q, rng):
        rng.normal(size=1)  # stand-in for build.assemble's own draws
        if combo != combo1:
            return cl_for(combo, True)
        r = sample.xtb.run_xtb(["Li", "F"], np.zeros((2, 3)), q, 0, workdir, mode="opt")
        return cl_for(combo, bool(r.get("ok")))

    monkeypatch.setattr(sample, "_combos", lambda ctx: [combo1, combo2])
    monkeypatch.setattr(build, "assemble", fake_assemble)
    real = sample.xtb.run_xtb

    def probe(ok):
        def fake(*a, mode="sp", **kw):
            if mode == "opt" and a[0] == ["Li", "F"]:
                return {"ok": ok, "err": "" if ok else "forced opt failure"}
            return real(*a, mode=mode, **kw)
        return fake

    cfg = {"seed": 0, "charge": {"range": [-2, 2]}}
    scfg = {"n_distort": 2, "rattle_A": 0.1, "dilate": []}

    def run(ok):
        with monkeypatch.context() as m:
            m.setattr(sample.xtb, "run_xtb", probe(ok))
            ctx = SimpleNamespace(cfg=cfg, cat=cat, workdir=tmp_path, clusters={},
                                  rng=np.random.default_rng(0))
            return [r for r in sample.cluster(ctx, scfg, 0)
                    if r["composition"] == build.composition(list(combo2))]

    both, alone = run(False), run(True)
    assert len(both) == len(alone) == 3  # opt + 2 rattles for the second combo
    for a, b in zip(both, alone):
        assert np.allclose(a["positions"], b["positions"])


def test_q_ok_li_dianion():
    """Li + 2 PF6 at q=-2 is Li0 plus two PF6-: an unbound gas-phase dianion."""
    import types
    from rmlipdata import sample
    ctx = types.SimpleNamespace(cat=species.catalog(), cfg={"charge": {"range": [-2, 2]}})
    assert not sample._q_ok(ctx, -2, ["Li", "PF6", "PF6"], -1)
    assert sample._q_ok(ctx, -1, ["Li", "PF6"], -1)
