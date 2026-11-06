"""Crash-safe deck writes, vertical-charge guards, parent ids and the JSON xtb cache."""
import csv
import json
import warnings
from pathlib import Path

import numpy as np
import pytest
import yaml

from rmlipdata import cli, registry, sample, xtb

H2O = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])


def _cfg(tmp_path, **charge):
    return {
        "out": str(tmp_path / "run"),
        "seed": 0,
        "species": {"cations": ["Li"], "anions": [], "solvents": ["EC"]},
        "charge": {"delta_q": [0], "mult_offsets": [0], "range": [-2, 2], "vertical_fraction": 0.0, **charge},
        "g16": {"nproc": 1, "mem": "1GB"},
        "limits": {"max_atoms": 60, "dedup_A": 0.0},
        "xtb": {"threads": 2, "timeout_s": 120, "n_max_process": 1},
        "cluster": {},
        "samplers": {"fake": {"on": True}},
    }


def _run(tmp_path, monkeypatch, recs, **charge):
    def fake(ctx, scfg, dq):
        yield from recs
    monkeypatch.setitem(sample.SAMPLERS, "fake", fake)
    cfg = _cfg(tmp_path, **charge)
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))
    assert cli.main(["generate", str(path)]) == 0
    return registry.connect(tmp_path / "run" / "registry.sqlite"), cfg


def _rec(symbols, positions, charge, job="sp", params=None):
    return dict(symbols=list(symbols), positions=np.asarray(positions, float), charge=int(charge), sampler="fake",
                composition="H2O1", job=job, bs=False, params=params or {})


def test_deck_write_crash_recovery(tmp_path, monkeypatch):
    """A crash while writing the deck leaves no half-row; a legacy NULL gjf_path row is repaired on rerun."""
    rec = _rec(["O", "H", "H"], H2O, 0)

    def fake(ctx, scfg, dq):
        yield rec
    monkeypatch.setitem(sample.SAMPLERS, "fake", fake)
    cfg = _cfg(tmp_path)
    path = tmp_path / "cfg.yaml"
    path.write_text(yaml.safe_dump(cfg))

    real = cli.g16.deck

    def flaky(*a, **k):
        monkeypatch.setattr(cli.g16, "deck", real)  # fail once, then behave
        raise RuntimeError("simulated crash during deck write")
    monkeypatch.setattr(cli.g16, "deck", flaky)
    with pytest.raises(RuntimeError):
        cli.main(["generate", str(path)])

    db = tmp_path / "run" / "registry.sqlite"
    con = registry.connect(db)
    sid, _ = registry.add(con, dict(rec, mult=1))  # the old crash order: row written, deck path still NULL
    con.commit()
    assert registry.gjf_path(con, sid) is None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cli.manifest(con, tmp_path / "run")
    assert any("gjf_path NULL" in str(w.message) for w in caught)
    con.close()

    assert cli.main(["generate", str(path)]) == 0
    con = registry.connect(db)
    row = con.execute("SELECT gjf_path FROM structures WHERE id=?", (sid,)).fetchone()
    con.close()
    assert row["gjf_path"] and Path(row["gjf_path"]).is_file()
    with open(tmp_path / "run" / "manifest.csv", newline="") as f:
        assert sid in {r["id"] for r in csv.DictReader(f)}


def test_vertical_never_beyond_formal_band(tmp_path, monkeypatch):
    """A delta_q=+1 record (formal charge 0) may not get a q=+2 vertical deck."""
    rec = _rec(["O", "H", "H"], H2O, 1)
    con, _ = _run(tmp_path, monkeypatch, [rec], delta_q=[1], vertical_fraction=1.0)
    rows = con.execute("SELECT charge, params FROM structures").fetchall()
    con.close()
    assert {r["charge"] for r in rows} == {1}


def test_vertical_skips_charge_independent_geometries(tmp_path, monkeypatch):
    """Bare cations and rcov-sum cation pairs are charge-independent: no vertical decks, so no id collisions."""
    recs = [_rec(["Li"], np.zeros((1, 3)), 1, params={"kind": "cation"}),
            _rec(["Li", "Li"], [[0.0, 0.0, 0.0], [2.56, 0.0, 0.0]], 2, params={"contact": "rcov_sum"})]
    con, _ = _run(tmp_path, monkeypatch, recs, vertical_fraction=1.0)
    rows = con.execute("SELECT charge, params FROM structures").fetchall()
    con.close()
    assert not [r for r in rows if json.loads(r["params"]).get("vertical_dq")]


def test_vertical_parent_is_lowest_multiplicity(tmp_path, monkeypatch):
    """mult_offsets [2, 0] emits the triplet first, but parent_id points at the singlet base deck."""
    rec = _rec(["O", "H", "H"], H2O, 0)
    con, _ = _run(tmp_path, monkeypatch, [rec], mult_offsets=[2, 0], vertical_fraction=1.0,
                  high_spin_fraction=1.0)
    rows = con.execute("SELECT id, parent_id, charge, mult, params FROM structures").fetchall()
    con.close()
    base = next(r for r in rows if r["charge"] == 0 and r["mult"] == 1)
    vert = [r for r in rows if json.loads(r["params"]).get("vertical_dq")]
    assert vert and all(r["parent_id"] == base["id"] for r in vert)


def test_xtb_cache_json_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(xtb, "CACHE", tmp_path / "cache")
    monkeypatch.setattr(xtb, "HITS", 0)
    r1 = xtb.run_xtb(["O", "H", "H"], H2O, 0, 0, tmp_path, mode="opt")
    assert r1["ok"], r1.get("err")
    files = list((tmp_path / "cache").rglob("*.*"))
    assert files and all(p.suffix == ".json" for p in files)
    r2 = xtb.run_xtb(["O", "H", "H"], H2O, 0, 0, tmp_path, mode="opt")
    assert xtb.HITS == 1
    assert np.array_equal(np.asarray(r1["positions"]), np.asarray(r2["positions"]))
    assert r1["E_eV"] == r2["E_eV"]


def test_gjf_info_without_route_line(tmp_path):
    p = tmp_path / "bad.gjf"
    p.write_text("%nprocshared=1\n%mem=1GB\n\ntitle\n\n0 1\nH 0.0 0.0 0.0\n\n\n")
    with pytest.raises(ValueError, match="no route line"):
        cli._gjf_info(p)


def test_vertical_reduction_stops_at_minus_one(tmp_path, monkeypatch):
    """A Li-containing record already at q=-1 gets no vertical q=-2 deck (unbound gas-phase dianion)."""
    li_f2 = [[0.0, 0.0, 0.0], [1.6, 0.0, 0.0], [-1.6, 0.0, 0.0]]
    con, _ = _run(tmp_path, monkeypatch, [_rec(["Li", "F", "F"], li_f2, -1)], vertical_fraction=1.0)
    assert sorted({r[0] for r in con.execute("SELECT charge FROM structures")}) == [-1, 0]


def test_stage_record_wins_over_vertical(tmp_path, monkeypatch):
    """When a vertical deck and a stage record are the same deck, the stage record's provenance is kept."""
    li_f = [[0.0, 0.0, 0.0], [1.6, 0.0, 0.0]]
    recs = [_rec(["Li", "F"], li_f, 0, params={"src": "a"}), _rec(["Li", "F"], li_f, 1, params={"src": "b"})]
    con, _ = _run(tmp_path, monkeypatch, recs, vertical_fraction=1.0)
    row = con.execute("SELECT parent_id, params FROM structures WHERE charge=1").fetchone()
    assert row[0] is None and json.loads(row[1])["src"] == "b" and "vertical_dq" not in json.loads(row[1])
