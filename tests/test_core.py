from pathlib import Path

import numpy as np
from ase.io import read, write

from rmlipdata import g16, registry
from rmlipdata.parse import HA_EV, log_atoms

FIX = Path(__file__).parent / "fixtures" / "li_ec3vc_red_trim.log"


def test_parse_opt_log(tmp_path):
    # trimmed copy of 2-2_sei_modeling .../li_ec3vc_red.log: UwB97XD opt, doublet, frames 1,2,72,73 + final pop
    fr = log_atoms(FIX)
    assert len(fr) == 4 and all(len(a) == 39 for a in fr)
    assert np.isclose(fr[0].get_potential_energy(), -1374.62554763 * HA_EV)
    assert fr[0].get_forces().shape == (39, 3)
    assert np.isclose(fr[0].get_forces()[0, 0], 0.000045136 * HA_EV / 0.529177210903)
    a, z = fr[0], fr[-1]
    assert a.info["charge"] == 0 and a.info["spin"] == 2 and a.info["job"] == "opt"
    assert np.isclose(a.arrays["mulliken_charges"][0], -0.926255) and a.info["q_sum_ok"] and a.info["spin_sum_ok"]
    assert np.isnan(fr[1].arrays["mulliken_charges"]).all()            # no population on intermediate steps
    assert not np.isnan(z.arrays["mulliken_charges"]).any() and z.info["converged"]  # merged final population
    assert np.isnan(a.arrays["lowdin_charges"]).all() and np.isnan(a.arrays["nbo_charges"]).all()
    write(tmp_path / "o.extxyz", fr)
    back = read(tmp_path / "o.extxyz", ":")
    assert len(back) == 4 and np.allclose(back[3].get_forces(), z.get_forces())


def test_spin_states_and_route():
    s = ["C", "O"]  # 14 electrons
    assert g16.spin_states(s, 0, bs=True) == [(1, False), (3, False), (1, True)]
    assert g16.spin_states(s, -1) == [(2, False), (4, False)]
    assert g16.spin_states(["H"], 0) == [(2, False)]  # quartet needs 3 electrons
    assert g16.route(1, False, "sp").startswith("#p RwB97XD/gen")
    r = g16.route(1, True, "sp")
    assert "UwB97XD" in r and "guess=mix" in r and "force" in r
    d = g16.deck(["Li", "F"], [[0, 0, 0], [0, 0, 1.6]], 0, 1, False, "sp", "t")
    assert d.count("****") == 2 and "\n0 1\nLi" in d


def test_registry_dedup(tmp_path):
    con = registry.connect(tmp_path / "r.sqlite")
    rec = dict(sampler="monomer", symbols=["H", "H"], positions=[[0, 0, 0], [0, 0, 0.74]], charge=0, mult=1,
               bs=False, job="sp", composition="H2", formula="H2", params={"a": 1})
    sid, new = registry.add(con, rec)
    assert new and sid.startswith("monomer-")
    rec["positions"] = [[0, 0, 0.0001], [0, 0, 0.74]]  # below 1e-3 A rounding
    assert registry.add(con, rec) == (sid, False)
    assert registry.set_status(con, sid, "parsed", "x.log")
    assert con.execute("SELECT status FROM structures").fetchone()[0] == "parsed"
