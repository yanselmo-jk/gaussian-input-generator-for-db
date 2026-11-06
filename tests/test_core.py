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
    assert g16.spin_states(s, 0, bs=True) == [(1, True), (3, False)]  # singlet only as BS (UKS guess=mix)
    assert g16.spin_states(s, -1) == [(2, False), (4, False)]
    assert g16.spin_states(["H"], 0) == [(2, False)]  # quartet needs 3 electrons
    assert g16.route(1, False, "sp").startswith("#p RwB97XD/gen")
    assert "pop=(npa,hirshfeld)" in g16.route(1, False, "sp") and "pop=(npa,hirshfeld,always)" in g16.route(2, False, "opt")
    assert "pop=(npa,hirshfeld) " in g16.route(2, False, "opt", pop_always=False)
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


def test_parse_hangang_test_decks():
    # test_decks/ t1-t3 run on Hangang (G16 c02, job 1212786, 2026-11-04): wB97XD/def2-TZVPPD, pop=npa, iop(6/80=1)
    t1 = log_atoms(FIX.parent / "t1_li_ec_q1_m1_opt.log")      # R closed shell opt, charge +1
    assert len(t1) == 8 and t1[-1].info["converged"]
    for a in (t1[0], t1[-1]):                                  # populations on first and converged step only
        for p in ("mulliken", "lowdin", "nbo"):
            assert np.isclose(a.arrays[f"{p}_charges"].sum(), 1.0, atol=1e-3)
            assert (a.arrays[f"{p}_spins"] == 0).all()        # restricted: spins are exactly 0
    assert np.isnan(t1[3].arrays["nbo_charges"]).all()
    t2 = log_atoms(FIX.parent / "t2_li_ec_q0_m2_sp.log")[0]    # U doublet
    t3 = log_atoms(FIX.parent / "t3_ec_co_stretch_bs_sp.log")[0]  # U broken-symmetry singlet (guess=mix)
    for a, nspin in ((t2, 1.0), (t3, 0.0)):
        assert a.info["q_sum_ok"] and a.info["spin_sum_ok"]
        assert np.isclose(a.arrays["nbo_spins"].sum(), nspin, atol=1e-3)
        assert not np.isnan(a.arrays["lowdin_charges"]).any()
        assert np.isnan(a.arrays["lowdin_spins"]).all()      # iop(6/80=1) prints no Lowdin spin for open shell
    assert np.isclose(t2.arrays["nbo_spins"][0], 0.926, atol=1e-3)   # extra electron on Li (NPA)
    assert np.isclose(t3.info["s2"], 0.3253) and np.abs(t3.arrays["nbo_spins"]).max() > 0.4  # symmetry broken
    from rmlipdata.parse import reject
    from ase.calculators.singlepoint import SinglePointCalculator
    assert reject(t1[0]) == reject(t2) == reject(t3) == ""       # OMol25 A.1 filters keep the test decks
    assert t3.info["bs"] and not t2.info["bs"]
    t3.info["s2_dev"] = 1.2
    assert reject(t3) == "s2"
    t2.info["s2_dev"] = 0.3                                      # contaminated doublet (non-BS bound 0.25)
    assert reject(t2) == ""                                      # unknown sampler: reactive bound 1.1
    t2.info["sampler"] = "cluster"
    assert reject(t2) == "s2"
    t2.info["s2_dev"], t2.info["homo_eV"] = 0.0, 0.5             # unbound electron
    assert reject(t2) == "homo"
    t2.info["homo_eV"] = -5.0
    big = t2.copy()
    big.calc = SinglePointCalculator(big, energy=0.0, forces=np.full((len(big), 3), 30.0))  # |F| = 52 eV/A
    assert reject(big) == "force"


def test_hirshfeld_moments_homo():
    """pop=(npa,hirshfeld) log of EC+ doublet (Hangang G16 c02 job 1224122, 2026-11-06)."""
    from rmlipdata.parse import log_atoms, DEBYE_EA
    (a,) = log_atoms(Path(__file__).parent / "fixtures" / "h1_ec_cation_hirsh.log")
    assert np.isclose(a.arrays["hirshfeld_charges"].sum(), 1.0, atol=1e-3)
    assert np.isclose(a.arrays["hirshfeld_spins"].sum(), 1.0, atol=1e-3)
    assert np.isclose(a.arrays["cm5_charges"][1], 0.409410) and np.isclose(a.arrays["cm5_charges"].sum(), 1.0, atol=1e-3)
    assert np.allclose(a.info["dipole_eA"], np.array([-4.2644, 0.5489, 0.1585]) * DEBYE_EA)
    assert np.isclose(a.info["quadrupole_eA2"][0], -2.1565 * DEBYE_EA)
    assert np.isclose(a.info["homo_eV"], -0.61773 * 27.211386245988, atol=1e-4)  # beta HOMO is the higher one
    assert not a.info["bs"]


def test_parse_cmd_adds_registry_group(tmp_path):
    """rmd parse with a registry writes sampler / composition / group (= composition) into the extxyz info."""
    import shutil
    from ase.io import read
    from rmlipdata import cli, registry
    log = tmp_path / "h1_ec_cation_hirsh.log"
    shutil.copy(Path(__file__).parent / "fixtures" / log.name, log)
    con = registry.connect(tmp_path / "registry.sqlite")
    con.execute("INSERT INTO structures (id, sampler, composition, status) VALUES (?, ?, ?, ?)",
                (log.stem, "monomer", "EC1", "submitted"))
    con.commit()
    con.close()
    out = tmp_path / "d.extxyz"
    cli.parse_cmd([str(log)], str(out), db=str(tmp_path / "registry.sqlite"))
    (a,) = read(out, ":")
    assert a.info["group"] == "EC1" and a.info["sampler"] == "monomer"
    assert np.isfinite(a.info["homo_eV"]) and a.arrays["cm5_charges"].shape == (len(a),)


def test_flag_paths_energy_and_spin(tmp_path):
    """Smooth harmonic path passes; an energy offset on one frame or a moved NPA spin marks its segments."""
    from ase import Atoms
    from ase.io import read, write
    from ase.calculators.singlepoint import SinglePointCalculator
    from rmlipdata.parse import flag_paths, path_key
    k, frames = 5.0, []
    for i, x in enumerate(np.linspace(1.0, 2.0, 6)):
        a = Atoms("HF", positions=[[0, 0, 0], [x, 0, 0]])
        e, f = 0.5 * k * (x - 1.2) ** 2, k * (x - 1.2)
        a.calc = SinglePointCalculator(a, energy=e, forces=[[f, 0, 0], [-f, 0, 0]])
        a.arrays["nbo_spins"] = np.array([0.9, 0.1])
        a.info.update(sampler="bond_scan1", composition="HF1", charge=0, spin=2, bs=False,
                      path="bond:0-1", path_pos=float(i))
        frames.append(a)
    assert flag_paths(frames) == 0 and all(a.info["path_flag"] == "ok" for a in frames)
    frames[3].calc.results["energy"] += 3.0                 # SCF landed on another solution
    frames[5].arrays["nbo_spins"] = np.array([0.1, 0.9])   # unpaired electron jumped
    assert flag_paths(frames) == 3
    assert [a.info["path_flag"] for a in frames] == ["ok", "ok", "energy", "energy", "energy,spin", "spin"]
    frames[3].info["rattled"] = True                       # rattled frames are checked like the others
    assert flag_paths(frames) == 3
    assert [a.info["path_flag"] for a in frames] == ["ok", "ok", "energy", "energy", "energy,spin", "spin"]
    write(tmp_path / "p.extxyz", frames, format="extxyz")
    assert read(tmp_path / "p.extxyz", ":")[0].info["path_flag"] == "ok"
    assert path_key("hop", {"kind": "hop", "t": 0.5, "cation_index": 0, "from_donor": 1, "to_donor": 2})[1] == 0.5
    assert path_key("bond_scan1", {"bond": [0, 1], "index": 3, "vertical_dq": 1}) is None
    assert path_key("bond_scan1", {"bond": [0, 1], "index": 3, "cfg_tag": "ab12cd34"}) == ("ab12cd34/bond:0-1", 3.0)


def test_parse_dir_layout_id_from_title(tmp_path):
    """<dir>/in.log from another machine: id comes from the deck title; unfinished logs leave the status alone."""
    from rmlipdata import cli, registry
    src = (FIX.parent / "t2_li_ec_q0_m2_sp.log").read_text()
    sid = "monomer-0123456789ab"
    done = src.replace(" t2 Li EC reduced doublet sp\n", f" {sid}\n", 1)
    cut = done[:done.index("Normal termination")]
    con = registry.connect(tmp_path / "registry.sqlite")
    for i, st in ((sid, "submitted"), ("monomer-ffffffffffff", "submitted"), ("monomer-eeeeeeeeeeee", "submitted")):
        con.execute("INSERT INTO structures (id, sampler, composition, status) VALUES (?, ?, ?, ?)", (i, "monomer", "Li1 EC1", st))
    con.commit()
    for d, text in ((sid, done), ("monomer-ffffffffffff", cut.replace(sid, "monomer-ffffffffffff")),
                    ("monomer-eeeeeeeeeeee", cut.replace(sid, "monomer-eeeeeeeeeeee") + " Error termination via Lnk1e\n")):
        (tmp_path / "back" / d).mkdir(parents=True)
        (tmp_path / "back" / d / "in.log").write_text(text)
    cli.parse_cmd([str(tmp_path / "back")], str(tmp_path / "d.extxyz"), db=str(tmp_path / "registry.sqlite"))
    st = dict(con.execute("SELECT id, status FROM structures"))
    assert st == {sid: "parsed", "monomer-ffffffffffff": "submitted", "monomer-eeeeeeeeeeee": "failed"}
    assert read(tmp_path / "d.extxyz", ":")[0].info["id"] == sid
    # duplicate copies (a failed one sorting last), CRLF, a run that died before the title, a slurm file
    (tmp_path / "back" / "z_dup").mkdir()
    (tmp_path / "back" / "z_dup" / "in.log").write_text(cut + " Error termination via Lnk1e\n")
    (tmp_path / "back" / sid / "in.log").write_text(done.replace("\n", "\r\n"))
    (tmp_path / "back" / "monomer-ffffffffffff" / "in.log").write_text(" Entering Gaussian System\n Error termination\n")
    (tmp_path / "back" / "slurm-1.out").write_text("job started\n")
    cli.parse_cmd([str(tmp_path / "back")], str(tmp_path / "d.extxyz"), db=str(tmp_path / "registry.sqlite"))
    st = dict(con.execute("SELECT id, status FROM structures"))
    assert st == {sid: "parsed", "monomer-ffffffffffff": "failed", "monomer-eeeeeeeeeeee": "failed"}
    assert len(read(tmp_path / "d.extxyz", ":")) == 2                 # sid once (normal copy) + eeee once
    import pytest
    with pytest.raises(SystemExit):
        cli.main(["status", "--registry", str(tmp_path / "typo.sqlite")])
    assert not (tmp_path / "typo.sqlite").exists()
    (tmp_path / "back" / sid / "in.gjf").write_text(f"%mem=1GB\n#p sp\n\n{sid}\n\n0 2\n")
    assert cli._ids([str(tmp_path / "back" / sid / "in.gjf")]) == [sid]
    cli.submit_cmd(str(tmp_path / "registry.sqlite"), ["monomer-ffffffffffff"])  # rerun of a failed deck
    assert con.execute("SELECT status FROM structures WHERE id='monomer-ffffffffffff'").fetchone()[0] == "submitted"
    import shutil
    shutil.rmtree(tmp_path / "back" / sid)
    shutil.rmtree(tmp_path / "back" / "monomer-eeeeeeeeeeee")
    shutil.rmtree(tmp_path / "back" / "z_dup")
    cli.parse_cmd([str(tmp_path / "back")], str(tmp_path / "d.extxyz"), db=str(tmp_path / "registry.sqlite"))
    assert (tmp_path / "d.extxyz").read_text() == ""                  # no stale frames from the previous parse
