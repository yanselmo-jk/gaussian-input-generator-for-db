"""Gaussian16 log -> frames (positions, energy, forces, Mulliken / Loewdin / NPA / Hirshfeld charges and spins, CM5
charges, dipole, traceless quadrupole, HOMO energy) -> extxyz.

A log is split into job steps (Link1), each step into frames at every 'Input orientation:' (route uses nosymm, so
forces are in the same frame as the coordinates). A frame needs 'SCF Done' and a force table. A trailing segment
without energy/forces at the same geometry (the population re-run after 'Stationary point found') is merged into
the previous frame. Populations that a frame does not print are NaN.
"""
import re
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import chemical_symbols

HA_EV = 27.211386245988
BOHR_A = 0.529177210903
NUM = r"-?\d+\.\d+"
POPS = ("mulliken", "lowdin", "nbo", "hirshfeld")
DEBYE_EA = 0.20819433622  # 1 Debye in e*A

STEP_SPLIT = re.compile(r"\n(?= Initial command:| Link1:  Proceeding to internal job step)")
GEOM = re.compile(r"Input orientation:.*?-{5,}\n.*?-{5,}\n(.*?)\n\s*-{5,}", re.S)
SCF = re.compile(rf"SCF Done:\s+E\(\S+\)\s*=\s*({NUM})")
# [^\n] (not .) so the match cannot run past the header under re.S (sei mecp.py pitfall)
FORCES = re.compile(r"Forces \(Hartrees/Bohr\)\n[^\n]*\n[^\n]*-{5,}\n(.*?)\n\s*-{5,}", re.S)
S2 = re.compile(rf"S\*\*2 before annihilation\s+({NUM})")
CHMULT = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
ROUTE = re.compile(r"\n -{10,}\n( #.*?)\n -{10,}\n", re.S)
ROW = re.compile(rf"^\s*\d+\s+[A-Z][a-z]?\s+({NUM})(?:\s+({NUM}))?\s*$")
MULLIKEN = re.compile(r"^ Mulliken charges( and spin densities)?:\s*$", re.M)
# header confirmed on G16 c02 (Hangang 2026-11-04): open shell prints charges only, no Lowdin spin
LOWDIN = re.compile(r"^ Lowdin (?:Atomic )?Charges( and Spin Densities)?\s*:\s*$", re.M | re.I)
NPA = re.compile(r"Summary of Natural Population Analysis:")
NPA_ROW = re.compile(rf"^\s*[A-Z][a-z]?\s+\d+\s+({NUM})\s+{NUM}\s+{NUM}\s+{NUM}\s+({NUM})\s*$")
SPIN_SECTION = re.compile(r"(Alpha|Beta)\s+spin\s+orbitals", re.I)
# pop=hirshfeld block, format confirmed on G16 c02 (Hangang job 1224122, 2026-11-06); S-H is 0 for closed shells
HIRSH = re.compile(r"Hirshfeld charges, spin densities, dipoles, and CM5 charges using[^\n]*\n[^\n]*Q-CM5[^\n]*\n")
HIRSH_ROW = re.compile(rf"^\s*\d+\s+[A-Z][a-z]?\s+({NUM})\s+({NUM})\s+{NUM}\s+{NUM}\s+{NUM}\s+({NUM})\s*$")
DIPOLE = re.compile(rf"Dipole moment \(field-independent basis, Debye\):\n\s+X=\s*({NUM})\s+Y=\s*({NUM})\s+Z=\s*({NUM})")
QUAD = re.compile(rf"Traceless Quadrupole moment \(field-independent basis, Debye-Ang\):\n"
                  rf"\s+XX=\s*({NUM})\s+YY=\s*({NUM})\s+ZZ=\s*({NUM})\n\s+XY=\s*({NUM})\s+XZ=\s*({NUM})\s+YZ=\s*({NUM})")
OCC = re.compile(r"^ +(?:Alpha|Beta) +occ\. eigenvalues --(.*)$", re.M)


def _table(text, start, n):
    """n rows of 'idx sym q [spin]' after position start -> (q, spin|None)."""
    q, s = [], []
    for line in text[start:].splitlines()[1:]:
        m = ROW.match(line)
        if not m:
            if q:
                break
            continue
        q.append(float(m[1]))
        s.append(None if m[2] is None else float(m[2]))
        if len(q) == n:
            break
    if len(q) != n:
        return None
    return np.array(q), (None if s[0] is None else np.array(s))


def _pop(text, header, n, restricted):
    """Last block of a Mulliken-style table; spins are 0 only for a restricted closed-shell calculation."""
    m = list(header.finditer(text))
    if not m:
        return None
    t = _table(text, m[-1].end(), n)
    if t is None:
        return None
    q, s = t
    if s is None and restricted:
        s = np.zeros(n)
    return q, s


def _npa(text, n):
    """Total-density NPA charges; spin = alpha - beta natural population when both sections are printed."""
    marks = [(m.start(), m[1].lower()) for m in SPIN_SECTION.finditer(text)]
    blocks = {}
    for m in NPA.finditer(text):
        kind = next((k for p, k in reversed(marks) if p < m.start()), "total")
        rows = []
        for line in text[m.end():].splitlines():
            r = NPA_ROW.match(line)
            if r:
                rows.append((float(r[1]), float(r[2])))
                if len(rows) == n:
                    break
            elif rows:
                break
        if len(rows) == n:
            blocks[kind] = np.array(rows)
    if "total" not in blocks:
        return None
    spin = blocks["alpha"][:, 1] - blocks["beta"][:, 1] if {"alpha", "beta"} <= set(blocks) else None
    return blocks["total"][:, 0], spin


def _hirsh(text, n):
    """Last Hirshfeld block: (charges, spins), CM5 charges."""
    m = list(HIRSH.finditer(text))
    if not m:
        return None, None
    rows = []
    for line in text[m[-1].end():].splitlines():
        r = HIRSH_ROW.match(line)
        if not r:
            break
        rows.append([float(x) for x in r.groups()])
    if len(rows) != n:
        return None, None
    a = np.array(rows)
    return (a[:, 0], a[:, 1]), a[:, 2]


def _moments(text):
    """Last dipole (e*A), traceless quadrupole (e*A^2: XX YY ZZ XY XZ YZ) and HOMO energy (eV, max of alpha/beta)."""
    d, q = DIPOLE.findall(text), QUAD.findall(text)
    occ = OCC.findall(text)
    homo = None
    if occ:  # the last population print holds alpha then beta lines; take the highest occupied level of both
        last = text.rfind(" Alpha  occ. eigenvalues")
        vals = [float(v) for line in OCC.findall(text[last:]) for v in re.findall(NUM, line)]
        homo = max(vals) * HA_EV if vals else None
    return (np.array([float(x) for x in d[-1]]) * DEBYE_EA if d else None,
            np.array([float(x) for x in q[-1]]) * DEBYE_EA if q else None, homo)


def _geom(block):
    rows = [l.split() for l in block.splitlines()]
    return [chemical_symbols[int(r[1])] for r in rows], np.array([[float(x) for x in r[3:6]] for r in rows])


def parse_step(text):
    """One job step -> dict(route, charge, mult, normal, stationary, frames=[...])."""
    r = ROUTE.search(text)
    route = "".join(l[1:] for l in r[1].split("\n")).strip() if r else ""
    cm = CHMULT.search(text)
    charge, mult = (int(cm[1]), int(cm[2])) if cm else (None, None)
    restricted = mult == 1 and not re.search(r"(^|[\s#])(p\s+)?U\w*/", route, re.I)
    frames = []
    cuts = [m.start() for m in re.finditer(r"Input orientation:", text)] + [len(text)]
    for a, b in zip(cuts, cuts[1:]):
        seg = text[a:b]
        g = GEOM.search(seg)
        if not g:
            continue
        symbols, pos = _geom(g[1])
        n = len(symbols)
        fr = {"symbols": symbols, "positions": pos}
        for key, pat in (("mulliken", MULLIKEN), ("lowdin", LOWDIN)):
            fr[key] = _pop(seg, pat, n, restricted)
        fr["nbo"] = _npa(seg, n)
        fr["hirshfeld"], fr["cm5"] = _hirsh(seg, n)
        fr["dipole"], fr["quadrupole"], fr["homo"] = _moments(seg)
        if restricted and fr["nbo"] and fr["nbo"][1] is None:  # same rule as Mulliken/Loewdin
            fr["nbo"] = (fr["nbo"][0], np.zeros(n))
        e, f, s2 = SCF.findall(seg), FORCES.search(seg), S2.findall(seg)
        if not e or not f:
            prev = frames[-1] if frames else None
            if prev is not None and np.allclose(prev["positions"], pos, atol=1e-4):
                for k in POPS:
                    prev[k] = fr[k] or prev[k]
                for k in ("cm5", "dipole", "quadrupole", "homo"):
                    prev[k] = fr[k] if fr[k] is not None else prev[k]
            continue
        fr["energy"] = float(e[-1]) * HA_EV
        fr["forces"] = np.array([[float(x) for x in l.split()[2:5]] for l in f[1].splitlines()]) * HA_EV / BOHR_A
        fr["s2"] = float(s2[-1]) if s2 else None
        frames.append(fr)
    return {"route": route, "charge": charge, "mult": mult, "normal": "Normal termination of Gaussian" in text,
            "stationary": "Stationary point found" in text, "frames": frames}


def parse_log(path):
    text = Path(path).read_text(errors="replace")
    return [parse_step(c) for c in STEP_SPLIT.split(text) if ROUTE.search(c)]


def to_atoms(step, sid, step_idx):
    """Frames of one step -> ase.Atoms with info/arrays in the 6_build_charge-aware-mlip convention."""
    out = []
    is_opt = bool(re.search(r"\bopt\b", step["route"], re.I) or re.search(r"\bopt=", step["route"], re.I))
    for k, fr in enumerate(step["frames"]):
        n = len(fr["symbols"])
        at = Atoms(fr["symbols"], positions=fr["positions"])
        nan = np.full(n, np.nan)
        for p in POPS:
            q, s = fr[p] if fr[p] else (None, None)
            at.arrays[f"{p}_charges"] = nan.copy() if q is None else q
            at.arrays[f"{p}_spins"] = nan.copy() if s is None else s
        at.arrays["cm5_charges"] = nan.copy() if fr["cm5"] is None else fr["cm5"]
        S = (step["mult"] - 1) / 2
        last = k == len(step["frames"]) - 1
        at.info.update(id=sid, step=step_idx, frame=k, charge=step["charge"], spin=step["mult"],
                       job="opt" if is_opt else "sp", normal_termination=step["normal"],
                       converged=bool(last and (step["stationary"] if is_opt else step["normal"])),
                       s2=np.nan if fr["s2"] is None else fr["s2"],
                       s2_dev=np.nan if fr["s2"] is None else fr["s2"] - S * (S + 1),
                       bs=bool(re.search(r"guess=mix", step["route"], re.I)),
                       homo_eV=np.nan if fr["homo"] is None else fr["homo"],
                       dipole=np.full(3, np.nan) if fr["dipole"] is None else fr["dipole"],
                       quadrupole=np.full(6, np.nan) if fr["quadrupole"] is None else fr["quadrupole"])
        mq = fr["mulliken"]
        if mq:
            at.info["q_sum_ok"] = bool(abs(mq[0].sum() - step["charge"]) < 1e-3)
            if mq[1] is not None:
                at.info["spin_sum_ok"] = bool(abs(mq[1].sum() - (step["mult"] - 1)) < 0.05)
        at.calc = SinglePointCalculator(at, energy=fr["energy"], forces=fr["forces"])
        out.append(at)
    return out


# OMol25 Appendix A.1 quality filters (arXiv 2505.08762). Their "S2 < 0.5 / < 1.1" is read as the deviation
# <S2> - S(S+1): an absolute bound would drop every triplet. "Metal-containing" is read as transition metal;
# the def2-TZVPPD file here covers H-Kr, so that means Sc-Zn. The 1.1 bound is kept for broken-symmetry
# singlets (full unpairing gives a deviation of ~1); other open-shell decks use a tighter 0.25 (our choice:
# e.g. a doublet with <S2> 1.8 is a contaminated state, not a label). A positive HOMO means the extra electron
# is unbound (it sits in the diffuse basis tail).
MAX_FORCE = 50.0      # eV/A, largest per-atom force
MAX_S2_DEV = 1.1      # broken-symmetry singlets
MAX_S2_DEV_PURE = 0.25
MAX_S2_DEV_TM = 0.5


def reject(at):
    """'force' or 's2' if the frame fails the OMol25 force / spin-contamination filters, else ''."""
    if np.linalg.norm(at.get_forces(), axis=1).max() > MAX_FORCE:
        return "force"
    d = at.info.get("s2_dev", np.nan)
    bound = MAX_S2_DEV if at.info.get("bs") else MAX_S2_DEV_PURE
    if any(21 <= z <= 30 for z in at.numbers):
        bound = min(bound, MAX_S2_DEV_TM)
    if np.isfinite(d) and d > bound:
        return "s2"
    if at.info.get("homo_eV", np.nan) > 0:
        return "homo"
    return ""


# Path continuity check. Frames of one generated path (same sampler, composition, charge, multiplicity, BS flag and
# path key) are ordered by their path position; for neighbours a, b the energy change is compared with the work of
# the forces along the straight segment, dE_pred = -1/2 (F_a + F_b) . (R_b - R_a) (trapezoid rule, exact for a
# quadratic PES; rigid translation / rotation do no work). A mismatch means the SCF landed on a different solution
# (or the PES is not smooth on that segment); a jump of an NPA atomic spin means the unpaired electron moved.
# Both frames of a bad segment are marked, nothing is dropped. Thresholds are our choice, not literature values.
# The trapezoid estimate is only good for short segments: on xTB paths of the c4 test config (2026-11-06, our
# calculation) exact neighbours differ by 0.42 eV median (max 9.2 eV) because atoms move 0.6 A median between
# frames, and rattled frames by 2.6 eV median. So the energy check uses only unrattled neighbours whose largest
# atomic displacement is below 0.5 A (xTB max residual there 1.66 eV) and flags only gross jumps above 2 eV.
PATH_JUMP_EV = 2.0
PATH_MAX_STEP_A = 0.5
PATH_SPIN_E = 0.5


def path_key(sampler, params):
    """(key, position) of a generated path frame from its registry params, or None (not on a path)."""
    p = params or {}
    if "vertical_dq" in p or "peak_rattle" in p:
        return None  # vertical charge decks and rattled peak copies branch off the path
    if sampler in ("bond_scan1", "bond_scan2") and "bond" in p:
        return f"bond:{p['bond'][0]}-{p['bond'][1]}", float(p.get("index", 0))
    if sampler in ("reaction", "reaction2") and "template" in p:
        return f"{p['template']}:{p.get('pairs')}", float(p.get("index", 0))
    if sampler == "hop" and p.get("kind") == "exchange_out":
        return f"exchange:{p.get('cation_index')}:{p.get('donor')}", float(p.get("index", 0))
    if sampler == "hop" and "t" in p:
        return f"hop:{p.get('cation_index')}:{p.get('from_donor')}>{p.get('to_donor')}", float(p["t"])
    return None


def flag_paths(frames):
    """Set info path_flag ('ok', 'energy', 'spin' or 'energy,spin'; extxyz cannot store '') on path frames."""
    groups = {}
    for a in frames:
        if "path" in a.info:
            a.info["path_flag"] = "ok"
            key = (a.info.get("sampler"), a.info.get("composition"), a.info["charge"], a.info["spin"],
                   bool(a.info.get("bs")), a.info["path"])
            groups.setdefault(key, []).append(a)
    n = 0
    for g in groups.values():
        g.sort(key=lambda a: a.info["path_pos"])
        for a, b in zip(g, g[1:]):
            if a.get_chemical_symbols() != b.get_chemical_symbols():
                continue
            bad = set()
            de = b.get_potential_energy() - a.get_potential_energy()
            pred = -0.5 * float(((a.get_forces() + b.get_forces()) * (b.positions - a.positions)).sum())
            short = np.linalg.norm(b.positions - a.positions, axis=1).max() < PATH_MAX_STEP_A
            if short and not (a.info.get("rattled") or b.info.get("rattled")) and abs(de - pred) > PATH_JUMP_EV:
                bad.add("energy")
            sa, sb = a.arrays.get("nbo_spins"), b.arrays.get("nbo_spins")
            if (a.info["spin"] > 1 or a.info.get("bs")) and sa is not None and sb is not None \
                    and np.isfinite(sa).all() and np.isfinite(sb).all() and np.abs(sb - sa).max() > PATH_SPIN_E:
                bad.add("spin")
            for x in (a, b) if bad else ():
                x.info["path_flag"] = ",".join(sorted(bad | (set(x.info["path_flag"].split(",")) - {"ok"})))
            n += bool(bad)
    return n


def log_atoms(path):
    """All frames of a log; id = file stem (= registry id)."""
    sid = Path(path).stem
    return [a for i, st in enumerate(parse_log(path)) for a in to_atoms(st, sid, i)]
