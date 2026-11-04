"""Gaussian16 log -> frames (positions, energy, forces, Mulliken / Loewdin / NPA charges and spins) -> extxyz.

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
POPS = ("mulliken", "lowdin", "nbo")

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
# UNVERIFIED header text for iop(6/80=1); confirm on the test decks
LOWDIN = re.compile(r"^ Lowdin (?:Atomic )?Charges( and Spin Densities)?\s*:\s*$", re.M | re.I)
NPA = re.compile(r"Summary of Natural Population Analysis:")
NPA_ROW = re.compile(rf"^\s*[A-Z][a-z]?\s+\d+\s+({NUM})\s+{NUM}\s+{NUM}\s+{NUM}\s+({NUM})\s*$")
SPIN_SECTION = re.compile(r"(Alpha|Beta)\s+spin\s+orbitals", re.I)


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
        e, f, s2 = SCF.findall(seg), FORCES.search(seg), S2.findall(seg)
        if not e or not f:
            prev = frames[-1] if frames else None
            if prev is not None and np.allclose(prev["positions"], pos, atol=1e-4):
                for k in POPS:
                    prev[k] = fr[k] or prev[k]
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
        S = (step["mult"] - 1) / 2
        last = k == len(step["frames"]) - 1
        at.info.update(id=sid, step=step_idx, frame=k, charge=step["charge"], spin=step["mult"],
                       job="opt" if is_opt else "sp", normal_termination=step["normal"],
                       converged=bool(last and (step["stationary"] if is_opt else step["normal"])),
                       s2=np.nan if fr["s2"] is None else fr["s2"],
                       s2_dev=np.nan if fr["s2"] is None else fr["s2"] - S * (S + 1))
        mq = fr["mulliken"]
        if mq:
            at.info["q_sum_ok"] = bool(abs(mq[0].sum() - step["charge"]) < 1e-3)
            if mq[1] is not None:
                at.info["spin_sum_ok"] = bool(abs(mq[1].sum() - (step["mult"] - 1)) < 0.05)
        at.calc = SinglePointCalculator(at, energy=fr["energy"], forces=fr["forces"])
        out.append(at)
    return out


def log_atoms(path):
    """All frames of a log; id = file stem (= registry id)."""
    sid = Path(path).stem
    return [a for i, st in enumerate(parse_log(path)) for a in to_atoms(st, sid, i)]
