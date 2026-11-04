"""Gaussian16 decks: wB97X-D/def2-TZVPPD (gen basis), gas phase, forces + Mulliken / NPA / Loewdin populations."""
import re
from pathlib import Path

from ase.data import atomic_numbers

BASIS_FILE = Path(__file__).parent / "data" / "basis" / "def2-TZVPPD.gbs"
# iop(6/80=1): Loewdin population; pop=npa: NBO natural population analysis.
# UNVERIFIED until the test decks run: Loewdin/NPA spin output format (see README).
ROUTE = "#p {ref}wB97XD/gen nosymm int=ultrafine scf=(xqc,maxcycle=300) pop=(npa) iop(6/80=1)"
JOB = {"sp": "force", "opt": "opt=(maxcycles=200)"}


def _basis():
    text = BASIS_FILE.read_text()
    return {e: f"{e}     0\n{b}\n****\n" for e, b in
            re.findall(r"^([A-Z][a-z]?)\s+0\s*\n(.*?)\n\*{4}", text, re.M | re.S)}


BASIS = _basis()


def n_electrons(symbols, charge):
    return sum(atomic_numbers[s] for s in symbols) - charge


def spin_states(symbols, charge, offsets=(0, 2), bs=False):
    """[(mult, bs)]: lowest multiplicity + offsets; even-electron systems add a broken-symmetry singlet if bs."""
    n = n_electrons(symbols, charge)
    if n < 1:
        return []
    low = 1 if n % 2 == 0 else 2
    out = [(low + o, False) for o in offsets if low + o - 1 <= n]
    if bs and low == 1 and n >= 2:
        out.append((1, True))
    return out


def route(mult, bs, job):
    """R for closed-shell singlet, U otherwise; BS singlet adds guess=mix."""
    ref = "U" if (mult > 1 or bs) else "R"
    r = ROUTE.format(ref=ref) + " " + JOB[job]
    return r + " guess=mix" if bs else r


def deck(symbols, positions, charge, mult, bs, job, title, nproc=16, mem="32GB"):
    missing = sorted(set(symbols) - set(BASIS))
    if missing:
        raise ValueError(f"no def2-TZVPPD block for {missing} in {BASIS_FILE}")
    ne = n_electrons(symbols, charge)
    if not 1 <= mult <= ne + 1 or (ne + mult - 1) % 2:
        raise ValueError(f"charge {charge} / multiplicity {mult} impossible for {ne} electrons")
    atoms = "".join(f"{s:<2} {x:14.8f} {y:14.8f} {z:14.8f}\n" for s, (x, y, z) in zip(symbols, positions))
    gen = "".join(BASIS[e] for e in dict.fromkeys(symbols))
    return (f"%nprocshared={nproc}\n%mem={mem}\n{route(mult, bs, job)}\n\n{title}\n\n"
            f"{charge} {mult}\n{atoms}\n{gen}\n")
