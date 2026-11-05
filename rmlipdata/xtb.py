"""GFN2-xTB driver (gas phase): single point / optimisation / Hessian / constrained scans and a nanoreactor.

Adapted from 9_build_automatic_reaction_search_agent/rxnagent/tools/nano/nano.py (ALPB removed: gas phase only).
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
from ase import Atoms

BOHR = 0.52917721
T_LOOSE, T_COMP = 3.0, 1.0       # ps per nanoreactor cycle phase
DUMP_FS = 10.0
ETEMP = (2000, 5000, 10000)      # electronic temperature ladder for SCC-failed segment retries
RHO = 1.2                        # g/cm3 liquid density for the packing sphere
HA_EV = 27.211386245988
XTB_MISSING = "xtb not found on PATH: run `conda activate rmlipdata` (or conda install -c conda-forge xtb)"


CALLS = 0  # xTB invocations since process start, shown in the generate progress line
HITS = 0   # run_xtb results served from CACHE
CACHE = None  # result cache directory (rmd generate sets <out>/xtb_cache); None = no cache


def find_xtb() -> str:
    """xtb from PATH; raises FileNotFoundError with an install hint otherwise."""
    path = shutil.which("xtb")
    if path is None:
        raise FileNotFoundError(XTB_MISSING)
    return path


# from 9_build_automatic_reaction_search_agent/rxnagent/tools/nano/nano.py
def xtb_env(threads=None) -> dict:
    """OMP threads from the argument, then RMLIPDATA_XTB_THREADS, else 2 (16 threads oversubscribes small clusters)."""
    n = str(threads or os.environ.get("RMLIPDATA_XTB_THREADS", "2"))
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = n
    env["MKL_NUM_THREADS"] = n
    env.setdefault("OMP_STACKSIZE", "1G")
    return env


def _meminfo() -> str:
    return Path("/proc/meminfo").read_text()


def wait_ram(poll_s: float = 30.0):
    """Block while used RAM exceeds RXNAGENT_RAM_STOP_MB (disabled when unset/0)."""
    stop = float(os.environ.get("RXNAGENT_RAM_STOP_MB", "0") or 0)
    if stop <= 0:
        return
    while True:
        kv = {}
        for line in _meminfo().splitlines():
            k, _, v = line.partition(":")
            if k in ("MemTotal", "MemAvailable"):
                kv[k] = float(v.split()[0]) / 1024
        if "MemTotal" not in kv or "MemAvailable" not in kv or kv["MemTotal"] - kv["MemAvailable"] <= stop:
            return
        time.sleep(poll_s)


class NanoFailed(RuntimeError):
    """xtb MD failed at every electronic temperature; traj.xyz written so far is kept."""


def r_liq(atoms: Atoms, rho: float = RHO) -> float:
    """radius (A) of a sphere holding the cluster at liquid density rho (g/cm3)"""
    v = atoms.get_masses().sum() / (0.6022 * rho)
    return (3 * v / (4 * np.pi)) ** (1 / 3)


def pack(molecules: list[Atoms], rng: np.random.Generator, rad: float) -> Atoms:
    """Place molecules at random orientations inside a sphere of radius rad, rejecting overlaps < 2 A."""
    out = Atoms()
    for m in molecules:
        m = m.copy()
        m.translate(-m.get_center_of_mass())
        for _ in range(20000):
            t = m.copy()
            t.rotate(rng.uniform(0, 360), rng.normal(size=3))
            v = rng.normal(size=3)
            t.translate(v / np.linalg.norm(v) * rad * rng.uniform() ** (1 / 3))
            if len(out) == 0 or np.min(np.linalg.norm(out.positions[:, None] - t.positions[None], axis=2)) > 2.0:
                out += t
                break
        else:
            raise RuntimeError(f"pack failed for {m.get_chemical_formula()}")
    return out


def write_xyz(path: Path, atoms: Atoms, comment: str = ""):
    lines = [str(len(atoms)), comment] + [f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in
                                          zip(atoms.get_chemical_symbols(), atoms.positions)]
    Path(path).write_text("\n".join(lines) + "\n")


def build(rundir: Path, molecules: list[tuple[list[str], np.ndarray]], label: str, chrg: int, uhf: int, T: int,
          seed: int, ncyc: int, rho: float = RHO) -> dict:
    """Pack molecules into a sphere -> rundir/init.xyz + meta.json. Deterministic per seed."""
    d = Path(rundir)
    d.mkdir(parents=True, exist_ok=True)
    mols = [Atoms(symbols=s, positions=np.asarray(x, float)) for s, x in molecules]
    rng = np.random.default_rng(1000 * seed + T + len(label))
    whole = Atoms(sum((list(m.symbols) for m in mols), []))
    at = pack(mols, rng, r_liq(whole, rho) + 1.0)
    rl = r_liq(at, rho)
    meta = dict(label=label, T=T, seed=seed, natoms=len(at), chrg=chrg, uhf=uhf, ncyc=ncyc,
                R1=round(rl + 1.5, 2), R2=round(0.7 * rl, 2))
    write_xyz(d / "init.xyz", at, label)
    (d / "meta.json").write_text(json.dumps(meta, indent=1))
    return meta


def _md_input(T: int, t_ps: float, R: float, restart: bool, dump: float = DUMP_FS) -> str:
    return f"""$md
   temp={T}
   time={t_ps}
   dump={dump}
   step=1.0
   hmass=2
   shake=0
   nvt=true
   restart={'true' if restart else 'false'}
$end
$wall
   potential=logfermi
   sphere: {R / BOHR:.4f}, all
$end
"""


def _read_scan_log(path: Path) -> tuple[list[float], list]:
    """xtbscan.log -> (energies in eV, frames as list of positions)."""
    lines, E, frames, k = path.read_text().splitlines(), [], [], 0
    while k < len(lines) and lines[k].strip():
        n = int(lines[k])
        E.append(float(re.search(r"energy:\s*(-?\d+\.\d+)", lines[k + 1]).group(1)) * HA_EV)
        frames.append([[float(v) for v in l.split()[1:4]] for l in lines[k + 2:k + 2 + n]])
        k += n + 2
    return E, frames


def _read_trj(path: Path) -> tuple[list[float], list]:
    """xtb.trj (multi-frame xyz) -> (energies in eV from the comment lines ([] when absent), frames)."""
    lines, E, frames, k = path.read_text().splitlines(), [], [], 0
    while k < len(lines):
        if not lines[k].strip():
            k += 1
            continue
        n = int(lines[k])
        m = re.search(r"energy:\s*(-?\d+\.\d+)", lines[k + 1])
        if m:
            E.append(float(m.group(1)) * HA_EV)
        frames.append([[float(v) for v in l.split()[1:4]] for l in lines[k + 2:k + 2 + n]])
        k += n + 2
    return E, frames


def _input(constraints, scan, fix) -> str:
    """$constrain / $scan / $fix blocks; xtb reads one $end per block."""
    def cblock():
        return "".join(f" distance: {i + 1}, {j + 1}, {r:.4f}\n" for i, j, r in constraints)

    if scan:
        lines = ["$constrain\n force constant=1.0\n", cblock(), "$scan\n mode=concerted\n"]
        lines += [f" {k + 1}: {a:.4f}, {b:.4f}, {int(scan['n'])}\n" for k, (a, b) in enumerate(scan["ranges"])]
        lines.append("$end\n")
        text = "".join(lines)
    elif constraints:
        text = "$constrain\n force constant=1.0\n" + cblock() + "$end\n"
    else:
        text = ""
    if fix:
        text += "$fix\n atoms: " + ",".join(str(i + 1) for i in fix) + "\n$end\n"
    return text


def _hessian(path: Path, n: int) -> np.ndarray:
    toks = [t for t in path.read_text().split() if not t.startswith("$")]
    vals = np.array([float(t) for t in toks], float)
    if vals.size != 9 * n * n:
        raise ValueError(f"hessian has {vals.size} values, expected {9 * n * n}")
    return vals.reshape(3 * n, 3 * n)


SCC_FAIL = "Self consistent charge iterator did not converge"
ETEMP_RETRY_K = 1000


def _xyz(symbols, positions):
    return f"{len(symbols)}\n\n" + "".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, (x, y, z) in zip(symbols, positions))


def run_xtb(symbols, positions, charge, uhf, workdir, mode="sp", constraints=(), scan=None, fix=(), timeout=1800,
            etemp=None, md=None):
    """_run_xtb behind a result cache: the key is a sha1 of the exact xtb input (mol.xyz text, job.inp / md.inp text,
    mode, charge, uhf, etemp, xtb binary path and mtime), so a rerun gets the same structures without running xtb.
    Only successful results are stored. md = dict(T_K, time_ps, dump_fs, wall_R_A) runs a gas-phase MD instead."""
    global HITS
    if not CACHE:
        return _run_xtb(symbols, positions, charge, uhf, workdir, mode, constraints, scan, fix, timeout, etemp, md)
    exe = shutil.which("xtb") or ""
    stamp = os.stat(exe).st_mtime_ns if exe else 0
    mdinp = _md_input(int(md["T_K"]), float(md["time_ps"]), float(md["wall_R_A"]), restart=False,
                      dump=float(md["dump_fs"])) + "$samerand\n" if md else ""
    text = (f"{mode}|{charge}|{uhf}|{etemp}|{exe}|{stamp}|{_xyz(symbols, positions)}|"
            f"{_input(list(constraints), scan, list(fix))}|{mdinp}")
    key = hashlib.sha1(text.encode()).hexdigest()
    path = Path(CACHE) / key[:2] / f"{key}.json"
    if path.exists():
        HITS += 1
        out = json.loads(path.read_text())
        if "hessian" in out:  # the only ndarray result; callers expect it back as an array
            out["hessian"] = np.asarray(out["hessian"], float)
        if out.get("etemp_K") and etemp is None:
            warnings.warn(f"xtb SCC did not converge (q={charge}, uhf={uhf}, {mode}); rerun with --etemp {ETEMP_RETRY_K}")
        return out
    out = _run_xtb(symbols, positions, charge, uhf, workdir, mode, constraints, scan, fix, timeout, etemp, md)
    if out["ok"]:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f"{key}.{os.getpid()}.tmp")
        tmp.write_text(json.dumps({k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in out.items()}))
        os.replace(tmp, path)  # atomic: parallel generate stages may write the same key
    return out


def _run_xtb(symbols, positions, charge, uhf, workdir, mode="sp", constraints=(), scan=None, fix=(), timeout=1800,
             etemp=None, md=None):
    """One GFN2 gas-phase job in a fresh temp dir under workdir, deleted afterwards. Never raises.

    mode: "sp" | "opt" | "hess"; scan = dict(n=int, ranges=[(r0, r1), ...]) implies an optimised scan over
    the constraints in order; md = dict(T_K, time_ps, dump_fs, wall_R_A) runs a $md + $wall MD instead.
    Returns {"ok", "E_eV", "positions" (opt), "frames"+"E_eV" (scan/md), "hessian" (hess),
    "charges" (sp/opt: per-atom partial charges), "err"}. If the SCC does not converge (gas-phase radical anions,
    dianions, near-degenerate states), the job is rerun once with Fermi smearing (--etemp 1000 K) and the result
    carries "etemp_K".
    """
    global CALLS
    CALLS += 1
    d = None
    try:
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)
        d = Path(tempfile.mkdtemp(dir=wd))
        (d / "mol.xyz").write_text(_xyz(symbols, positions))
        args = [find_xtb(), "mol.xyz", "--gfn", "2", "--chrg", str(charge), "--uhf", str(uhf)]
        inp = _input(list(constraints), scan, list(fix))
        if md:
            # $samerand: fixed RNG seed, so the same step gives the same trajectories on a rerun (cache, test ids)
            (d / "md.inp").write_text(_md_input(int(md["T_K"]), float(md["time_ps"]), float(md["wall_R_A"]),
                                               restart=False, dump=float(md["dump_fs"])) + "$samerand\n")
            args += ["--input", "md.inp", "--md"]
        elif inp:
            (d / "job.inp").write_text(inp)
            args += ["--input", "job.inp"]
        if mode == "opt" or scan:
            args.append("--opt")
        if mode == "hess":
            args.append("--hess")
        if etemp:
            args += ["--etemp", str(etemp)]
        wait_ram()
        try:
            r = subprocess.run(args, cwd=d, capture_output=True, text=True, timeout=timeout, env=xtb_env())
        except subprocess.TimeoutExpired:
            return {"ok": False, "err": f"timeout after {timeout} s"}
        e = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+) Eh", r.stdout)
        if md:
            bad = r.returncode or not (d / "xtbmdok").exists() or "did not converge" in r.stdout \
                or "MD is unstable" in r.stdout
        else:
            bad = r.returncode or not e
        if bad:
            if etemp is None and SCC_FAIL in r.stdout:
                warnings.warn(f"xtb SCC did not converge (q={charge}, uhf={uhf}, {mode}); rerun with --etemp {ETEMP_RETRY_K}")
                return run_xtb(symbols, positions, charge, uhf, workdir, mode, constraints, scan, fix, timeout,
                               etemp=ETEMP_RETRY_K, md=md)
            return {"ok": False, "err": (r.stdout + r.stderr)[-400:]}
        out = {"ok": True, "E_eV": float(e[-1]) * HA_EV if e else None}
        if etemp:
            out["etemp_K"] = etemp
        if mode in ("sp", "opt") and (d / "charges").exists():
            try:
                ch = [float(v) for v in (d / "charges").read_text().split()]
                if len(ch) == len(symbols):
                    out["charges"] = ch
            except ValueError:
                pass  # malformed file: charges stay unknown
        if md:
            if not (d / "xtb.trj").exists():
                return {"ok": False, "err": "no xtb.trj"}
            E, frames = _read_trj(d / "xtb.trj")
            out["frames"] = frames
            out["E_eV"] = E or [out["E_eV"]]  # no energies in the comment lines: the final total energy
        if scan:
            if not (d / "xtbscan.log").exists():
                return {"ok": False, "err": "no xtbscan.log"}
            out["E_eV"], out["frames"] = _read_scan_log(d / "xtbscan.log")
            if len(out["E_eV"]) != int(scan["n"]):
                return {"ok": False, "err": f"scan has {len(out['E_eV'])} of {scan['n']} points"}
        if mode == "opt":
            if not (d / "xtbopt.xyz").exists():
                return {"ok": False, "err": "no xtbopt.xyz"}
            lines = (d / "xtbopt.xyz").read_text().splitlines()
            out["positions"] = [[float(v) for v in l.split()[1:4]] for l in lines[2:2 + int(lines[0])]]
        if mode == "hess":
            try:
                out["hessian"] = _hessian(d / "hessian", len(symbols))
            except Exception as exc:  # missing or malformed file
                return {"ok": False, "err": str(exc)}
        return out
    except Exception as exc:  # missing binary, bad directory, ...
        return {"ok": False, "err": str(exc)}
    finally:
        if d is not None:
            shutil.rmtree(d, ignore_errors=True)


def run(rundir: Path, ps: tuple[float, float] = (T_LOOSE, T_COMP), xtb_bin: str | None = None,
        deadline: float = 0, s_per_ps: float = 0.3, segment_timeout: float = 4 * 3600,
        max_segments: int | None = None) -> Path:
    """Run the loose/compressed nanoreactor cycles; returns rundir/traj.xyz (resumable via segments_done)."""
    d = Path(rundir).resolve()
    xtb = xtb_bin or find_xtb()
    m = json.loads((d / "meta.json").read_text())
    seg_file = d / "segments_done"
    done, nbytes = 0, 0
    if seg_file.exists():
        parts_ = seg_file.read_text().split()
        done, nbytes = int(parts_[0]), int(parts_[1]) if len(parts_) > 1 else 0
    traj = d / "traj.xyz"
    if traj.exists() and nbytes and traj.stat().st_size > nbytes:
        with open(traj, "r+b") as f:
            f.truncate(nbytes)
    if done and (d / "mdrestart.committed").exists():
        shutil.copy(d / "mdrestart.committed", d / "mdrestart")
    plan = [(R, t) for _ in range(m["ncyc"]) for R, t in ((m["R1"], ps[0]), (m["R2"], ps[1]))]
    ie_file = d / "etemp_idx"
    for i, (R, t) in enumerate(plan):
        if i < done:
            continue
        if max_segments is not None and i >= max_segments:
            return traj
        if deadline and time.time() + 3600 * t * s_per_ps > deadline:
            print("deadline reached; resume later", flush=True)
            return d / "traj.xyz"
        (d / "md.inp").write_text(_md_input(m["T"], t, R, restart=bool(i)))
        ie = int(ie_file.read_text()) if ie_file.exists() else 0
        if i and (d / "mdrestart").exists():
            shutil.copy(d / "mdrestart", d / "mdrestart.prev")
        while True:
            et = ETEMP[ie]
            timed_out = False
            wait_ram()
            with open(d / "seg.out", "w") as fo:
                try:
                    rc = subprocess.call([xtb, "init.xyz", "--input", "md.inp", "--md", "--gfn", "2",
                                          "--etemp", str(et), "--chrg", str(m["chrg"]), "--uhf", str(m["uhf"])],
                                         cwd=d, stdout=fo, stderr=subprocess.STDOUT, env=xtb_env(),
                                         timeout=segment_timeout)
                except subprocess.TimeoutExpired:
                    rc, timed_out = 1, True
            out = (d / "seg.out").read_text() + ("\n[rmlipdata] segment timeout\n" if timed_out else "")
            bad = rc or not (d / "xtbmdok").exists() or "did not converge" in out or "MD is unstable" in out
            with open(d / "xtb.out", "a") as log:
                log.write(f"\n##### segment {i} R={R} A t={t} ps etemp={et} {'FAILED' if bad else 'ok'}\n" + out)
            if not bad:
                break
            (d / "xtbmdok").unlink(missing_ok=True)
            if i and (d / "mdrestart.prev").exists():
                shutil.copy(d / "mdrestart.prev", d / "mdrestart")
            ie += 1
            if ie == len(ETEMP):
                raise NanoFailed(f"xtb failed in segment {i} at all etemps {ETEMP}"
                                 + (" (segment timeout)" if timed_out else ""))
            ie_file.write_text(str(ie))
        with open(traj, "a") as f, open(d / "xtb.trj") as g:
            shutil.copyfileobj(g, f)
        (d / "xtbmdok").unlink()
        if (d / "mdrestart").exists():
            shutil.copy(d / "mdrestart", d / "mdrestart.committed")
        seg_file.write_text(f"{i + 1} {traj.stat().st_size}")
    (d / "done").touch()
    return d / "traj.xyz"
