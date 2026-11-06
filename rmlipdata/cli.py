"""rmd command line: generate Gaussian16 decks from xTB-sampled geometries, parse logs, inspect the registry."""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import time
import warnings
import zlib
from pathlib import Path

import numpy as np
import yaml
from ase import Atoms
from ase.io import read as ase_read
from ase.io import write as ase_write
from rdkit import Chem

from . import build, g16, parse, registry, sample, species, xtb


def _ctx(cfg, cat, out):
    from types import SimpleNamespace
    return SimpleNamespace(cat=cat, cfg=cfg, workdir=out, rng=np.random.default_rng(int(cfg.get("seed", 0))),
                           clusters={})


def _max_atoms(cfg):
    return int(cfg.get("limits", {}).get("max_atoms", 10 ** 9))


def _enabled(sc):
    """YAML 1.1 turns the unquoted key `on` into True, so accept both spellings."""
    return bool(sc.get("on", sc.get(True, True)))


def _combos(cfg, cat):
    return build.combos(cat, cfg["species"], cfg.get("cluster", {}), _max_atoms(cfg),
                        np.random.default_rng(int(cfg.get("seed", 0)) + zlib.crc32(b"cluster")))


def _plan(cfg, cat, names):
    combos = _combos(cfg, cat)
    sp = [n for k in ("cations", "anions", "solvents") for n in cfg["species"].get(k, [])]
    npairs = len(sp) * (len(sp) + 1) // 2  # all pairs, cation pairs included
    plan = {}
    for name in names:
        scfg = cfg["samplers"][name]
        if name == "monomer":
            plan[name] = len(sp)
        elif name == "atoms":
            dq0 = 0 in [int(q) for q in cfg.get("charge", {}).get("delta_q", [0])]
            elems = {a.GetSymbol() for k in ("cations", "anions", "solvents")
                     for n in cfg["species"].get(k, [])
                     for m in [Chem.MolFromSmiles(cat[n].smiles or "")] if m is not None
                     for a in Chem.AddHs(m).GetAtoms()}
            plan[name] = len(elems) if dq0 else 0  # the atoms sampler only runs the delta_q 0 stage
        elif name == "pair_scan":
            plan[name] = min(npairs, int(scfg.get("max_pairs", npairs)))
        elif name == "cluster_md":
            plan[name] = min(len(sample.md_compositions(cfg, cat, scfg)), int(scfg.get("max_clusters", 20)))
        elif name == "cluster":
            plan[name] = len([c for c in combos if len(c) >= 2])
        else:
            plan[name] = len(combos)
    return plan


def _dry(cfg, cat, names):
    for name, k in _plan(cfg, cat, names).items():
        print(f"{name}: {k}")


def _stage(cfg, cat, out, dq, name, put):
    """One (delta_q, sampler) stage with its own ctx. The rng depends only on seed and sampler name and clusters
    are seeded per (combo, q), so a stage gives the same records alone, in sequence or in parallel."""
    ctx = _ctx(cfg, cat, out)
    ctx.combos = _combos(cfg, cat)
    ctx.rng = np.random.default_rng(int(cfg.get("seed", 0)) + zlib.crc32(name.encode()))
    ctx.progress = lambda label: put(("unit", label, xtb.CALLS, xtb.HITS))
    cap = int(cfg["samplers"][name].get("max_atoms", _max_atoms(cfg)))
    for rec in sample.SAMPLERS[name](ctx, cfg["samplers"][name], dq):
        if len(rec["symbols"]) <= cap:
            put(("rec", rec))


def _child(i, args, q):
    """Process target: stream stage i through q; always ends with ("end", warnings, calls, hits, traceback)."""
    import traceback
    put = lambda m: q.put((i,) + m)  # noqa: E731
    xtb.CALLS = xtb.HITS = 0
    err = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            _stage(*args, put)
        except BaseException:
            err = traceback.format_exc()
    put(("end", [str(w.message) for w in caught], xtb.CALLS, xtb.HITS, err))


def generate(cfg, only=None, dry=False):
    import multiprocessing
    import queue
    from collections import Counter
    cat = species.catalog()
    names = [n for n, sc in cfg["samplers"].items()
             if _enabled(sc) and (not only or n in only)]
    unknown = [n for n in names if n not in sample.SAMPLERS]
    if unknown:
        raise SystemExit(f"unknown sampler(s) {unknown}; known: {list(sample.SAMPLERS)} "
                         "(bond_scan became bond_scan1 = xtb relaxed scan and bond_scan2 = IDPP interpolation)")
    if dry:
        _dry(cfg, cat, names)
        return
    if shutil.which("xtb") is None:  # without xtb every structure would silently stay at its start geometry
        raise SystemExit(xtb.XTB_MISSING)
    out = Path(cfg["out"])
    out.mkdir(parents=True, exist_ok=True)
    xcfg = cfg.get("xtb", {})
    threads = int(xcfg.get("threads", 2))
    nmax = max(1, int(xcfg.get("n_max_process", 2)))
    os.environ["RMLIPDATA_XTB_THREADS"] = str(threads)
    xtb.CACHE = out / "xtb_cache"
    xtb.TIMEOUT = int(cfg.get("xtb", {}).get("timeout_s", 1800))  # also species / cluster opts
    con = registry.connect(out / "registry.sqlite")
    gcfg = cfg.get("g16", {})
    nproc, mem = int(gcfg.get("nproc", 16)), gcfg.get("mem", "32GB")
    pop_always = bool(gcfg.get("pop_always", True))
    offsets = tuple(cfg.get("charge", {}).get("mult_offsets", (0, 2)))
    redox = float(cfg.get("charge", {}).get("redox_fraction", 1.0))
    high_spin = float(cfg.get("charge", {}).get("high_spin_fraction", 1.0))
    vertical = float(cfg.get("charge", {}).get("vertical_fraction", 0.2))
    dedup = float(cfg.get("limits", {}).get("dedup_A", 0.03))
    qlo, qhi = cfg.get("charge", {}).get("range", (-2, 2))
    added = {n: 0 for n in names}
    hsfrac = {n: float(cfg["samplers"][n].get("high_spin_fraction", high_spin)) for n in names}
    # settings that shape the geometries: paths from another seed or n_points in the same out dir stay apart in
    # rmd parse (path_key); high_spin_fraction only picks multiplicities
    def canon(x):  # str keys: YAML reads 'on' as True, and int keys (max_components) cannot be sorted with str
        return {str(k): canon(v) for k, v in x.items()} if isinstance(x, dict) else \
            [canon(v) for v in x] if isinstance(x, (list, tuple)) else x

    tags = {n: hashlib.sha1(json.dumps(canon([cfg.get("seed", 0), cfg.get("species"), cfg.get("cluster"),
                                              {k: v for k, v in cfg["samplers"][n].items() if k != "high_spin_fraction"}]),
                                       sort_keys=True, default=str).encode()).hexdigest()[:8] for n in names}
    dedup_dropped, seen = Counter(), {}

    dq_list = cfg.get("charge", {}).get("delta_q", [0])
    plan = _plan(cfg, cat, names)
    todo = [(i, dq, name) for i, (dq, name) in enumerate(((dq, n) for dq in dq_list for n in names), 1)]
    n_stages = len(todo)
    extra = "".join(f", {k} {v:g}" for k, v in (("redox_fraction", redox), ("high_spin_fraction", high_spin),
                                                ("vertical_fraction", vertical)) if v < 1.0)
    if dedup:
        extra += f", dedup_A {dedup:g}"
    for name in names:
        print(f"  {name}: {plan[name]} units x {len(dq_list)} delta_q{extra}", file=sys.stderr)
    print(f"{n_stages} stages (delta_q x sampler), {nmax} at a time, xtb threads {threads} each", file=sys.stderr)

    tty = sys.stderr.isatty()
    start = time.time()
    stage, running, msgs, errors = {}, {}, Counter(), []
    st = {"done": 0, "new": 0, "last": 0.0}

    def status():
        s = int(time.time() - start)
        h, s = divmod(s, 3600)
        mm, ss = divmod(s, 60)
        clock = f"{h}:{mm:02d}:{ss:02d}" if h else f"{mm:02d}:{ss:02d}"
        act = ", ".join(f"dq={float(s['dq']):+g} {s['name']} {s['u']}/{s['k']} {s['label']}".rstrip()
                        for i, s in stage.items() if i in running)
        calls, hits = sum(s["calls"] for s in stage.values()), sum(s["hits"] for s in stage.values())
        return f"[{st['done']}/{n_stages}] {act} | new {st['new']} | xtb {calls} (cached {hits}) | {clock}"

    def draw(force=False):
        if not tty:
            if force:
                print(status(), file=sys.stderr, flush=True)
            return
        now = time.time()
        if now - st["last"] < 0.5 and not force:
            return
        st["last"] = now
        w = max(shutil.get_terminal_size((80, 24)).columns - 1, 1)
        sys.stderr.write("\r" + status()[:w].ljust(w))
        sys.stderr.flush()

    def add(name, rec, dq, hsfrac):
        rec = dict(rec, params=dict(rec.get("params") or {}, cfg_tag=tags[name]))
        x = registry.make_id(rec["sampler"], rec["symbols"], rec["positions"], rec["charge"], 0, False,
                             rec["job"]).encode()
        if dedup:
            key = (name, int(dq), rec["composition"], rec["charge"], tuple(rec["symbols"]))
            pos = np.asarray(rec["positions"], float)
            idx = range(len(pos))
            if any(np.linalg.norm(sample._kabsch(pos, ref, idx) - ref, axis=1).max() < dedup
                   for ref in seen.get(key, ())):
                dedup_dropped[name] += 1
                return 0
            seen.setdefault(key, []).append(pos)
        u = int(hashlib.sha1(x).hexdigest()[:8], 16) / 2 ** 32
        if int(dq) != 0 and u >= redox:
            return 0
        u_hs = int(hashlib.sha1(x + b"hs").hexdigest()[:8], 16) / 2 ** 32

        def emit(q, params, parent=None):
            n, low_sid = 0, None
            fixed = rec.get("mults") is not None  # the sampler fixed the multiplicity (atom ground states)
            states = g16.spin_states(rec["symbols"], q, offsets=tuple(range(0, 20, 2)) if fixed else offsets,
                                     bs=rec["bs"])
            low = min((mult for mult, _ in states), default=0)
            for mult, bs in states:
                if fixed:
                    if mult not in rec["mults"]:
                        continue
                elif mult > low and not bs and u_hs >= hsfrac:
                    continue
                row = dict(rec, charge=q, mult=mult, bs=bs, params=params, parent_id=parent,
                           formula=Atoms(symbols=rec["symbols"]).get_chemical_formula())
                sid = registry.make_id(row["sampler"], row["symbols"], row["positions"], q, mult, bs, row["job"])
                if low_sid is None or mult == low:
                    low_sid = sid
                if registry.gjf_path(con, sid):
                    if parent is None:  # same deck reached as a vertical one first: the stage record wins
                        con.execute("UPDATE structures SET parent_id=NULL, params=? WHERE id=? AND parent_id IS NOT NULL",
                                    (json.dumps(params or {}), sid))
                    continue  # deck already written by an earlier run
                path = out / "gjf" / rec["job"] / name / f"{sid}.gjf"
                path.parent.mkdir(parents=True, exist_ok=True)
                tmp = path.with_name(path.name + ".tmp")
                tmp.write_text(g16.deck(rec["symbols"], rec["positions"], q, mult, bs,
                                        rec["job"], sid, nproc, mem, pop_always))
                os.replace(tmp, path)  # deck first: a crash never leaves a registry row without its deck
                row["gjf_path"] = str(path)
                _, is_new = registry.add(con, row)
                if not is_new:  # row from a crashed earlier run: fill in the deck path
                    registry.set_gjf(con, sid, path)
                added[name] += 1
                n += 1
            return n, low_sid

        q0 = int(rec["charge"])
        n, base = emit(q0, rec["params"])
        if (rec["job"] != "sp" or vertical <= 0.0 or rec["params"].get("kind") in ("cation", "atom")
                or rec["params"].get("contact") == "rcov_sum"):
            return n  # vertical decks only for geometries that depend on charge
        u_v = int(hashlib.sha1(x + b"vert").hexdigest()[:8], 16) / 2 ** 32
        if u_v >= vertical:
            return n
        symbols = rec["symbols"]
        formal = q0 - int(dq)
        if q0 - 1 >= max(qlo, -1) and abs(q0 - 1 - formal) <= 1 and any(s in ("Li", "Na") for s in symbols):
            n += emit(q0 - 1, dict(rec["params"], vertical_dq=-1), base)[0]
        if q0 + 1 <= qhi and abs(q0 + 1 - formal) <= 1 and any(s not in ("Li", "Na") for s in symbols):
            n += emit(q0 + 1, dict(rec["params"], vertical_dq=1), base)[0]
        return n

    mp = multiprocessing.get_context("fork")  # children inherit catalog, config and monkeypatched samplers
    q = mp.Queue()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            while todo or running:
                while todo and len(running) < nmax:
                    i, dq, name = todo.pop(0)
                    stage[i] = dict(dq=dq, name=name, u=0, k=plan[name], label="", new=0, calls=0, hits=0,
                                    t0=time.time())
                    running[i] = mp.Process(target=_child, args=(i, (cfg, cat, out, dq, name), q), daemon=True)
                    running[i].start()
                try:
                    i, kind, *rest = q.get(timeout=1)
                except queue.Empty:
                    for i, p in list(running.items()):
                        if p.exitcode not in (None, 0):  # killed (e.g. out of memory) before sending "end"
                            del running[i]
                            st["done"] += 1
                            s = stage[i]
                            errors.append(f"dq={float(s['dq']):+g} {s['name']}: process died, exit code {p.exitcode}")
                    draw()
                    continue
                s = stage[i]
                if kind == "unit":
                    s["u"] += 1
                    s["label"], s["calls"], s["hits"] = rest
                    draw(force=True)
                elif kind == "rec":
                    n = add(s["name"], rest[0], s["dq"], hsfrac[s["name"]])
                    s["new"] += n
                    st["new"] += n
                    draw()
                else:
                    w, s["calls"], s["hits"], err = rest
                    msgs |= Counter(w)  # max per message: stages rebuild the same clusters, warn once
                    running.pop(i).join()
                    st["done"] += 1
                    con.commit()
                    if err:
                        errors.append(f"dq={float(s['dq']):+g} {s['name']}:\n{err}")
                    if tty:
                        sys.stderr.write("\r" + " " * max(shutil.get_terminal_size((80, 24)).columns - 1, 1) + "\r")
                    print(f"[{i}/{n_stages}] dq={float(s['dq']):+g} {s['name']} {'FAILED' if err else 'done'}: "
                          f"{s['new']} new decks, xtb {s['calls']} (cached {s['hits']}), "
                          f"{time.time() - s['t0']:.1f} s", file=sys.stderr, flush=True)
            con.commit()
            manifest(con, out)
            for name in names:
                total = con.execute("SELECT COUNT(*) FROM structures WHERE sampler=?", (name,)).fetchone()[0]
                more = f", {dedup_dropped[name]} near-duplicates dropped" if dedup_dropped[name] else ""
                print(f"{name}: {added[name]} new, {total} total{more}")
        finally:
            for p in running.values():
                p.terminate()
            con.commit()
            xtb.CACHE = None
            _log_warnings(out, list(msgs.elements()) + [str(w.message) for w in caught])
    con.close()
    if errors:
        print("\n".join(errors), file=sys.stderr)
        raise SystemExit(f"{len(errors)} stage(s) failed; decks from the other stages are written")


def _log_warnings(out, caught):
    """warnings.log holds one warning per line; stderr gets a short summary grouped by kind."""
    from collections import Counter
    msgs = [re.sub(r"[ \t]+", " ", re.sub(r"\s*\n\s*", " | ", str(m))).strip() for m in caught]
    log = Path(out) / "warnings.log"
    log.write_text("".join(m + "\n" for m in msgs))
    if not caught:
        return

    def kind(m):
        k = re.sub(r"\[[^]]*\]", "[..]", m.split("(", 1)[0])
        return re.sub(r"[+-]?\d+(?:\.\d+)?", "#", k).strip()[:90]

    print(f"{len(caught)} warnings, details in {log}", file=sys.stderr)
    for k, n in Counter(map(kind, msgs)).most_common():
        print(f"  {n}  {k}", file=sys.stderr)


LINK0 = {"nproc": re.compile(r"^%nprocshared=(.*)$", re.M), "mem": re.compile(r"^%mem=(.*)$", re.M)}


def _link0(path):
    """(nproc, mem) as written in the deck itself, so the manifest always matches the files."""
    text = Path(path).read_text()
    return tuple((m[1].strip() if m else "") for m in (LINK0["nproc"].search(text), LINK0["mem"].search(text)))


def manifest(con, out):
    cols = ["id", "gjf_path", "sampler", "n_atoms", "charge", "mult", "bs", "job", "nproc", "mem"]
    rows = con.execute("SELECT " + ",".join(cols[:8]) + " FROM structures WHERE status='generated' "
                       "ORDER BY sampler, id").fetchall()
    missing = [r for r in rows if not r["gjf_path"]]
    if missing:
        warnings.warn(f"manifest: skipped {len(missing)} row(s) without a deck (gjf_path NULL)")
    with open(Path(out) / "manifest.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            if r["gjf_path"]:
                w.writerow([r[k] for k in cols[:8]] + list(_link0(r["gjf_path"])))


def resources_cmd(db, nproc=None, mem=None, samplers=None):
    """Rewrite %nprocshared / %mem of not-yet-run decks (status 'generated') in place, then the manifest."""
    if nproc is None and mem is None:
        raise SystemExit("give --nproc and/or --mem")
    con = registry.connect(db)
    sql = "SELECT gjf_path FROM structures WHERE status='generated'"
    if samplers:
        sql += f" AND sampler IN ({','.join('?' * len(samplers))})"
    n = 0
    for (path,) in con.execute(sql, list(samplers or [])).fetchall():
        p = Path(path)
        text = p.read_text()
        for key, val in (("nproc", nproc), ("mem", mem)):
            if val is not None:
                text, k = LINK0[key].subn(f"%{'nprocshared' if key == 'nproc' else 'mem'}={val}", text, count=1)
                if not k:
                    raise SystemExit(f"{p}: no %{key} line")
        tmp = p.with_suffix(".gjf.tmp")
        tmp.write_text(text)
        os.replace(tmp, p)  # atomic: a crash never leaves a half-written deck
        n += 1
    manifest(con, Path(db).parent)
    con.close()
    print(f"{n} decks updated")


def parse_cmd(paths, out, db=None, keep_all=False):
    logs = []
    for p in map(Path, paths):
        logs += (sorted(p.rglob("*.log")) + sorted(p.rglob("*.out"))) if p.is_dir() else [p]
    logs = list(dict.fromkeys(logs))
    from collections import Counter
    frames, failed, missing, incomplete, other, dropped, updates = [], 0, 0, 0, 0, Counter(), []
    con = registry.connect(db) if db else None
    runs = {}
    for log in logs:
        steps = parse.parse_log(log)
        normal = bool(steps) and all(s["normal"] for s in steps)
        err = any(s["error"] for s in steps)
        if not steps:  # nothing past the route echo: died early (bad %mem, ...), or not a Gaussian log (slurm-*.out)
            text = Path(log).read_text(errors="replace")
            err = "Error termination" in text
            if not err and "Entering Gaussian System" not in text:
                other += 1
                continue
        if not normal and not err:
            incomplete += 1  # still running or killed: no frames, registry status unchanged
            continue
        runs.setdefault(parse.log_id(log, steps), []).append((log, steps, normal))  # id from the deck title
    for sid, cands in runs.items():
        log, steps, normal = next((c for c in cands if c[2]), cands[-1])  # a normal-terminated copy wins
        if len(cands) > 1:
            print(f"warning: {sid}: {len(cands)} logs, using {log}", file=sys.stderr)
        fr = [a for i, st in enumerate(steps) for a in parse.to_atoms(st, sid, i)]
        row = con.execute("SELECT sampler, composition, params FROM structures WHERE id=?", (sid,)).fetchone() \
            if con else None
        try:
            prm = json.loads(row["params"] or "{}") if row else {}
        except ValueError:
            print(f"warning: {sid}: unreadable registry params; no path information", file=sys.stderr)
            prm = {}
        pk = parse.path_key(row["sampler"], prm) if row else None
        for a in fr:  # group = composition: split train/test by group so no cluster or path is on both sides
            if row:
                a.info.update(sampler=row["sampler"], composition=row["composition"], group=row["composition"])
            if pk and a.info.get("job") == "sp":
                a.info.update(path=pk[0], path_pos=pk[1], rattled="rattle_A" in prm)
        ok = normal and bool(fr)
        failed += 0 if ok else 1
        frames += fr
        updates.append((sid, "parsed" if ok else "failed", log))
    n_bad = parse.flag_paths(frames)  # before filtering: neighbours are the generated path, not the survivors
    kept = []
    for a in frames:
        why = "" if keep_all else parse.reject(a)
        dropped[why] += bool(why)
        if not why:
            kept.append(a)
    frames = kept
    tmp = Path(f"{out}.tmp")  # write first: a failed write must not leave 'parsed' rows without data
    if frames:
        ase_write(tmp, frames, format="extxyz")
    else:
        tmp.write_text("")  # an empty result replaces an older file: no stale frames
    os.replace(tmp, out)
    if con:  # one short write transaction after the extxyz: generate on the same registry is blocked only briefly
        for sid, status, log in updates:
            if not registry.set_status(con, sid, status, log):
                print(f"warning: {sid} is not in the registry", file=sys.stderr)
                missing += 1
        con.commit()
        con.close()
    print(f"{len(logs)} logs, {len(frames)} frames, {failed} failed"
          + (f", {incomplete} without a termination line skipped (running or killed)" if incomplete else "")
          + (f", {other} non-Gaussian files ignored" if other else ""))
    if n_bad:
        print(f"path check: {n_bad} neighbour pairs marked (info path_flag): energy change differs from the force "
              f"work by > {parse.PATH_JUMP_EV:g} eV (unrattled, step < {parse.PATH_MAX_STEP_A:g} A), or an NPA "
              f"atomic spin changes by > {parse.PATH_SPIN_E:g} e")
    if +dropped:
        print(f"dropped by quality filters: {dropped['force']} with a force > {parse.MAX_FORCE:g} eV/A, "
              f"{dropped['s2']} with <S2> - S(S+1) above {parse.MAX_S2_DEV_PURE:g} ({parse.MAX_S2_DEV:g} for BS "
              f"singlets, {parse.MAX_S2_DEV_TM:g} with Sc-Zn), {dropped['homo']} with HOMO > 0 eV")
    return missing


def _ids(targets):
    """ids from registry ids, .gjf files, directories of .gjf decks, or a manifest .csv."""
    ids = []
    for t in targets:
        p = Path(t)
        if p.is_dir():
            ids += [_gjf_id(g) for g in sorted(p.rglob("*.gjf"))]
        elif p.suffix == ".csv":
            with open(p, newline="") as f:
                ids += [r["id"] for r in csv.DictReader(f)]
        else:
            ids.append(_gjf_id(p) if p.suffix == ".gjf" else t)
    return ids


def _gjf_id(path):
    """Registry id of a deck: its stem, else its title line (a copy renamed to <id>/in.gjf)."""
    p = Path(path)
    if parse.ID_RE.fullmatch(p.stem):
        return p.stem
    return next((l.strip() for l in p.read_text(errors="replace").splitlines() if parse.ID_RE.fullmatch(l.strip())),
                p.stem)


def submit_cmd(db, targets, undo=False):
    """Mark decks submitted (generated or failed -> submitted: a rerun) so the manifest no longer lists them;
    undo sets submitted -> generated."""
    old, new = (("submitted",), "generated") if undo else (("generated", "failed"), "submitted")
    ids = _ids(targets)
    con = registry.connect(db)
    n = sum(con.execute(f"UPDATE structures SET status=? WHERE id=? AND status IN ({','.join('?' * len(old))})",
                        (new, i, *old)).rowcount for i in ids)
    con.commit()
    manifest(con, Path(db).parent)
    con.close()
    print(f"{n} of {len(ids)} decks {'/'.join(old)} -> {new}")


def list_cmd(db, status="submitted", samplers=None):
    con = registry.connect(db)
    q, args = "SELECT id, sampler, status, gjf_path, log_path FROM structures WHERE status=?", [status]
    if samplers:
        q += f" AND sampler IN ({','.join('?' * len(samplers))})"
        args += samplers
    rows = con.execute(q + " ORDER BY sampler, id", args).fetchall()
    con.close()
    print("id\tsampler\tstatus\tgjf_path\tlog_path")
    for r in rows:
        print("\t".join(str(r[k] or "") for k in ("id", "sampler", "status", "gjf_path", "log_path")))
    print(f"{len(rows)} {status}", file=sys.stderr)


def status_cmd(db):
    con = registry.connect(db)
    rows = registry.counts(con)
    print(f"{'sampler':<14}{'status':<12}{'n':>6}")
    for r in rows:
        print(f"{r['sampler']:<14}{r['status']:<12}{r['n']:>6}")
    con.close()


def export_cmd(db, out):
    con = registry.connect(db)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=registry.COLUMNS)
        w.writeheader()
        for r in con.execute("SELECT * FROM structures ORDER BY id"):
            w.writerow(dict(r))
    con.close()


def add_species_cmd(path, name, charge=None, kind=None, force=False):
    try:
        sp = species.add_user_species(path, name, charge=charge, kind=kind, force=force)
    except ValueError as exc:  # user input problem: message, not a traceback
        raise SystemExit(f"add-species: {exc}")
    mol = Chem.MolFromMolFile(sp.molfile, removeHs=False)
    formula = Atoms(symbols=[a.GetSymbol() for a in mol.GetAtoms()]).get_chemical_formula()
    for key, val in (("name", sp.name), ("formula", formula), ("charge", sp.charge),
                     ("kind", sp.kind), ("smiles", sp.smiles), ("path", sp.molfile)):
        print(f"{key}: {val}")


def species_cmd(kind=None, group=None):
    cat = species.catalog()
    rows = [sp for sp in cat.values()
            if (kind is None or sp.kind == kind) and (group is None or sp.group == group)]
    rows.sort(key=lambda sp: (sp.kind, sp.name))
    cols = ("name", "kind", "charge", "n_atoms", "group", "smiles")
    data = [{"name": sp.name, "kind": sp.kind, "charge": str(sp.charge),
             "n_atoms": str(build.n_atoms(sp)), "group": sp.group, "smiles": sp.smiles}
            for sp in rows]
    width = {c: max([len(c)] + [len(r[c]) for r in data]) for c in cols}
    print("  ".join(f"{c:<{width[c]}}" for c in cols))
    for r in data:
        print("  ".join(f"{r[c]:<{width[c]}}" for c in cols))
    print(f"{len(data)} species")


def _display():
    if not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")):
        raise SystemExit("no display found. On WSL2 run: export DISPLAY=:0 (WSLg), "
                         "over ssh use ssh -X, or use --save FILE.xyz")


def _gjf_info(path):
    """(id, charge, mult, bs) from a deck: title after the route's blank line, first 'q m' line, guess=mix."""
    lines = Path(path).read_text().splitlines()
    try:
        start = next(i for i, l in enumerate(lines) if l.startswith("#"))
    except StopIteration:
        raise ValueError(f"{path}: no route line") from None
    gap = next(i for i in range(start, len(lines)) if not lines[i].strip())
    title = next(l.strip() for l in lines[gap + 1:] if l.strip())
    m = re.search(r"^-?\d+ \d+$", "\n".join(lines), re.M)
    if m is None:
        raise SystemExit(f"{path}: no charge/mult line")
    q, mult = m.group().split()
    return title, int(q), int(mult), "guess=mix" in "\n".join(lines[:gap]).lower()


def view_cmd(targets, opt=False, charge=None, save=None):
    cat = species.catalog()
    if len(targets) == 1 and targets[0] in cat:
        return _view_species(cat[targets[0]], opt, charge, save)
    if opt or charge is not None:
        raise SystemExit("--opt and --charge are only valid for a species name")
    files = []
    for target in targets:
        p = Path(target)
        if p.is_dir():
            files += sorted(p.rglob("*.gjf"))
        elif p.is_file():
            files.append(p)
        else:
            raise SystemExit(f"unknown species or file: {target}")
    if not files:
        raise SystemExit("no .gjf files found")
    frames = []
    for i, path in enumerate(files):
        atoms = ase_read(path, format="gaussian-in")
        try:
            title, q, mult, bs = _gjf_info(path)
        except ValueError as exc:  # deck without a route line: message, not a traceback
            raise SystemExit(str(exc)) from None
        atoms.info.update(id=title, charge=q, spin=mult, bs=bs)
        print(f"{i}: id={title} charge={q} mult={mult} bs={bs} n_atoms={len(atoms)} {path}")
        frames.append(atoms)
    if save:
        ase_write(save, frames)
        return
    _display()
    from ase.visualize import view
    view(frames)


def _view_species(sp, opt=False, charge=None, save=None):
    charge = sp.charge if charge is None else int(charge)
    if sp.molfile:
        mol = Chem.MolFromMolFile(sp.molfile, removeHs=False)
        if mol is None:
            raise SystemExit(f"cannot read {sp.molfile}")
    else:
        mol = species.embed(sp.smiles)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    positions = np.zeros((1, 3)) if len(symbols) == 1 else mol.GetConformer().GetPositions()
    if opt:
        with tempfile.TemporaryDirectory() as d:
            r = xtb.run_xtb(symbols, positions, charge, g16.n_electrons(symbols, charge) % 2, d, mode="opt")
        if not r["ok"]:
            raise SystemExit(f"xtb opt failed: {r.get('err', '')[:200]}")
        positions = np.asarray(r["positions"], float)
        print(f"xtb energy: {r['E_eV']:.6f} eV")
    atoms = Atoms(symbols=symbols, positions=positions)
    if save:
        ase_write(save, atoms)
        return
    _display()
    from ase.visualize import view
    view(atoms)


def main(argv=None):
    p = argparse.ArgumentParser(prog="rmd")
    sub = p.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate", help="sample geometries and write Gaussian16 decks")
    g.add_argument("config")
    g.add_argument("--dry", action="store_true", help="print counts without running xtb")
    g.add_argument("--only", nargs="+", default=None)
    q = sub.add_parser("parse", help="parse Gaussian16 logs to extxyz")
    q.add_argument("paths", nargs="+")
    q.add_argument("-o", required=True)
    q.add_argument("--registry")
    q.add_argument("--keep-all", action="store_true", help="skip the OMol25 force / <S2> filters")
    r = sub.add_parser("resources", help="change %%nprocshared / %%mem of not-yet-run decks")
    r.add_argument("--registry", required=True)
    r.add_argument("--nproc", type=int)
    r.add_argument("--mem", help="e.g. 64GB")
    r.add_argument("--sampler", nargs="+", default=None)
    s = sub.add_parser("status", help="registry counts")
    s.add_argument("--registry", required=True)
    m = sub.add_parser("submit", help="mark decks submitted so the manifest no longer lists them")
    m.add_argument("targets", nargs="+", help="ids, .gjf files, directories of .gjf decks, or a manifest .csv")
    m.add_argument("--registry", required=True)
    m.add_argument("--undo", action="store_true", help="submitted -> generated (e.g. cancelled jobs)")
    ls = sub.add_parser("list", help="list decks with one status (default submitted)")
    ls.add_argument("--registry", required=True)
    ls.add_argument("--status", default="submitted", choices=("generated", "submitted", "parsed", "failed"))
    ls.add_argument("--sampler", nargs="+", default=None)
    e = sub.add_parser("export", help="registry rows to CSV")
    e.add_argument("--registry", required=True)
    e.add_argument("-o", required=True)
    u = sub.add_parser("add-species", help="add a user species from .mol/.mol2/.xyz/.pdb")
    u.add_argument("file")
    u.add_argument("--name", required=True)
    u.add_argument("--charge", type=int)
    u.add_argument("--kind", choices=("anion", "solvent"))
    u.add_argument("--force", action="store_true")
    sp = sub.add_parser("species", help="list the species catalog")
    sp.add_argument("--kind", choices=("cation", "anion", "solvent"))
    sp.add_argument("--group")
    v = sub.add_parser("view", help="show or save a species geometry or Gaussian16 decks")
    v.add_argument("targets", nargs="+",
                   help="species name, or .gjf files / directories of .gjf decks")
    v.add_argument("--opt", action="store_true")
    v.add_argument("--charge", type=int)
    v.add_argument("--save")
    a = p.parse_args(argv)
    if a.cmd != "generate" and getattr(a, "registry", None) and not Path(a.registry).is_file():
        raise SystemExit(f"{a.registry}: no such registry")  # sqlite would silently create an empty one
    if a.cmd == "generate":
        generate(yaml.safe_load(Path(a.config).read_text()), only=a.only, dry=a.dry)
    elif a.cmd == "parse":
        parse_cmd(a.paths, a.o, a.registry, a.keep_all)
    elif a.cmd == "status":
        status_cmd(a.registry)
    elif a.cmd == "submit":
        submit_cmd(a.registry, a.targets, a.undo)
    elif a.cmd == "list":
        list_cmd(a.registry, a.status, a.sampler)
    elif a.cmd == "resources":
        resources_cmd(a.registry, a.nproc, a.mem, a.sampler)
    elif a.cmd == "add-species":
        add_species_cmd(a.file, a.name, a.charge, a.kind, a.force)
    elif a.cmd == "species":
        species_cmd(a.kind, a.group)
    elif a.cmd == "view":
        view_cmd(a.targets, a.opt, a.charge, a.save)
    else:
        export_cmd(a.registry, a.o)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
