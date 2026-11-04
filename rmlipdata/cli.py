"""rmd command line: generate Gaussian16 decks from xTB-sampled geometries, parse logs, inspect the registry."""
import argparse
import csv
import os
import re
import sys
import zlib
from pathlib import Path

import numpy as np
import yaml
from ase import Atoms
from ase.io import write as ase_write

from . import build, g16, parse, registry, sample, species


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


def _dry(cfg, cat, names):
    combos = _combos(cfg, cat)
    sp = [n for k in ("cations", "anions", "solvents") for n in cfg["species"].get(k, [])]
    npairs = len([(a, b) for i, a in enumerate(sp) for b in sp[i:]
                  if not (cat[a].kind == "cation" and cat[b].kind == "cation")])
    for name in names:
        scfg = cfg["samplers"][name]
        if name == "monomer":
            k = len(sp)
        elif name == "pair_scan":
            k = min(npairs, int(scfg.get("max_pairs", npairs)))
        elif name == "cluster":
            k = len([c for c in combos if len(c) >= 2])
        else:
            k = len(combos)
        print(f"{name}: {k}")


def generate(cfg, only=None, dry=False):
    cat = species.catalog()
    names = [n for n, sc in cfg["samplers"].items()
             if _enabled(sc) and (not only or n in only)]
    if dry:
        _dry(cfg, cat, names)
        return
    out = Path(cfg["out"])
    out.mkdir(parents=True, exist_ok=True)
    os.environ["RMLIPDATA_XTB_THREADS"] = str(cfg.get("xtb", {}).get("threads", 4))
    con = registry.connect(out / "registry.sqlite")
    ctx = _ctx(cfg, cat, out)
    ctx.combos = _combos(cfg, cat)
    nproc = int(cfg.get("g16", {}).get("nproc", 16))
    mem = cfg.get("g16", {}).get("mem", "32GB")
    offsets = tuple(cfg.get("charge", {}).get("mult_offsets", (0, 2)))
    added = {n: 0 for n in names}
    for dq in cfg.get("charge", {}).get("delta_q", [0]):
        for name in names:
            scfg = cfg["samplers"][name]
            ctx.rng = np.random.default_rng(int(cfg.get("seed", 0)) + zlib.crc32(name.encode()))
            for rec in sample.SAMPLERS[name](ctx, scfg, dq):
                if len(rec["symbols"]) > _max_atoms(cfg):
                    continue
                for mult, bs in g16.spin_states(rec["symbols"], rec["charge"], offsets=offsets, bs=rec["bs"]):
                    row = dict(rec, mult=mult, bs=bs,
                               formula=Atoms(symbols=rec["symbols"]).get_chemical_formula())
                    sid, is_new = registry.add(con, row)
                    if is_new:
                        path = out / "gjf" / name / f"{sid}.gjf"
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_text(g16.deck(rec["symbols"], rec["positions"], rec["charge"], mult, bs,
                                                 rec["job"], sid, nproc, mem))
                        registry.set_gjf(con, sid, path)
                        added[name] += 1
            con.commit()
    manifest(con, out)
    for name in names:
        total = con.execute("SELECT COUNT(*) FROM structures WHERE sampler=?", (name,)).fetchone()[0]
        print(f"{name}: {added[name]} new, {total} total")
    con.close()


LINK0 = {"nproc": re.compile(r"^%nprocshared=(.*)$", re.M), "mem": re.compile(r"^%mem=(.*)$", re.M)}


def _link0(path):
    """(nproc, mem) as written in the deck itself, so the manifest always matches the files."""
    text = Path(path).read_text()
    return tuple((m[1].strip() if m else "") for m in (LINK0["nproc"].search(text), LINK0["mem"].search(text)))


def manifest(con, out):
    cols = ["id", "gjf_path", "sampler", "n_atoms", "charge", "mult", "bs", "job", "nproc", "mem"]
    rows = con.execute("SELECT " + ",".join(cols[:8]) + " FROM structures WHERE status='generated' "
                       "ORDER BY sampler, id").fetchall()
    with open(Path(out) / "manifest.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
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


def parse_cmd(paths, out, db=None):
    logs = []
    for p in map(Path, paths):
        logs += (sorted(p.glob("*.log")) + sorted(p.glob("*.out"))) if p.is_dir() else [p]
    frames, failed, missing = [], 0, 0
    con = registry.connect(db) if db else None
    for log in logs:
        steps = parse.parse_log(log)
        fr = [a for i, st in enumerate(steps) for a in parse.to_atoms(st, log.stem, i)]
        ok = bool(fr) and bool(steps) and all(s["normal"] for s in steps)
        failed += 0 if ok else 1
        frames += fr
        if con and not registry.set_status(con, log.stem, "parsed" if ok else "failed", log):
            print(f"warning: {log.stem} is not in the registry", file=sys.stderr)
            missing += 1
    if frames:  # write first: a failed write must not leave 'parsed' rows without data
        tmp = Path(f"{out}.tmp")
        ase_write(tmp, frames, format="extxyz")
        os.replace(tmp, out)
    if con:
        con.commit()
        con.close()
    print(f"{len(logs)} logs, {len(frames)} frames, {failed} failed")
    return missing


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
    r = sub.add_parser("resources", help="change %nprocshared / %mem of not-yet-run decks")
    r.add_argument("--registry", required=True)
    r.add_argument("--nproc", type=int)
    r.add_argument("--mem", help="e.g. 64GB")
    r.add_argument("--sampler", nargs="+", default=None)
    s = sub.add_parser("status", help="registry counts")
    s.add_argument("--registry", required=True)
    e = sub.add_parser("export", help="registry rows to CSV")
    e.add_argument("--registry", required=True)
    e.add_argument("-o", required=True)
    a = p.parse_args(argv)
    if a.cmd == "generate":
        generate(yaml.safe_load(Path(a.config).read_text()), only=a.only, dry=a.dry)
    elif a.cmd == "parse":
        parse_cmd(a.paths, a.o, a.registry)
    elif a.cmd == "status":
        status_cmd(a.registry)
    elif a.cmd == "resources":
        resources_cmd(a.registry, a.nproc, a.mem, a.sampler)
    else:
        export_cmd(a.registry, a.o)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
