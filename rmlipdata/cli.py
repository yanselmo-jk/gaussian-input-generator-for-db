"""rmd command line: generate Gaussian16 decks from xTB-sampled geometries, parse logs, inspect the registry."""
import argparse
import csv
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
    npairs = len([(a, b) for i, a in enumerate(sp) for b in sp[i:]
                  if not (cat[a].kind == "cation" and cat[b].kind == "cation")])
    plan = {}
    for name in names:
        scfg = cfg["samplers"][name]
        if name == "monomer":
            plan[name] = len(sp)
        elif name == "pair_scan":
            plan[name] = min(npairs, int(scfg.get("max_pairs", npairs)))
        elif name == "cluster":
            plan[name] = len([c for c in combos if len(c) >= 2])
        else:
            plan[name] = len(combos)
    return plan


def _dry(cfg, cat, names):
    for name, k in _plan(cfg, cat, names).items():
        print(f"{name}: {k}")


def generate(cfg, only=None, dry=False):
    cat = species.catalog()
    names = [n for n, sc in cfg["samplers"].items()
             if _enabled(sc) and (not only or n in only)]
    if dry:
        _dry(cfg, cat, names)
        return
    if shutil.which("xtb") is None:  # without xtb every structure would silently stay at its start geometry
        raise SystemExit(xtb.XTB_MISSING)
    out = Path(cfg["out"])
    out.mkdir(parents=True, exist_ok=True)
    threads = cfg.get("xtb", {}).get("threads", 4)
    os.environ["RMLIPDATA_XTB_THREADS"] = str(threads)
    con = registry.connect(out / "registry.sqlite")
    ctx = _ctx(cfg, cat, out)
    ctx.combos = _combos(cfg, cat)
    nproc = int(cfg.get("g16", {}).get("nproc", 16))
    mem = cfg.get("g16", {}).get("mem", "32GB")
    offsets = tuple(cfg.get("charge", {}).get("mult_offsets", (0, 2)))
    added = {n: 0 for n in names}

    dq_list = cfg.get("charge", {}).get("delta_q", [0])
    plan = _plan(cfg, cat, names)
    n_stages = len(dq_list) * len(names)
    for name in names:
        print(f"  {name}: {plan[name]} units x {len(dq_list)} delta_q", file=sys.stderr)
    print(f"{n_stages} stages (delta_q x sampler), xtb threads {threads}", file=sys.stderr)

    tty = sys.stderr.isatty()
    start = time.time()
    calls0 = xtb.CALLS
    st = {"i": 0, "dq": 0, "name": "", "u": 0, "k": 0, "label": "", "new": 0, "last": 0.0, "shown": None}

    def status():
        s = int(time.time() - start)
        h, s = divmod(s, 3600)
        mm, ss = divmod(s, 60)
        clock = f"{h}:{mm:02d}:{ss:02d}" if h else f"{mm:02d}:{ss:02d}"
        return (f"[{st['i']}/{n_stages}] dq={float(st['dq']):+g} {st['name']} {st['u']}/{st['k']} {st['label']}"
                f" | new {st['new']} | xtb {xtb.CALLS - calls0} | {clock}")

    def draw():
        if tty:
            now = time.time()
            if now - st["last"] < 0.5:
                return
            st["last"] = now
            w = max(shutil.get_terminal_size((80, 24)).columns - 1, 1)
            sys.stderr.write("\r" + status()[:w].ljust(w))
            sys.stderr.flush()
        elif st["shown"] != (st["i"], st["u"]):
            st["shown"] = (st["i"], st["u"])
            print(status(), file=sys.stderr, flush=True)

    def progress(label):
        st["u"] += 1
        st["label"] = label
        draw()

    ctx.progress = progress
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            i = 0
            for dq in dq_list:
                for name in names:
                    i += 1
                    scfg = cfg["samplers"][name]
                    st.update(i=i, dq=dq, name=name, u=0, k=plan[name], label="", new=0, last=0.0, shown=None)
                    stage_start = time.time()
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
                                st["new"] += 1
                        draw()
                    con.commit()
                    if tty:
                        sys.stderr.write("\r" + " " * max(shutil.get_terminal_size((80, 24)).columns - 1, 1) + "\r")
                    print(f"[{i}/{n_stages}] dq={float(dq):+g} {name} done: {st['new']} new decks, "
                          f"{time.time() - stage_start:.1f} s", file=sys.stderr, flush=True)
            manifest(con, out)
            for name in names:
                total = con.execute("SELECT COUNT(*) FROM structures WHERE sampler=?", (name,)).fetchone()[0]
                print(f"{name}: {added[name]} new, {total} total")
        finally:
            _log_warnings(out, caught)
    con.close()


def _log_warnings(out, caught):
    """warnings.log holds one warning per line; stderr gets a short summary grouped by kind."""
    from collections import Counter
    msgs = [re.sub(r"[ \t]+", " ", re.sub(r"\s*\n\s*", " | ", str(w.message))).strip() for w in caught]
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
    start = next(i for i, l in enumerate(lines) if l.startswith("#"))
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
        title, q, mult, bs = _gjf_info(path)
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
    r = sub.add_parser("resources", help="change %%nprocshared / %%mem of not-yet-run decks")
    r.add_argument("--registry", required=True)
    r.add_argument("--nproc", type=int)
    r.add_argument("--mem", help="e.g. 64GB")
    r.add_argument("--sampler", nargs="+", default=None)
    s = sub.add_parser("status", help="registry counts")
    s.add_argument("--registry", required=True)
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
    if a.cmd == "generate":
        generate(yaml.safe_load(Path(a.config).read_text()), only=a.only, dry=a.dry)
    elif a.cmd == "parse":
        parse_cmd(a.paths, a.o, a.registry)
    elif a.cmd == "status":
        status_cmd(a.registry)
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
