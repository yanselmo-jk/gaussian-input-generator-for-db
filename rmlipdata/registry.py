"""SQLite registry of generated structures: one row per Gaussian16 deck, id = file stem."""
import hashlib
import json
import sqlite3
import time

COLUMNS = ("id", "parent_id", "sampler", "composition", "formula", "n_atoms", "charge", "mult", "bs", "job",
           "params", "gjf_path", "log_path", "status", "created")
SCHEMA = """CREATE TABLE IF NOT EXISTS structures (
  id TEXT PRIMARY KEY, parent_id TEXT, sampler TEXT NOT NULL, composition TEXT, formula TEXT,
  n_atoms INTEGER, charge INTEGER, mult INTEGER, bs INTEGER, job TEXT, params TEXT,
  gjf_path TEXT, log_path TEXT, status TEXT NOT NULL DEFAULT 'generated', created TEXT)"""


def connect(path):
    con = sqlite3.connect(path)
    con.row_factory = sqlite3.Row
    con.execute(SCHEMA)
    return con


def make_id(sampler, symbols, positions, charge, mult, bs, job):
    """sampler + 12 hex of sha1 over symbols, coordinates rounded to 1e-3 A, charge, mult, bs, job."""
    key = json.dumps([list(symbols), [[round(float(x), 3) + 0.0 for x in p] for p in positions],
                      int(charge), int(mult), bool(bs), job])
    return f"{sampler}-{hashlib.sha1(key.encode()).hexdigest()[:12]}"


def add(con, rec):
    """Insert rec (dict with COLUMNS minus id/created/status, plus symbols/positions); returns (id, is_new)."""
    sid = make_id(rec["sampler"], rec["symbols"], rec["positions"], rec["charge"], rec["mult"], rec["bs"],
                  rec["job"])
    row = {k: rec.get(k) for k in COLUMNS}
    row.update(id=sid, n_atoms=len(rec["symbols"]), bs=int(bool(rec["bs"])), status="generated",
               params=json.dumps(rec.get("params") or {}), created=time.strftime("%Y-%m-%dT%H:%M:%S"))
    cur = con.execute(f"INSERT OR IGNORE INTO structures ({','.join(COLUMNS)}) VALUES "
                      f"({','.join('?' * len(COLUMNS))})", [row[k] for k in COLUMNS])
    return sid, cur.rowcount == 1


def gjf_path(con, sid):
    """Deck path of sid, or None if the row is missing or has no deck yet."""
    row = con.execute("SELECT gjf_path FROM structures WHERE id=?", (sid,)).fetchone()
    return row["gjf_path"] if row is not None else None


def set_gjf(con, sid, path):
    con.execute("UPDATE structures SET gjf_path=? WHERE id=?", (str(path), sid))


def set_status(con, sid, status, log_path=None):
    """Returns False if sid is not in the registry."""
    cur = con.execute("UPDATE structures SET status=?, log_path=COALESCE(?, log_path) WHERE id=?",
                      (status, None if log_path is None else str(log_path), sid))
    return cur.rowcount == 1


def counts(con):
    return con.execute("SELECT sampler, status, COUNT(*) AS n FROM structures GROUP BY sampler, status "
                       "ORDER BY sampler, status").fetchall()
