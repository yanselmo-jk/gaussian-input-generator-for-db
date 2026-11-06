from collections import Counter

import numpy as np

from rmlipdata import build, species


def _n_cations(cat, combo):
    return sum(cat[n].kind == "cation" for n in combo)


def _cfg(**kw):
    out = {"min_components": 1, "max_cations": 2}
    out.update(kw)
    return out


def test_dict_max_components():
    """max_components may be an int or a dict keyed by the number of cations."""
    cat = species.catalog()
    sp = {"cations": ["Li"], "anions": ["PF6"], "solvents": ["EC"]}
    out = build.combos(cat, sp, _cfg(max_components={0: 1, 1: 4, 2: 2}, max_combos=10 ** 6),
                       10 ** 6, np.random.default_rng(0))
    longest = {}
    for combo in out:
        nc = _n_cations(cat, combo)
        longest[nc] = max(longest.get(nc, 0), len(combo))
    assert longest == {0: 1, 1: 5, 2: 4}

    flat = build.combos(cat, sp, _cfg(max_components=2, max_combos=10 ** 6), 10 ** 6,
                        np.random.default_rng(0))
    assert all(len(c) - _n_cations(cat, c) <= 2 for c in flat)  # int: same component cap


def test_two_cations_need_anion():
    cat = species.catalog()
    sp = {"cations": ["Li", "Na"], "anions": ["PF6"], "solvents": ["EC"]}
    out = build.combos(cat, sp, _cfg(max_components={0: 2, 1: 2, 2: 3}, max_combos=10 ** 6),
                       10 ** 6, np.random.default_rng(0))
    multi = [c for c in out if _n_cations(cat, c) >= 2]
    assert multi
    assert all(any(cat[n].kind == "anion" for n in c) for c in multi)


def test_cation_weights():
    """Weights [0, 1, 0]: every combo has exactly one cation."""
    cat = species.catalog()
    sp = {"cations": ["Li", "Na"], "anions": ["PF6"], "solvents": ["EC"]}
    out = build.combos(cat, sp, _cfg(max_components={0: 2, 1: 2, 2: 2}, max_combos=50,
                                     cation_weights=[0, 1, 0]),
                       10 ** 6, np.random.default_rng(0))
    assert out and all(_n_cations(cat, c) == 1 for c in out)


def test_combos_deterministic():
    cat = species.catalog()
    sp = {"cations": ["Li", "Na"], "anions": ["PF6", "BF4"], "solvents": ["EC", "H2O"]}
    cfg = _cfg(max_components={0: 2, 1: 2, 2: 2}, max_combos=6, cation_weights=[0.2, 0.7, 0.1])
    a = build.combos(cat, sp, cfg, 60, np.random.default_rng(7))
    b = build.combos(cat, sp, cfg, 60, np.random.default_rng(7))
    assert a == b and len(a) == 6


def test_group_shortfall_redistributed():
    """Group 0 has 2 combos but asks for 6 of 12: its share goes to the higher-weight groups."""
    cat = species.catalog()
    sp = {"cations": ["Li"], "anions": ["PF6"], "solvents": ["EC"]}
    cfg = _cfg(max_components={0: 1, 1: 2, 2: 3}, max_combos=12, cation_weights=[0.5, 0.3, 0.2])
    out = build.combos(cat, sp, cfg, 10 ** 6, np.random.default_rng(0))
    assert len(out) == 12
    assert Counter(_n_cations(cat, c) for c in out) == {0: 2, 1: 5, 2: 5}
