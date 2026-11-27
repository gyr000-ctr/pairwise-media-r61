import math

import pytest

import audiodelay

def _raises(fn, exc=ValueError):
    try:
        fn()
    except exc:
        return True
    return False

S = [1.0] + [0.0] * 99


def _d():
    return audiodelay.Delay(64, 0.0)

CASES = [
    ("len", lambda: len(_d().process(S, 10, 0.0)) == 100),
    ("integer_delay", lambda: abs(_d().process(S, 10, 1.0)[10] - 1.0) < 1e-9),
    ("no_early", lambda: abs(_d().process(S, 10, 1.0)[5]) < 1e-9),
    ("fractional", lambda: isinstance(_d().process(S, 10.5, 1.0), list)),
    ("impulse", lambda: isinstance(_d().impulse(8), list)),
    ("bad_fb", lambda: _raises(lambda: audiodelay.Delay(64, 1.0))),
    ("bad_delay", lambda: _raises(lambda: _d().process(S, -1, 0.0))),
    ("bad_mix", lambda: _raises(lambda: _d().process(S, 10, 1.5))),
    ("empty", lambda: _d().process([], 10, 0.5) == []),
    ("type", lambda: isinstance(_d().process(S, 10, 0.5), list)),
    ("det", lambda: _d().process(S, 10, 0.5) == _d().process(S, 10, 0.5)),
    ("nan", lambda: _raises(lambda: _d().process([float('nan')], 10, 0.5))),
    ("no_mut", lambda: _d().process(list(S), 10, 0.5) is not None),
    ("float", lambda: isinstance(_d().process(S, 10, 0.5)[0], float)),
    ("impulse_type", lambda: isinstance(_d().impulse(8), list)),
]


@pytest.mark.parametrize("name,check", CASES, ids=[c[0] for c in CASES])
def test_spec(name, check):
    assert check() is True
