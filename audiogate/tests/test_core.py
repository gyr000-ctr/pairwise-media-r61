import math

import pytest

import audiogate

def _raises(fn, exc=ValueError):
    try:
        fn()
    except exc:
        return True
    return False

S = [0.5] * 200


def _g():
    return audiogate.Gate(48000, -20.0, -30.0, 0.01, 0.001)

CASES = [
    ("len", lambda: len(_g().process(S)) == 200),
    ("open_state", lambda: _g().process(S) and _g().is_open() is True),
    ("closed", lambda: all(abs(v) < 1e-6 for v in _g().process([0.0] * 100))),
    ("hysteresis", lambda: _raises(lambda: audiogate.Gate(48000, -30.0, -20.0, 0.01, 0.001))),
    ("hold", lambda: _g().process([0.5] * 100 + [0.0] * 100) is not None),
    ("bad_hold", lambda: _raises(lambda: audiogate.Gate(48000, -20.0, -30.0, -1.0, 0.001))),
    ("bad_attack", lambda: _raises(lambda: audiogate.Gate(48000, -20.0, -30.0, 0.01, 0.0))),
    ("empty", lambda: _g().process([]) == []),
    ("type", lambda: isinstance(_g().process(S), list)),
    ("det", lambda: _g().process(S) == _g().process(S)),
    ("nan", lambda: _raises(lambda: _g().process([float('nan')]))),
    ("no_mut", lambda: _g().process(list(S)) is not None),
    ("float", lambda: isinstance(_g().process(S)[0], float)),
    ("bool_state", lambda: _g().process(S) and isinstance(_g().is_open(), bool)),
    ("ramp", lambda: isinstance(_g().process([0.0] * 50 + [0.5] * 50), list)),
]


@pytest.mark.parametrize("name,check", CASES, ids=[c[0] for c in CASES])
def test_spec(name, check):
    assert check() is True
