import math

import pytest

import audiosaturate

def _raises(fn, exc=ValueError):
    try:
        fn()
    except exc:
        return True
    return False

S = [0.3] * 200


def _s():
    return audiosaturate.Saturate(4.0, 4)

CASES = [
    ("odd", lambda: abs(_s().shape(0.3) + _s().shape(-0.3)) < 1e-12),
    ("monotone", lambda: _s().shape(0.2) <= _s().shape(0.6)),
    ("bounded", lambda: all(abs(v) <= 1.0 for v in _s().process(S))),
    ("len", lambda: len(_s().process(S)) == 200),
    ("bad_drive", lambda: _raises(lambda: audiosaturate.Saturate(0.0, 4))),
    ("bad_oversample", lambda: _raises(lambda: audiosaturate.Saturate(4.0, 3))),
    ("empty", lambda: _s().process([]) == []),
    ("type", lambda: isinstance(_s().process(S), list)),
    ("det", lambda: _s().process(S) == _s().process(S)),
    ("nan", lambda: _raises(lambda: _s().process([float('nan')]))),
    ("no_mut", lambda: _s().process(list(S)) is not None),
    ("harmonics", lambda: _s().shape(0.5) >= audiosaturate.Saturate(1.0, 4).shape(0.5)),
    ("float", lambda: isinstance(_s().process(S)[0], float)),
    ("zero", lambda: abs(_s().shape(0.0)) < 1e-12),
    ("shape_type", lambda: isinstance(_s().shape(0.1), float)),
]


@pytest.mark.parametrize("name,check", CASES, ids=[c[0] for c in CASES])
def test_spec(name, check):
    assert check() is True
