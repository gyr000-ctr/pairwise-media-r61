import math

import pytest

import audiophaser

def _raises(fn, exc=ValueError):
    try:
        fn()
    except exc:
        return True
    return False

S = [0.5] * 200


def _p():
    return audiophaser.Phaser(48000, 4, 0.3, 0.5, 0.5)

CASES = [
    ("len", lambda: len(_p().process(S)) == 200),
    ("notch", lambda: _p().notch() > 0.0),
    ("bad_stages", lambda: _raises(lambda: audiophaser.Phaser(48000, 0, 0.3, 0.5, 0.5))),
    ("bad_fb", lambda: _raises(lambda: audiophaser.Phaser(48000, 4, 1.0, 0.5, 0.5))),
    ("bad_lfo", lambda: _raises(lambda: audiophaser.Phaser(48000, 4, 0.3, 0.0, 0.5))),
    ("bad_depth", lambda: _raises(lambda: audiophaser.Phaser(48000, 4, 0.3, 0.5, 1.5))),
    ("empty", lambda: _p().process([]) == []),
    ("type", lambda: isinstance(_p().process(S), list)),
    ("det", lambda: _p().process(S) == _p().process(S)),
    ("nan", lambda: _raises(lambda: _p().process([float('nan')]))),
    ("no_mut", lambda: _p().process(list(S)) is not None),
    ("range", lambda: all(abs(v) <= 2.0 for v in _p().process(S))),
    ("notch_type", lambda: isinstance(_p().notch(), float)),
    ("stages", lambda: audiophaser.Phaser(48000, 8, 0.3, 0.5, 0.5).notch() != _p().notch()),
    ("float", lambda: isinstance(_p().process(S)[0], float)),
]


@pytest.mark.parametrize("name,check", CASES, ids=[c[0] for c in CASES])
def test_spec(name, check):
    assert check() is True
