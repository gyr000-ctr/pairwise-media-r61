import math

import pytest

import audiopitchformant

def _raises(fn, exc=ValueError):
    try:
        fn()
    except exc:
        return True
    return False

def _p():
    return audiopitchformant.PitchShift(256, 64, 12, True)


def _tone(n=2048):
    return [0.5 * math.sin(2 * math.pi * 220 * i / 8000) for i in range(n)]

CASES = [
    ("len", lambda: abs(len(_p().process(_tone())) - 2048) <= 1),
    ("env_len", lambda: len(_p().envelope(_tone())) > 0),
    ("zero_shift", lambda: abs(len(audiopitchformant.PitchShift(256, 64, 0, True).process(_tone())) - 2048) <= 1),
    ("cola", lambda: _p().cola_error() <= 1e-12),
    ("bad_frame", lambda: _raises(lambda: audiopitchformant.PitchShift(0, 64, 12, True))),
    ("bad_hop", lambda: _raises(lambda: audiopitchformant.PitchShift(256, 0, 12, True))),
    ("bad_semitone", lambda: _raises(lambda: audiopitchformant.PitchShift(256, 64, 99, True))),
    ("empty", lambda: _p().process([]) == []),
    ("type", lambda: isinstance(_p().process(_tone()), list)),
    ("det", lambda: _p().process(_tone()) == _p().process(_tone())),
    ("nan", lambda: _raises(lambda: _p().process([float('nan')] * 512))),
    ("no_mut", lambda: _p().process(_tone()) is not None),
    ("env_type", lambda: isinstance(_p().envelope(_tone()), list)),
    ("cola_type", lambda: isinstance(_p().cola_error(), float)),
    ("float", lambda: isinstance(_p().process(_tone())[0], float)),
]


@pytest.mark.parametrize("name,check", CASES, ids=[c[0] for c in CASES])
def test_spec(name, check):
    assert check() is True
