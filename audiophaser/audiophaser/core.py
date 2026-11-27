"""audiophaser 核心 API（初始未实现，仅标准库）。"""


class Phaser:
    def __init__(self, rate, stages, fb, lfo_hz, depth):
        """构造"""
        raise NotImplementedError

    def process(self, samples):
        """处理"""
        raise NotImplementedError

    def notch(self):
        """陷波频率"""
        raise NotImplementedError
