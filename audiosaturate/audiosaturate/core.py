"""audiosaturate 核心 API（初始未实现，仅标准库）。"""


class Saturate:
    def __init__(self, drive, oversample):
        """构造"""
        raise NotImplementedError

    def shape(self, x):
        """曲线"""
        raise NotImplementedError

    def process(self, samples):
        """处理"""
        raise NotImplementedError
