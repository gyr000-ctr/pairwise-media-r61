"""audiodelay 核心 API（初始未实现，仅标准库）。"""


class Delay:
    def __init__(self, max_delay, fb):
        """构造"""
        raise NotImplementedError

    def process(self, samples, delay, mix):
        """处理"""
        raise NotImplementedError

    def impulse(self, delay):
        """脉冲响应"""
        raise NotImplementedError
