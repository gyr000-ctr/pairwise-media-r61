"""audiogate 核心 API（初始未实现，仅标准库）。"""


class Gate:
    def __init__(self, rate, open_db, close_db, hold, attack):
        """构造"""
        raise NotImplementedError

    def process(self, samples):
        """处理"""
        raise NotImplementedError

    def is_open(self):
        """当前状态"""
        raise NotImplementedError
