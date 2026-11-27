"""audiopitchformant 核心 API（初始未实现，仅标准库）。"""


class PitchShift:
    def __init__(self, frame, hop, semitones, preserve_formant):
        """构造"""
        raise NotImplementedError

    def process(self, samples):
        """处理"""
        raise NotImplementedError

    def envelope(self, samples):
        """谱包络"""
        raise NotImplementedError
