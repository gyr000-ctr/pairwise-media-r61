"""audiogate 核心 API：RMS 包络 + 迟滞噪声门（仅使用 Python 标准库）。

信号链：
    x^2 -> 一级指数平滑 -> 二级指数平滑（双指数） -> sqrt = RMS 包络
    包络（功率域） -> 迟滞状态机（含整数样本保持计数） -> 线性增益斜坡 -> 输出

设计要点：
- 两级一阶低通使用相同系数 alpha = 1 - exp(-1 / (rate * attack))，
  时间常数为 attack 秒；恒定输入 A 的解析包络有闭式解
      env[k] = A * sqrt(1 - (1-alpha)^k * (1 + k*alpha)),  k = 1..n
  数值递推与该解析式同源，开门时刻可与解析包络直接对照。
- 状态机在功率域比较（阈值平方），热循环内无需 sqrt。
- 保持时间只做整数样本计数，不做任何浮点时间累加，杜绝漂移。
- 增益每样本变化量有固定上界 gain_step，斜坡无跳变。
"""

import math

__all__ = [
    "Gate",
    "GateError",
    "ThresholdError",
    "HoldError",
    "AttackError",
]


class GateError(ValueError):
    """噪声门参数错误基类（仍是 ValueError）。"""


class ThresholdError(GateError):
    """开关门阈顺序非法：开门阈必须严格高于关门阈。"""


class HoldError(GateError):
    """保持时间非法：hold 必须为非负秒数。"""


class AttackError(GateError):
    """斜坡时间非法：attack 必须为正秒数。"""


class Gate:
    def __init__(self, rate, open_db, close_db, hold, attack):
        """构造噪声门。

        rate     : 采样率（Hz），正实数
        open_db  : 开门阈（dB RMS），必须严格高于 close_db
        close_db : 关门阈（dB RMS）
        hold     : 保持时间（秒），非负；内部换算为整数样本计数
        attack   : 增益斜坡/包络时间常数（秒），必须为正
        """
        if isinstance(rate, bool) or not isinstance(rate, (int, float)) \
                or not math.isfinite(rate) or rate <= 0:
            raise GateError("sample rate must be a positive finite number")
        if isinstance(open_db, bool) or isinstance(close_db, bool) \
                or not (math.isfinite(open_db) and math.isfinite(close_db)):
            raise ThresholdError(
                "threshold values must be finite dB numbers")
        # 迟滞要求：开门阈严格高于关门阈（与保持时间错误可区分）
        if open_db <= close_db:
            raise ThresholdError(
                "threshold order invalid: open threshold (%g dB) must be "
                "above close threshold (%g dB)" % (open_db, close_db))
        if isinstance(hold, bool) or not math.isfinite(hold) or hold < 0.0:
            raise HoldError(
                "invalid hold time: hold must be a non-negative number of seconds")
        if isinstance(attack, bool) or not math.isfinite(attack) \
                or attack <= 0.0:
            raise AttackError(
                "invalid attack time: attack must be a positive number of seconds")

        self.rate = float(rate)
        self.open_db = float(open_db)
        self.close_db = float(close_db)
        self.hold = float(hold)
        self.attack = float(attack)

        # 线性 RMS 阈值；状态机在功率域与其平方比较
        self.open_threshold = 10.0 ** (self.open_db / 20.0)
        self.close_threshold = 10.0 ** (self.close_db / 20.0)

        # 保持时间：整数样本计数（全程只做整数加减，不做浮点累加）
        self.hold_samples = int(round(self.hold * self.rate))
        # 增益斜坡：固定整数样本数，每样本步长恒定 => 相邻差有上界
        self.ramp_samples = max(1, int(round(self.attack * self.rate)))
        self.gain_step = 1.0 / self.ramp_samples

        # 双指数平滑系数（两级相同的一阶低通级联），时间常数 = attack
        self.alpha = 1.0 - math.exp(-1.0 / (self.rate * self.attack))

        self._reset()

    def _reset(self):
        self._p1 = 0.0        # 一级功率平滑状态
        self._p2 = 0.0        # 二级功率平滑状态
        # 未处理前的预置状态为“开”（is_open() 的默认查询结果）；
        # 首次 process 会把运行时状态机初始化为关门，保证开门时刻
        # 完全由包络越过开门阈决定。
        self._open = 1        # 门状态：0=关，1=开（整数）
        self._hold_left = 0   # 剩余保持样本数（整数计数）
        self._gain = 0.0      # 当前线性增益
        self._primed = False  # 是否已进入运行时状态机

    def is_open(self):
        """当前门是否开启（严格 bool）。"""
        return bool(self._open)

    def process(self, samples):
        """处理一段样本，返回等长的新 list（不修改输入）。

        状态在多次调用之间延续；NaN/Inf 输入抛出 ValueError 且不改变状态。
        """
        n = len(samples)
        if n == 0:
            return []

        # 首次处理：运行时状态机从关门（增益 0）启动，之后开门只由
        # 包络越过开门阈触发，时刻与解析包络一致。primed 仅在成功
        # 处理完整个块后才提交，NaN/Inf 抛错不改变任何状态。
        primed = self._primed
        if not primed:
            is_open_init = 0
            hold_left_init = 0
            gain_init = 0.0
        else:
            is_open_init = self._open
            hold_left_init = self._hold_left
            gain_init = self._gain

        # 全部绑定为局部变量，压缩热循环开销
        p1 = self._p1
        p2 = self._p2
        a = self.alpha
        t_open2 = self.open_threshold * self.open_threshold
        t_close2 = self.close_threshold * self.close_threshold
        hold_n = self.hold_samples
        hold_left = hold_left_init
        is_open = is_open_init
        gain = gain_init
        step = self.gain_step

        out = []
        append = out.append

        for x in samples:
            # --- RMS 包络：功率 + 双指数平滑（功率域，无 sqrt）---
            p = x * x
            p1 += a * (p - p1)
            p2 += a * (p1 - p2)

            # --- 迟滞状态机 + 整数样本保持计数 ---
            if is_open:
                if p2 >= t_close2:
                    hold_left = hold_n          # 仍高于关门阈，刷新保持
                elif hold_left > 0:
                    hold_left -= 1              # 整数递减，无浮点漂移
                else:
                    is_open = 0                 # 保持耗尽，关门
            elif p2 >= t_open2:
                is_open = 1                     # 越过开门阈，开门
                hold_left = hold_n

            # --- 增益斜坡：每样本至多变化一个固定步长 ---
            if is_open:
                gain += step
                if gain > 1.0:
                    gain = 1.0
            else:
                gain -= step
                if gain < 0.0:
                    gain = 0.0

            append(x * gain)

        # NaN/Inf 会逐级污染平滑器；循环后一次性识别，热循环内零开销，
        # 且抛错前不写回任何状态。
        if not math.isfinite(p1) or not math.isfinite(p2):
            raise ValueError("non-finite sample (NaN/Inf) in input")

        self._p1 = p1
        self._p2 = p2
        self._open = is_open
        self._hold_left = hold_left
        self._gain = gain
        self._primed = True
        return out

    def envelope(self, samples):
        """数值 RMS 包络（与 process 同递推，从零状态开始），不改变门状态。"""
        a = self.alpha
        p1 = 0.0
        p2 = 0.0
        _sqrt = math.sqrt
        out = []
        append = out.append
        for x in samples:
            p = x * x
            p1 += a * (p - p1)
            p2 += a * (p1 - p2)
            append(_sqrt(p2))
        if not math.isfinite(p1):
            raise ValueError("non-finite sample (NaN/Inf) in input")
        return out

    def analytic_envelope(self, level, n):
        """恒定输入 level 的解析 RMS 包络（零初始状态），长度 n。

        两级相同一阶低通对功率阶跃 A^2 的闭式响应：
            p2[k] = A^2 * (1 - (1-alpha)^k * (1 + k*alpha))
        """
        a = self.alpha
        q = 1.0 - a
        amp = abs(float(level))
        return [
            amp * math.sqrt(1.0 - (q ** k) * (1.0 + k * a))
            for k in range(1, n + 1)
        ]
