"""audiophaser：全通级联相位器 + LFO（仅使用 Python 标准库）。

信号模型
--------
每一级是一阶全通滤波器

    H(z) = (a + z^-1) / (1 + a z^-1),        |a| < 1

在单位圆上 ``|H(e^{jw})| == 1``，即幅度恒为 1；其相位为

    phi(w) = -2 * atan( tan(w/2) / t0 ),      t0 = (1 + a) / (1 - a) > 0

对 w in (0, pi) 单调下降（0 -> -pi）。系数与目标中心频率 f0 的关系为

    a = (t0 - 1) / (t0 + 1),                 t0 = tan(pi * f0 / fs)

级联 N 级后总相位为 ``N * phi(w)``。与干信号等增益混合后

    M(w) = (1 + exp(j*N*phi(w))) / 2

在 ``N * phi(w) = -(2k-1) * pi`` 处形成陷波；取 k = 1 得到首个（最低）
陷波频率的解析解：

    f_notch = fs/pi * atan( t0 * tan(pi / (2N)) )

LFO 以固定相位增量 ``2*pi*lfo_hz/fs`` 逐样本累积并在 ``[0, 2*pi)`` 内
回卷，正弦值连续，因此不存在相位跳变。调制通过线性插值全通系数实现，
任一时刻的系数仍满足 ``|a| < 1``，故每级幅度仍恒为 1。

反馈把级联输出经 ``fb`` 送回级联输入。由于每级幅度为 1，环路增益恰为
``|fb|``，仅当 ``|fb| < 1`` 时系统稳定；构造函数会拒绝 ``|fb| >= 1``。
"""

import math

__all__ = [
    "Phaser",
    "PhaserError",
    "InvalidRateError",
    "InvalidStagesError",
    "InvalidFeedbackError",
    "InvalidLFOError",
    "InvalidDepthError",
    "InvalidSampleError",
    "allpass_coefficient",
    "allpass_phase",
    "notch_frequency",
]

_TAU = 2.0 * math.pi


# ---------------------------------------------------------------------------
# 可区分的非法参数异常（全部继承 ValueError，便于调用方统一捕获）
# ---------------------------------------------------------------------------
class PhaserError(ValueError):
    """相位器参数 / 输入错误的基类。"""


class InvalidRateError(PhaserError):
    """采样率非法。"""


class InvalidStagesError(PhaserError):
    """级数非法（必须为正整数）。"""


class InvalidFeedbackError(PhaserError):
    """反馈系数非法（必须满足 |fb| < 1）。"""


class InvalidLFOError(PhaserError):
    """LFO 速率非法（必须为有限正数）。"""


class InvalidDepthError(PhaserError):
    """调制深度非法（必须在 [0, 1] 内）。"""


class InvalidSampleError(PhaserError):
    """输入样本非法（NaN / 无穷）。"""


# ---------------------------------------------------------------------------
# 解析辅助函数
# ---------------------------------------------------------------------------
def allpass_coefficient(freq, rate):
    """由中心频率 ``freq`` 与采样率 ``rate`` 给出的一阶全通系数。

    ``0 < freq < rate/2`` 时返回值严格位于 ``(-1, 1)``。
    """
    freq = float(freq)
    rate = float(rate)
    if not math.isfinite(rate) or rate <= 0.0:
        raise InvalidRateError(f"采样率必须为有限正数，收到 {rate!r}")
    if not math.isfinite(freq) or not 0.0 < freq < 0.5 * rate:
        raise ValueError(f"频率必须位于 (0, rate/2) 内，收到 {freq!r}")
    t0 = math.tan(math.pi * freq / rate)
    return (t0 - 1.0) / (t0 + 1.0)


def allpass_phase(freq, rate, coefficient):
    """单级全通在 ``freq`` 处的相位（弧度），在 ``(0, rate/2]`` 上单调下降。"""
    freq = float(freq)
    rate = float(rate)
    coefficient = float(coefficient)
    if not math.isfinite(rate) or rate <= 0.0:
        raise InvalidRateError(f"采样率必须为有限正数，收到 {rate!r}")
    if not -1.0 < coefficient < 1.0:
        raise ValueError(f"全通系数必须满足 |a| < 1，收到 {coefficient!r}")
    if not math.isfinite(freq) or not 0.0 <= freq <= 0.5 * rate:
        raise ValueError(f"频率必须位于 [0, rate/2] 内，收到 {freq!r}")
    t0 = (1.0 + coefficient) / (1.0 - coefficient)
    return -2.0 * math.atan(math.tan(math.pi * freq / rate) / t0)


def notch_frequency(coefficient, stages, rate):
    """首个陷波频率的解析解：``fs/pi * atan(t0 * tan(pi/(2N)))``。"""
    if isinstance(stages, bool) or not isinstance(stages, int):
        raise InvalidStagesError(f"级数必须为正整数，收到 {stages!r}")
    if stages < 1:
        raise InvalidStagesError(f"级数必须为正整数，收到 {stages!r}")
    coefficient = float(coefficient)
    rate = float(rate)
    if not math.isfinite(rate) or rate <= 0.0:
        raise InvalidRateError(f"采样率必须为有限正数，收到 {rate!r}")
    if not -1.0 < coefficient < 1.0:
        raise ValueError(f"全通系数必须满足 |a| < 1，收到 {coefficient!r}")
    t0 = (1.0 + coefficient) / (1.0 - coefficient)
    return float((rate / math.pi) * math.atan(t0 * math.tan(math.pi / (2.0 * stages))))


# ---------------------------------------------------------------------------
# 相位器
# ---------------------------------------------------------------------------
class Phaser:
    """N 级一阶全通级联相位器，带 LFO 调制与反馈。

    参数
    ----
    rate:   采样率（Hz，有限正数）。
    stages: 全通级数（正整数，>= 1）。
    fb:     反馈系数，稳定条件 ``|fb| < 1``。
    lfo_hz: LFO 速率（Hz，有限正数）。
    depth:  调制深度，位于 ``[0, 1]``；0 表示不调制。

    ``process`` 为流式接口：内部状态（反馈记忆、各级状态、LFO 相位）
    在多次调用之间保持连续。
    """

    #: LFO 未调制时的基准中心频率（Hz）。
    CENTER_HZ = 1000.0

    def __init__(self, rate, stages, fb, lfo_hz, depth):
        rate = self._to_float(rate, InvalidRateError, "采样率")
        if not math.isfinite(rate) or rate <= 0.0:
            raise InvalidRateError(f"采样率必须为有限正数，收到 {rate!r}")

        if isinstance(stages, bool) or not isinstance(stages, int):
            raise InvalidStagesError(f"级数必须为正整数，收到 {stages!r}")
        if stages < 1:
            raise InvalidStagesError(f"级数必须为正整数，收到 {stages!r}")

        fb = self._to_float(fb, InvalidFeedbackError, "反馈系数")
        if not math.isfinite(fb) or abs(fb) >= 1.0:
            raise InvalidFeedbackError(f"反馈系数必须满足 |fb| < 1，收到 {fb!r}")

        lfo_hz = self._to_float(lfo_hz, InvalidLFOError, "LFO 速率")
        if not math.isfinite(lfo_hz) or lfo_hz <= 0.0:
            raise InvalidLFOError(f"LFO 速率必须为有限正数，收到 {lfo_hz!r}")

        depth = self._to_float(depth, InvalidDepthError, "调制深度")
        if not math.isfinite(depth) or not 0.0 <= depth <= 1.0:
            raise InvalidDepthError(f"调制深度必须位于 [0, 1] 内，收到 {depth!r}")

        self._rate = rate
        self._stages = stages
        self._fb = fb
        self._lfo_hz = lfo_hz
        self._depth = depth

        # LFO 相位增量（弧度 / 样本）。
        self._inc = _TAU * lfo_hz / rate

        # 调制范围：中心频率按倍频程对称展开，并保证严格低于奈奎斯特频率。
        center = min(self.CENTER_HZ, 0.225 * rate)
        center = max(center, 1e-3)
        span = 2.0 ** depth
        f_lo = min(max(center / span, 1e-6), 0.49 * rate)
        f_hi = min(max(center * span, 1e-6), 0.49 * rate)
        if f_hi <= f_lo:  # 仅 depth == 0 时可能相等
            f_hi = f_lo

        a_lo = allpass_coefficient(f_lo, rate)
        a_hi = allpass_coefficient(f_hi, rate)
        self._a_mid = 0.5 * (a_lo + a_hi)
        self._a_half = 0.5 * (a_hi - a_lo)

        self._reset_state()

    # -- 构造辅助 ---------------------------------------------------------
    @staticmethod
    def _to_float(value, error, name):
        try:
            return float(value)
        except (TypeError, ValueError):
            raise error(f"{name}必须为数值，收到 {value!r}") from None

    def _reset_state(self):
        self._phase = 0.0
        self._fb_state = 0.0
        # 每级全通的状态变量（转置直接 II 型），见 process 注释。
        self._state = [0.0] * self._stages

    # -- 只读属性 ---------------------------------------------------------
    @property
    def rate(self):
        return self._rate

    @property
    def stages(self):
        return self._stages

    @property
    def feedback(self):
        return self._fb

    @property
    def lfo_hz(self):
        return self._lfo_hz

    @property
    def depth(self):
        return self._depth

    @property
    def coefficient(self):
        """LFO 相位为 0 时的全通系数（调制区间的中点）。"""
        return self._a_mid

    def reset(self):
        """清零滤波器与 LFO 状态（参数保持不变）。"""
        self._reset_state()

    # -- 处理 -------------------------------------------------------------
    def process(self, samples):
        """逐样本处理，返回与输入等长的 ``list[float]``。

        非法样本（NaN / 正负无穷）会抛出 :class:`InvalidSampleError`，
        且在预先校验完成前不会改动内部状态。
        """
        if samples is None:
            raise InvalidSampleError("输入样本不得为 None")
        vals = samples if isinstance(samples, list) else list(samples)

        # 先整体校验，避免处理到一半才失败而污染流式状态。
        inf = math.inf
        for v in vals:
            if v != v or v == inf or v == -inf:  # NaN / +inf / -inf
                raise InvalidSampleError(f"输入样本必须为有限值，收到 {v!r}")

        n_out = len(vals)
        if n_out == 0:
            return []

        stages = self._stages
        inc = self._inc
        a_mid = self._a_mid
        a_half = self._a_half
        fb = self._fb
        state = self._state
        phase = self._phase
        fb_state = self._fb_state
        sin = math.sin
        tau = _TAU
        indices = range(stages)

        out = [0.0] * n_out
        for i, v in enumerate(vals):
            # LFO：相位连续累积，sin 连续；系数严格落在 (a_lo, a_hi) 内。
            a = a_mid + a_half * sin(phase)

            # 反馈（环路增益 |fb| < 1）。
            u = v + fb * fb_state

            # N 级全通级联，每级采用转置直接 II 型：
            #   y[n] = a*x[n] + s[n-1]
            #   s[n] = x[n] - a*y[n]
            # 其传递函数为 (a + z^-1)/(1 + a*z^-1)，单位圆上幅度恒为 1。
            for s_i in indices:
                s = state[s_i]
                y = a * u + s
                state[s_i] = u - a * y
                u = y

            fb_state = u
            out[i] = 0.5 * (v + u)

            phase += inc
            if phase >= tau:
                phase %= tau

        self._phase = phase
        self._fb_state = fb_state
        return out

    # -- 陷波 -------------------------------------------------------------
    def notch(self):
        """首个陷波频率（Hz），由级数与 LFO 中点处系数解析给出。"""
        return notch_frequency(self._a_mid, self._stages, self._rate)
