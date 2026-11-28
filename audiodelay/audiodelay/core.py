"""audiodelay 核心 API：分数延迟线 + 反馈（仅使用 Python 标准库）。

特性
----
* 整数延迟走环形缓冲直读路径，逐样本精确（无插值运算）。
* 分数延迟采用四点三次 Hermite（Catmull-Rom）插值，系数在整个
  处理块开始前一次性预算。
* 反馈系数要求 ``|fb| < 1``，否则构造时抛出 :class:`FeedbackError`。
* 干湿信号等功率混合：``dry = cos(mix*pi/2)``、``wet = sin(mix*pi/2)``，
  满足 ``dry**2 + wet**2 == 1``。
* 各类非法参数使用可区分的 :class:`ValueError` 子类报错。
"""

import math
import numbers

__all__ = [
    "Delay",
    "DelayError",
    "FeedbackError",
    "MixError",
    "SignalError",
]


class DelayError(ValueError):
    """非法延迟参数（``max_delay`` 或 ``delay``）。"""


class FeedbackError(ValueError):
    """反馈系数不满足稳定条件 ``|fb| < 1``。"""


class MixError(ValueError):
    """非法干湿混合系数（``mix`` 必须落在 ``[0, 1]``）。"""


class SignalError(ValueError):
    """输入信号含有 NaN 等非法样本。"""


def _is_real(value):
    return isinstance(value, numbers.Real) and not isinstance(value, bool)


class Delay:
    """带反馈的（分数）延迟线。

    参数
    ----
    max_delay:
        允许的最大延迟（样本数），正整数。
    fb:
        反馈系数，必须为有限实数且 ``|fb| < 1``。
    """

    def __init__(self, max_delay, fb):
        if isinstance(max_delay, bool) or not isinstance(max_delay, int) \
                or max_delay < 1:
            raise DelayError("max_delay 必须为正整数，得到 %r" % (max_delay,))
        if not _is_real(fb) or not math.isfinite(fb) or abs(fb) >= 1.0:
            raise FeedbackError(
                "反馈系数必须为有限实数且满足 |fb| < 1，得到 %r" % (fb,))

        self.max_delay = max_delay
        self.fb = float(fb)

        # 分数插值最多回看 k+2 个样本，额外留 3 个槽位余量；
        # 取 2 的幂以便用按位与完成环形索引。
        need = max_delay + 4
        size = 8
        while size < need:
            size <<= 1
        self._size = size
        self._mask = size - 1
        self._buf = [0.0] * size
        self._w = 0

    # -- 参数校验 --------------------------------------------------------

    def _check_delay(self, delay):
        if not _is_real(delay) or not math.isfinite(delay):
            raise DelayError("delay 必须为有限实数，得到 %r" % (delay,))
        if delay < 0.0 or delay > self.max_delay:
            raise DelayError(
                "delay 须满足 0 <= delay <= max_delay(%d)，得到 %r"
                % (self.max_delay, delay))
        return float(delay)

    def _check_mix(self, mix):
        if not _is_real(mix) or not math.isfinite(mix):
            raise MixError("mix 必须为有限实数，得到 %r" % (mix,))
        if mix < 0.0 or mix > 1.0:
            raise MixError("mix 须在 [0, 1] 区间内，得到 %r" % (mix,))
        return float(mix)

    # -- 处理 ------------------------------------------------------------

    def process(self, samples, delay, mix):
        """对一个样本块做延迟处理，返回等长的 ``list[float]``。

        延迟线保持流式状态：多次调用之间缓冲与写指针连续。
        """
        d = self._check_delay(delay)
        m = self._check_mix(mix)
        if not isinstance(samples, list):
            samples = list(samples)
        if not samples:
            return []

        # 等功率干湿混合增益；端点直接取精确的 0/1，保证纯干/纯湿通路
        # 不引入多余的浮点乘法（整数延迟因而逐样本精确）。
        if m == 0.0:
            wet, dry = 0.0, 1.0
        elif m == 1.0:
            wet, dry = 1.0, 0.0
        else:
            wet = math.sin(m * (math.pi / 2.0))
            dry = math.cos(m * (math.pi / 2.0))

        buf = self._buf
        mask = self._mask
        w = self._w
        fb = self.fb

        k = int(math.floor(d))
        f = d - k

        out = []
        append = out.append

        if f == 0.0:
            # 整数延迟：直读历史样本，延迟通路不含任何插值算术。
            if k == 0:
                # 零延迟：反馈环方程 z = x + fb*z，即 z = x/(1-fb)。
                # 先写入当前样本再读出；fb=0 时纯湿输出与输入逐位相同。
                if fb == 0.0:
                    for x in samples:
                        if x != x:  # NaN 判定（NaN 是唯一不等于自身的值）
                            raise SignalError("输入样本不能为 NaN")
                        xf = float(x)
                        buf[w] = xf
                        append(dry * xf + wet * xf)
                        w = (w + 1) & mask
                else:
                    scale = 1.0 / (1.0 - fb)
                    for x in samples:
                        if x != x:
                            raise SignalError("输入样本不能为 NaN")
                        xf = float(x)
                        z = xf * scale
                        buf[w] = z
                        append(dry * xf + wet * z)
                        w = (w + 1) & mask
            else:
                for x in samples:
                    if x != x:  # NaN 判定（NaN 是唯一不等于自身的值）
                        raise SignalError("输入样本不能为 NaN")
                    xf = float(x)
                    dv = buf[(w - k) & mask]
                    buf[w] = xf + fb * dv
                    append(dry * xf + wet * dv)
                    w = (w + 1) & mask
        else:
            # 分数延迟：三次 Hermite 插值。
            # u 为 p1(时间 n-k-1) 与 p2(时间 n-k) 之间的位置参数。
            u = 1.0 - f
            u2 = u * u
            u3 = u2 * u
            h00 = 2.0 * u3 - 3.0 * u2 + 1.0      # p1 基函数
            h10 = u3 - 2.0 * u2 + u             # m1 基函数
            h01 = -2.0 * u3 + 3.0 * u2          # p2 基函数
            h11 = u3 - u2                       # m2 基函数
            # 切线 m1=(p2-p0)/2, m2=(p3-p1)/2，合并为四点线性系数，
            # 每样本只需 4 次乘法。
            c0 = -0.5 * h10
            c1 = h00 - 0.5 * h11
            c2 = h01 + 0.5 * h10
            c3 = 0.5 * h11

            if k == 0:
                # 零延迟分数情形（0 < d < 1）：p2 对应槽位在写入前是
                # 上一块残留，对当前插值无效，按零样本处理；p3 对应
                # “未来”样本，同样零填充。先算插值 d[n]，再写入
                # z[n] = x[n] + fb*d[n]，输出用真实 z[n] 重算。
                for x in samples:
                    if x != x:
                        raise SignalError("输入样本不能为 NaN")
                    xf = float(x)
                    p0 = buf[(w - 2) & mask]
                    p1 = buf[(w - 1) & mask]
                    dv_pre = c0 * p0 + c1 * p1
                    z = xf + fb * dv_pre
                    buf[w] = z
                    dv = dv_pre + c2 * z
                    append(dry * xf + wet * dv)
                    w = (w + 1) & mask
            else:
                for x in samples:
                    if x != x:
                        raise SignalError("输入样本不能为 NaN")
                    xf = float(x)
                    base = w - k
                    p0 = buf[(base - 2) & mask]
                    p1 = buf[(base - 1) & mask]
                    p2 = buf[base & mask]
                    p3 = buf[(base + 1) & mask]
                    dv = c0 * p0 + c1 * p1 + c2 * p2 + c3 * p3
                    buf[w] = xf + fb * dv
                    append(dry * xf + wet * dv)
                    w = (w + 1) & mask

        self._w = w
        return out

    def impulse(self, delay):
        """返回延迟网络的零状态脉冲响应（长度 ``max_delay + 1``）。

        使用同参数的临时延迟线计算，不改变本对象的流式状态。
        """
        d = self._check_delay(delay)
        probe = Delay(self.max_delay, self.fb)
        block = [1.0] + [0.0] * self.max_delay
        return probe.process(block, d, 1.0)
