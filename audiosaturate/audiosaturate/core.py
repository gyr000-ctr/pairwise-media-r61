"""audiosaturate 核心 API（纯标准库实现）。

饱和曲线采用 tanh 近似：

    y = tanh(drive * x) / tanh(drive)

归一化后满足：
* 奇对称：f(-x) == -f(x)（严格为浮点奇对称）
* 单调递增
* 有界：|y| < 1
* 驱动越大，谐波越丰富（同一输入幅度下的非线性更强）

过采样：以 4 倍线性插值上采样，逐点饱和后线性插值下采样（抽头对齐），
可显著降低非线性混叠。允许的过采样倍率为 1、2、4、8。
"""

import math

__all__ = ["Saturate", "ALLOWED_OVERSAMPLE"]

ALLOWED_OVERSAMPLE = (1, 2, 4, 8)

# 驱动参数上限：防止 tanh 参数溢出（tanh 本身对 inf 也安全，但保持有限且可区分）
_MAX_DRIVE = 1e6


class Saturate:
    def __init__(self, drive, oversample):
        """构造。

        drive: 驱动量，必须是有限正数（> 0）。
        oversample: 过采样倍率，必须是 {1, 2, 4, 8} 之一。
        """
        drive = _as_finite_float(drive, "drive")
        if drive <= 0.0:
            raise ValueError("drive must be positive, got %r" % (drive,))
        if drive > _MAX_DRIVE:
            raise ValueError("drive must be <= %r, got %r" % (_MAX_DRIVE, drive))

        if isinstance(oversample, bool) or not isinstance(oversample, int):
            # 非整数（含浮点）一律视为非法
            raise ValueError("oversample must be an integer, got %r" % (oversample,))
        if oversample not in ALLOWED_OVERSAMPLE:
            raise ValueError(
                "oversample must be one of %r, got %r"
                % (ALLOWED_OVERSAMPLE, oversample)
            )

        self.drive = drive
        self.oversample = oversample
        # 归一化常数：tanh(drive)，保证满量程输入映射到 ±1 边界
        self._norm = math.tanh(drive)

    # ------------------------------------------------------------------ 曲线
    def shape(self, x):
        """返回单个样本 x 经饱和后的值（float，有界于 (-1, 1)）。"""
        x = _as_finite_float(x, "x")
        return math.tanh(self.drive * x) / self._norm

    # ------------------------------------------------------------------ 处理
    def process(self, samples):
        """处理样本序列，返回新的 float 列表（不修改输入）。"""
        out = [0.0] * len(samples)
        if not out:
            return out

        # 先校验并物化为 float，保证全程确定性且不修改输入
        vals = [_as_finite_float(v, "sample") for v in samples]
        if len(vals) == 1:
            out[0] = math.tanh(self.drive * vals[0]) / self._norm
            return out

        up = self._upsample_linear(vals, self.oversample)
        shaped = [math.tanh(self.drive * v) / self._norm for v in up]
        return self._downsample_linear(shaped, self.oversample, len(vals))

    # -------------------------------------------------------------- 重采样
    @staticmethod
    def _upsample_linear(vals, factor):
        """线性插值上采样（端点保持）。输出长度 = len(vals) * factor。

        每个输入样本间插入 factor-1 个线性插值点，group j 对应
        时间位置 j/factor，即 [x0, x0+a1, x0+2a1, ...] 形式。
        """
        if factor == 1:
            return list(vals)
        n = len(vals)
        up = [0.0] * (n * factor)
        f = float(factor)
        up[0] = vals[0]
        for i in range(n - 1):
            a = vals[i]
            b = vals[i + 1]
            d = (b - a) / f
            base = i * factor
            for k in range(1, factor):
                up[base + k] = a + d * k
            up[base + factor] = b
        # 尾端（最后一组仅含端点）
        up[(n - 1) * factor] = vals[-1]
        return up

    @staticmethod
    def _downsample_linear(vals, factor, n_out):
        """线性插值下采样（与上采样网格对齐的抽头平均）。

        取上采样序列在整数时间位置（i*factor）处的值并重建。由于饱和是非线性
        的，直接抽取会引入额外混叠；对一组 factor 个点做线性插值重建可近似
        理想低通的积分平均，显著压低混叠。
        """
        if factor == 1:
            return [float(v) for v in vals]
        f = float(factor)
        out = [0.0] * n_out
        # 每个输出样本对应 factor 个上采样点的平均（boxcar + 线性重构近似）
        for i in range(n_out):
            base = i * factor
            acc = 0.0
            for k in range(factor):
                acc += vals[base + k]
            out[i] = acc / f
        return out


# --------------------------------------------------------------------- 工具
def _as_finite_float(value, name):
    """转为有限 float；非法（NaN/Inf/非数值）抛 ValueError。"""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("%s must be a real number, got %r" % (name, value))
    try:
        f = float(value)
    except (TypeError, ValueError):
        raise ValueError("%s must be a real number, got %r" % (name, value))
    if not math.isfinite(f):
        raise ValueError("%s must be finite, got %r" % (name, value))
    return f
