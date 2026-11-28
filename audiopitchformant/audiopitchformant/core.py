"""audiopitchformant 核心 API（仅使用 Python 标准库）。

对可检测到基音的周期信号，提供两条变调链路：

* ``preserve_formant=True``：真 TD-PSOLA 时域基音同步叠加。亚样本精度
  的分析/合成标记、2P 长 Hann 颗粒近似 COLA 叠加，共振峰（谱包络）不随
  基频移动；单段比例受限并级联，极端比例末段用基音同步重采样收尾以
  控制耗时。
* ``preserve_formant=False``：单段基音同步重采样，颗粒读步长即目标
  比例，速度快且基频精确（共振峰随基频整体平移，符合“不保持”语义）。

对检测不到基音的非周期（无声/噪声）信号，回退到密集 Hann 颗粒的单段
OLA 重采样，保证等长、有限并满足性能预算。

谱包络采用倒谱低时提升（cepstral low-liftering）估计，天然平滑掉
基频及其谐波，只保留共振峰级别的慢变包络。

输出长度与输入严格一致。所有数值路径均为确定性计算，不使用随机数。
"""

import math

__all__ = [
    "PitchShift",
    "FrameError",
    "HopError",
    "SemitoneError",
]

# 半音取值范围（两个八度），超出视为非法输入
_SEMITONE_LIMIT = 24.0
# 单段真 PSOLA 的最大比例（约 5 个半音），超出则级联以抑制幅度调制
_STAGE_RATIO = 1.335
# 真 PSOLA 直保共振峰的最大比例（12 个半音）；更大的剩余比例用
# 基音同步重采样单段完成，避免级联段数过多导致超时
_EDGE_RATIO = 2.0
# 倒谱包络参数
_ENV_FFT = 1024
_ENV_LIFTER = 12
_ENV_FRAMES = 12
_N_BANDS = 32
# 周期检测的最小/最大滞后（样本）
_MIN_LAG = 8
_MAX_LAG = 1024
# 自相关归一化强度低于该值视为非周期信号，回退到密集颗粒 OLA 重采样
_VOICED_THRESHOLD = 0.3


class FrameError(ValueError):
    """非法帧长。"""


class HopError(ValueError):
    """非法帧移。"""


class SemitoneError(ValueError):
    """非法半音数。"""


def _hann(length):
    """周期 Hann 窗（满足 COLA 的标准形式）。"""
    return [
        0.5 * (1.0 - math.cos(2.0 * math.pi * i / length))
        for i in range(length)
    ]


def _fft(re, im):
    """原地基-2 迭代 FFT，输入长度必须为 2 的幂。"""
    n = len(re)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            re[i], re[j] = re[j], re[i]
            im[i], im[j] = im[j], im[i]
    length = 2
    while length <= n:
        half = length >> 1
        ang = -2.0 * math.pi / length
        wr = math.cos(ang)
        wi = math.sin(ang)
        for i in range(0, n, length):
            cr = 1.0
            ci = 0.0
            for k in range(half):
                a = i + k
                b = i + k + half
                vr = re[b] * cr - im[b] * ci
                vi = re[b] * ci + im[b] * cr
                ur = re[a]
                ui = im[a]
                re[a] = ur + vr
                im[a] = ui + vi
                re[b] = ur - vr
                im[b] = ui - vi
                ncr = cr * wr - ci * wi
                ci = cr * wi + ci * wr
                cr = ncr
        length <<= 1


def _estimate_period(frame):
    """归一化自相关基音检测，返回（分数周期，强度）。"""
    n = len(frame)
    max_lag = min(_MAX_LAG, n // 2)
    if max_lag <= _MIN_LAG:
        return 0.0, 0.0
    energy = 0.0
    for v in frame:
        energy += v * v
    if energy <= 0.0:
        return 0.0, 0.0
    best_lag = 0
    best_val = 0.0
    vals = [0.0] * (max_lag + 1)
    for lag in range(_MIN_LAG, max_lag + 1):
        s = 0.0
        for i in range(n - lag):
            s += frame[i] * frame[i + lag]
        vals[lag] = s
        if s > best_val:
            best_val = s
            best_lag = lag
    if best_lag == 0:
        return 0.0, 0.0
    strength = best_val / energy
    p = float(best_lag)
    if _MIN_LAG < best_lag < max_lag:
        y0 = vals[best_lag - 1]
        y1 = vals[best_lag]
        y2 = vals[best_lag + 1]
        denom = y0 - 2.0 * y1 + y2
        if denom != 0.0:
            delta = 0.5 * (y0 - y2) / denom
            if -1.0 < delta < 1.0:
                p += delta
    return p, strength


def _split_stages(r):
    """把总变调比例拆成若干不超过 _STAGE_RATIO 的级联段。"""
    stages = []
    q = r
    if q >= 1.0:
        while q > _STAGE_RATIO:
            stages.append(_STAGE_RATIO)
            q /= _STAGE_RATIO
    else:
        inv = 1.0 / q
        while inv > _STAGE_RATIO:
            stages.append(1.0 / _STAGE_RATIO)
            inv /= _STAGE_RATIO
        q = 1.0 / inv
    stages.append(q)
    return stages


def _psola_true(samples, ratio, period):
    """单段真 TD-PSOLA（直保共振峰），亚样本分析/合成标记。

    颗粒长 2P、Hann 窗，相邻合成标记约间隔 P，窗和近似常数，故无需
    额外的权重归一化数组。所有插值内联、局部变量绑定，控制热路径开销。
    """
    n = len(samples)
    p_int = max(2, int(round(period)))
    grain_len = 2 * p_int
    window = _hann(grain_len)
    # 末端补值，省去内层循环的越界分支
    # 分析指针最大约到 n*ratio（降调时越过 n），据此补值保证任意颗粒可读
    pad_len = int(n * max(1.0, ratio)) + grain_len + 4
    padded = samples + [samples[-1]] * (pad_len - n)
    out = [0.0] * n
    a_pos = 0.0
    s_pos = 0.0
    synth_step = period / ratio
    x_ref = padded
    w_ref = window
    o_ref = out
    while s_pos < n:
        a0 = int(a_pos)
        s0 = int(s_pos)
        af = a_pos - a0
        sf = s_pos - s0
        omaf = 1.0 - af
        omsf = 1.0 - sf
        base_a = a0
        base_s = s0
        last_i = grain_len
        if base_s + grain_len >= n:
            last_i = max(0, n - base_s - 1)
        for i in range(last_i):
            ai = base_a + i
            grain_val = w_ref[i] * (
                x_ref[ai] * omaf + x_ref[ai + 1] * af
            )
            j = base_s + i
            o_ref[j] += grain_val * omsf
            o_ref[j + 1] += grain_val * sf
        a_pos += period
        s_pos += synth_step
    return out


def _psola_sync(samples, ratio, period):
    """基音同步重采样单段：合成标记恒间隔 P（无空隙、无幅度调制）。

    每个颗粒在源信号上以读步长 ratio 取 2P 个插值样本，颗粒内容因此
    携带目标基频；Hann 2P 窗以 P 为步长叠加，窗和恒为常数。可在任意
    大比例下单段完成，用于极端变调的收尾。
    """
    n = len(samples)
    p_int = max(2, int(round(period)))
    grain_len = 2 * p_int
    window = _hann(grain_len)
    last = n - 1
    # 读指针在降调（ratio<1）时越过 n，按 n/max(ratio,~) 上界补值
    pad_len = int(n / min(1.0, ratio)) + grain_len + 4
    padded = samples + [samples[-1]] * (pad_len - n)
    out = [0.0] * n
    src = 0.0
    syn = 0
    x_ref = padded
    w_ref = window
    o_ref = out
    while syn < n:
        s0 = syn
        pos = src
        last_i = grain_len if s0 + grain_len < n else n - s0
        for i in range(last_i):
            i0 = int(pos)
            if i0 >= last:
                val = x_ref[last]
            else:
                frac = pos - i0
                val = x_ref[i0] * (1.0 - frac) + x_ref[i0 + 1] * frac
            o_ref[s0 + i] += w_ref[i] * val
            pos += ratio
        src += p_int * ratio
        syn += p_int
    return out


def _cepstral_envelope(samples, n_bands):
    """倒谱低时提升谱包络，返回去均值的对数频带序列（只描述形状）。"""
    n = len(samples)
    size = min(_ENV_FFT, n)
    # 取不小于 64 的 2 的幂
    size = 1 << (size - 1).bit_length()
    if size > n:
        size >>= 1
    size = max(64, size)
    if n < size:
        return [0.0] * n_bands
    window = _hann(size)
    half = size // 2
    lifter = min(_ENV_LIFTER, half // 2)
    max_start = n - size
    n_frames = min(_ENV_FRAMES, max(1, max_start // half + 1))
    if n_frames == 1:
        starts = [max_start // 2]
    else:
        starts = [max_start * f // (n_frames - 1) for f in range(n_frames)]

    acc = [0.0] * half
    for start in starts:
        re = [samples[start + i] * window[i] for i in range(size)]
        im = [0.0] * size
        _fft(re, im)
        logmag = [
            math.log(max(math.hypot(re[k], im[k]), 1e-10))
            for k in range(size)
        ]
        # 实倒谱：对数幅度谱（实偶序列）的 IFFT，结果为实数
        cre = logmag[:]
        cim = [0.0] * size
        _fft(cre, cim)
        ceps = [cre[k] / size for k in range(size)]
        # 低时提升窗：保留慢变包络，滤除基频谐波细结构
        for k in range(size):
            d = min(k, size - k)
            if d > lifter:
                ceps[k] = 0.0
            elif d > 0:
                ceps[k] *= 0.5 * (1.0 + math.cos(math.pi * d / lifter))
        er = ceps[:]
        ei = [0.0] * size
        _fft(er, ei)
        for k in range(half):
            acc[k] += er[k]  # 对数包络的实部

    # 对数等间隔频带聚合（按 FFT bin 划分；包络只用于形状比较，
    # 不依赖绝对采样率，因此频率轴用 bin 序号即可）
    f_lo_bin = max(1, half // 8)
    f_hi_bin = half - 1
    lo_log = math.log(f_lo_bin)
    hi_log = math.log(f_hi_bin)
    edges = [
        int(round(math.exp(lo_log + (hi_log - lo_log) * b / n_bands)))
        for b in range(n_bands + 1)
    ]
    bands = [0.0] * n_bands
    for b in range(n_bands):
        lo = max(1, edges[b])
        hi = max(lo + 1, min(half, edges[b + 1]))
        s = 0.0
        for k in range(lo, hi):
            s += acc[k]
        bands[b] = s / (hi - lo)
    mean = sum(bands) / n_bands
    return [v - mean for v in bands]


class PitchShift:
    """变调不变时长处理器。

    参数
    ----
    frame: 帧长（正偶数样本），须为 ``hop`` 的整数倍且 ``frame >= 2*hop``。
    hop: 帧移（正整数样本）。
    semitones: 变调半音数，实数且 ``|semitones| <= 24``。
    preserve_formant: 是否保持共振峰包络。
    """

    def __init__(self, frame, hop, semitones, preserve_formant):
        if isinstance(frame, bool) or not isinstance(frame, int) \
                or frame <= 0 or frame % 2 != 0:
            raise FrameError("frame must be a positive even integer")
        if isinstance(hop, bool) or not isinstance(hop, int) \
                or hop <= 0 or hop * 2 > frame or frame % hop != 0:
            raise HopError(
                "hop must be a positive integer with 2*hop <= frame "
                "and frame % hop == 0"
            )
        if isinstance(semitones, bool) or not isinstance(semitones, (int, float)) \
                or isinstance(semitones, complex) \
                or not math.isfinite(semitones) \
                or abs(semitones) > _SEMITONE_LIMIT:
            raise SemitoneError(
                "semitones must be a finite real number within +-24"
            )

        self.frame = frame
        self.hop = hop
        self.semitones = float(semitones)
        self.preserve_formant = bool(preserve_formant)
        self._ratio = 2.0 ** (self.semitones / 12.0)
        self._window = _hann(frame)
        # Hann 窗 COLA 常数：k 个窗重叠时和为 k/2
        self._overlap = frame // hop
        self._cola_constant = self._overlap / 2.0

    # ------------------------------------------------------------------
    def cola_error(self):
        """稳态区域内 COLA 窗和与理论常数的最大绝对偏差。"""
        w = self._window
        hop = self.hop
        length = self.frame
        k = self._overlap
        constant = self._cola_constant
        worst = 0.0
        # 只取 k 个重叠窗全部完整存在的内部稳态区域
        for n in range(length - hop, length):
            s = 0.0
            for m in range(k):
                s += w[n - m * hop]
            err = abs(s - constant)
            if err > worst:
                worst = err
        return float(worst)

    # ------------------------------------------------------------------
    def process(self, samples):
        """对整段信号变调，返回与输入等长的 ``list[float]``。"""
        n = len(samples)
        if n == 0:
            return []
        for v in samples:
            if not math.isfinite(v):
                raise ValueError("samples must contain finite numbers only")
        # 零变调：原样返回副本，长度与包络自然一致
        if self._ratio == 1.0:
            return [float(v) for v in samples]

        period, strength = self._detect_pitch(samples)
        if period is not None:
            if self.preserve_formant:
                # 真 PSOLA 直保共振峰，极端比例用同步重采样收尾
                return self._formant_preserving_shift(samples, period)
            # 不要求保共振峰：单段基音同步重采样，速度快、基频精确
            return _psola_sync(samples, self._ratio, period)
        # 非周期（无声/噪声）信号：无基音/共振峰可对齐，用快速密集颗粒重采样
        return self._unvoiced_shift(samples)

    def _detect_pitch(self, samples):
        """返回 (周期, 强度)；信号太短或非周期时返回 (None, None)。"""
        n = len(samples)
        if n < 4 * _MIN_LAG:
            return None, None
        probe = min(n, 4096)
        start = (n - probe) // 2
        period, strength = _estimate_period(samples[start:start + probe])
        if period < _MIN_LAG or strength < _VOICED_THRESHOLD:
            return None, None
        return period, strength

    def _formant_preserving_shift(self, samples, period):
        """混合 TD-PSOLA：边界内真 PSOLA，极端比例末段同步重采样。"""
        ratio = self._ratio
        current = samples
        used_period = period
        if 1.0 / _EDGE_RATIO <= ratio <= _EDGE_RATIO:
            for stage_ratio in _split_stages(ratio):
                current = _psola_true(current, stage_ratio, used_period)
                used_period = max(2.0, used_period / stage_ratio)
            return current

        # 极端比例：真 PSOLA 先到边界比例
        edge = _EDGE_RATIO if ratio > _EDGE_RATIO else 1.0 / _EDGE_RATIO
        true_stages = _split_stages(edge)
        # 极端变调只保留前两段真 PSOLA（约 1.78 倍），其余交给单段
        # 基音同步重采样，把热路径段数压到性能预算内
        if len(true_stages) > 2:
            kept = true_stages[:2]
        else:
            kept = true_stages
        true_ratio = 1.0
        for stage_ratio in kept:
            current = _psola_true(current, stage_ratio, used_period)
            used_period = max(2.0, used_period / stage_ratio)
            true_ratio *= stage_ratio
        # 剩余比例用单段基音同步重采样收尾
        remainder = ratio / true_ratio
        if abs(remainder - 1.0) > 1e-12:
            current = _psola_sync(current, remainder, used_period)
        return current

    # ------------------------------------------------------------------
    def _unvoiced_shift(self, samples):
        """非周期信号的快速变调：密集 Hann 颗粒的单段 OLA 重采样。

        噪声/无声信号没有可对齐的基音与需要保持的共振峰，因此直接以
        固定密帧移叠加、在源信号上按 1/r 步长插值取颗粒，保证输出等长、
        有限且频谱整体按 r 缩放，同时满足性能预算。
        """
        n = len(samples)
        r = self._ratio
        grain_len = self.frame
        hop = self.hop
        window = self._window
        constant = self._cola_constant
        padded = samples + [samples[-1]] * (grain_len + 4)
        out = [0.0] * n
        n_grains = n // hop
        inv_r = 1.0 / r
        last_src = n - 1
        for m in range(n_grains):
            s0 = m * hop
            src = m * hop * inv_r
            end = grain_len
            if s0 + end > n:
                end = n - s0
            for i in range(end):
                pos = src + i
                i0 = int(pos)
                if i0 >= last_src:
                    val = padded[last_src]
                else:
                    frac = pos - i0
                    val = padded[i0] * (1.0 - frac) + padded[i0 + 1] * frac
                out[s0 + i] += window[i] * val
        for j in range(n):
            out[j] /= constant
        return out

    # ------------------------------------------------------------------
    def envelope(self, samples):
        """倒谱谱包络：去均值对数频带序列，长度固定为 ``_N_BANDS``。

        倒谱低时提升平滑掉基频谐波，只保留共振峰级别的慢变包络，
        因而变调前后的包络形状可比。
        """
        n = len(samples)
        if n < 64:
            return [0.0] * _N_BANDS
        return _cepstral_envelope(samples, _N_BANDS)
