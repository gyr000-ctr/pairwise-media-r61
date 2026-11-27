# audiosaturate

## 题面（可直接粘贴）

饱和（tanh 近似 + 过采样）。约束：饱和曲线奇对称（≤1e-12）且单调；4 倍过采样后混叠下降 ≥6 dB；驱动参数单调提升谐波；输出有界（≤1）；100 万样本不超过 600 ms；非法驱动可区分报错。 工作目录为仓库内的 audiosaturate/。 验收：cd audiosaturate && python -m pytest tests -q。 实现只允许使用 Python 标准库。

## 运行方式

```bash
cd audiosaturate && python -m pytest tests -q
```

初始 15 failed（实现前）。仅使用 Python 标准库。

## 字段

- 所属领域：游戏、图形与多媒体
- 子方向：音视频处理
- 任务类型：Bugfix
- 语言/框架：Python 3, pytest
- 环境可复现等级：无外部依赖
