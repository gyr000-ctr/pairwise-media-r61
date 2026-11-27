# audiogate

## 题面（可直接粘贴）

噪声门（RMS 包络 + 迟滞）。约束：包络用 RMS 加双指数平滑，开门时刻与解析包络 ≤1e-6；开关门阈必须构成迟滞（开门阈高于关门阈），双向判定与状态机一致（≤1e-12）；保持时间用整数样本计数（不得用浮点累加，避免漂移）；增益斜坡无跳变（相邻差有上界）；100 万样本不超过 500 ms；非法阈值顺序与非法保持时间可区分报错。 工作目录为仓库内的 audiogate/。 验收：cd audiogate && python -m pytest tests -q。 实现只允许使用 Python 标准库。

## 运行方式

```bash
cd audiogate && python -m pytest tests -q
```

初始 15 failed（实现前）。仅使用 Python 标准库。

## 字段

- 所属领域：游戏、图形与多媒体
- 子方向：音视频处理
- 任务类型：Bugfix
- 语言/框架：Python 3, pytest
- 环境可复现等级：无外部依赖
