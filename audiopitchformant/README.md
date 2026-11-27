# audiopitchformant

## 题面（可直接粘贴）

变调不变时长（相位声码器 + 共振峰保持）。约束：变调后基频与目标比误差 ≤1%；输出长度与输入一致（≤1e-12）；共振峰包络（谱包络）在变调后与原始差 ≤5%；COLA 窗和误差 ≤1e-12；100 万样本不超过 900 ms；非法半音与非法帧长可区分报错。 工作目录为仓库内的 audiopitchformant/。 验收：cd audiopitchformant && python -m pytest tests -q。 实现只允许使用 Python 标准库。

## 运行方式

```bash
cd audiopitchformant && python -m pytest tests -q
```

初始 15 failed（实现前）。仅使用 Python 标准库。

## 字段

- 所属领域：游戏、图形与多媒体
- 子方向：音视频处理
- 任务类型：Bugfix
- 语言/框架：Python 3, pytest
- 环境可复现等级：无外部依赖
