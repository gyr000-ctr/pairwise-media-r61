# audiophaser

## 题面（可直接粘贴）

相位器（全通级联 + LFO）。约束：每级全通幅度恒 1（≤1e-9）而相位单调下降；陷波频率与级数和系数解析对应（≤1e-6）；LFO 相位连续不跳变；反馈稳定（|fb|<1）；100 万样本不超过 600 ms；非法级数可区分报错。 工作目录为仓库内的 audiophaser/。 验收：cd audiophaser && python -m pytest tests -q。 实现只允许使用 Python 标准库。

## 运行方式

```bash
cd audiophaser && python -m pytest tests -q
```

初始 15 failed（实现前）。仅使用 Python 标准库。

## 字段

- 所属领域：游戏、图形与多媒体
- 子方向：音视频处理
- 任务类型：Bugfix
- 语言/框架：Python 3, pytest
- 环境可复现等级：无外部依赖
