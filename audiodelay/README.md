# audiodelay

## 题面（可直接粘贴）

延迟线（分数延迟 + 反馈）。约束：整数延迟逐样本精确（≤1e-12）；分数延迟用插值与精确延迟对照 ≤1e-3；反馈系数 |fb|<1 保证稳定（否则报错）；干湿等功率混合；100 万样本不超过 500 ms；非法延迟可区分报错。 工作目录为仓库内的 audiodelay/。 验收：cd audiodelay && python -m pytest tests -q。 实现只允许使用 Python 标准库。

## 运行方式

```bash
cd audiodelay && python -m pytest tests -q
```

初始 15 failed（实现前）。仅使用 Python 标准库。

## 字段

- 所属领域：游戏、图形与多媒体
- 子方向：音视频处理
- 任务类型：Bugfix
- 语言/框架：Python 3, pytest
- 环境可复现等级：无外部依赖
