# gseek 引擎（可选）

gapit 可以把对齐子进程切到 [gseek](https://github.com/indexofire/gseek)，但原有
子进程实现仍是默认值，并且始终可用。

## 开关

```bash
GAPIT_ENGINE=legacy   # 默认：原有 blastn/blastx/minimap2 调用
GAPIT_ENGINE=gseek    # 由 gseek 规划并解析同样的对齐
GAPIT_ENGINE=shadow   # 两套都跑，记录粗粒度差异，返回 legacy
```

`GAPIT_ENGINE` 每次调用时读取。非法值会在任何对齐启动前报错。`shadow` 是验证
模式，不是第三条生产路径：它始终返回 legacy 结果。

## 影响范围

| 路径 | legacy | gseek |
| --- | --- | --- |
| contig 筛查 BLAST | `blastn`/`blastx` + 表格解析 | `gseek.align(backend="blast")` |
| 蛋白筛查 | `blastx -task blastx-fast` | `gseek.align(backend="blast", task="blastx-fast")` |
| reads 筛查 | 流式 `minimap2` PAF | `gseek.align(backend="minimap2")` |
| cluster 筛查 | `minimap2 -x asm20 --cs` | `gseek.align(backend="minimap2", preset="asm20", cs=True)` |

gseek 模式保留 gapit 的阈值、坐标约定和输出契约，只改变命令规划与行解析。
gseek 模式下 minimap2 PAF 由 gseek API 按 lane 物化，legacy 路径则是流式的。
gseek 返回确定性的排序视图，因此原始 PAF 行序可能与 minimap2 不同；gapit 的
reads/cluster 消费者会聚合/排序，最终报告不变。bridge 还保留两个消费者兼容
处理：VFDB 的重复 subject ID 由 gseek 全量保留；CARD 中超 255 字符的 target
名改用直接 FASTA target（`index="none"`），因为 minimap2 的 `.mmi` 索引会把
这类名称替换成数字序号。

## 安装

gapit 要求 Python >= 3.12，所有支持的环境都包含 gseek：

```bash
pip install "gapit[gseek]"
```

pixi 的 `default` 环境（启用 `py312` feature）使用 Python 3.12，并包含 `gseek>=0.9.1`。
没有 gseek 时，`GAPIT_ENGINE=gseek` 会抛类型化的 `MISSING_DEPENDENCY` 错误；
`legacy` 和所有现有命令不受影响。

## 为什么保留 legacy

切换本身使用 `gseek.compat.Engine`，与其他消费者共享 legacy/gseek/shadow
选择与回退逻辑；gapit 只保留自己的 legacy 实现、gseek 调用和结果语义。

这是有意的双路径设计：

- `legacy` 是所有安装的默认值；
- 每条已迁移路径都保留原实现，仍可正常使用；
- `shadow` 同时跑两套，记录粗粒度差异，并返回 legacy 结果；
- legacy 是 parity 基准，也是 gseek 不可用时的回退；
- 数据库构建、typing、格式、MCP 与 CLI 行为 **不** 经过 gseek。

## 验证

bridge 单测覆盖 `GAPIT_ENGINE`、坐标/链方向映射、蛋白任务选择、PAF/NM/cs
标签、分派和 shadow 回退。真实 parity 测试会构建小型核酸/蛋白 BLAST 数据库
和 minimap2 target，逐行比较 legacy 与 gseek；缺少 gseek 或外部二进制时 skip。
