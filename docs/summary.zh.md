# 汇总报告

`gapit summary` 把一份或多份 abricate 格式报告表（包括 gapit screen 的输出）折叠成
基因存在/缺失矩阵：每文件一行，每基因一列。它取代
`abricate --summary`。

输入报告是带标准 15 列表头（`#FILE  SEQUENCE  ...`）的 TSV 或 CSV，正是 `gapit
screen` 写出的格式。如何产生报告见 [./screen.md](./screen.md)；输出 schema 见
[./outputs.md](./outputs.md)。

## 选项

转写自 `gapit summary --help`（gapit 0.5.2）：

| 参数 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `FILE...` | 路径 | stdin | 要汇总的 abricate 格式报告文件。至少一个输入：文件、管道送入的表或 `-` 表示 stdin。 |
| `--identity` / `-i` | 开关 | 关闭 | 单元格显示 %IDENTITY 数值（默认 `+`/`-`）。 |
| `--coverage` / `-c` | 开关 | 关闭 | 单元格显示 %COVERAGE 数值（默认 `+`/`-`）。 |
| `--nopath` | 开关 | 关闭 | 行键只保留文件名（FILE 值 / 输入文件名）。 |
| `--quiet` | 开关 | 关闭 | 静默 stderr 诊断。 |
| `--format` | tsv\|csv\|json\|md | `tsv` | 输出格式。 |

## dutch 模式与多文件模式

行键取决于输入文件的数量，与 abricate 的行为一致：

- **恰好一份报告（dutch 模式）：** 行按该报告的 `FILE` 列取键，报告内每个 assembly
  一行。
- **多于一份报告：** 行按输入文件名的原样取键，每文件一行。零命中的文件仍会出现，
  带 `NUM_FOUND 0`。

```console
$ gapit summary tests/data/summary/multi_sample.tsv
#FILE	NUM_FOUND	feature_a	feature_b
aa_assembly.fa	1	99.00	.
mm_assembly.fa	1	.	50.00
zz_assembly.fa	1	91.00	.
```

同样的三个 assembly 若拆成三份独立报告再汇总，则按文件名取键：

```console
$ gapit summary tests/data/summary/sample_a.tsv tests/data/summary/sample_b.tsv tests/data/summary/empty.tsv
#FILE	NUM_FOUND	feature_a	feature_b
tests/data/summary/empty.tsv	0	.	.
tests/data/summary/sample_a.tsv	2	99.50;52.00	76.00
tests/data/summary/sample_b.tsv	2	90.00	100.00
```

注意那个零命中文件：它以 `NUM_FOUND 0` 和全为 `.` 的基因列出现。同一文件里同一基因
的多个命中按报告顺序用 `;` 串接（见 `99.50;52.00`）。

## 经管道送入 screen 输出

只要 stdin 不是终端（管道或重定向），`gapit summary` 就从 stdin 读取表，因此筛查
阶段可以直接喂给它——规范的批量管道。一次 screen 调用筛查多个 assembly 会写出一张
批量表（每行携带同一个 `DATABASE`；`FILE` 列横跨整个批次），而管道送入的表算作**一个
输入**，于是 dutch 模式生效——每个 FILE 值一行，得到基因×文件存在矩阵：

```console
$ gapit screen -d ecoli_dec *.fna --nopath --quiet | gapit summary
#FILE	NUM_FOUND	aggR	astA	escV	pic	stx2a	uidA
dec_s2_pic_astA_uidA.fasta	3	.	100.00	.	100.00	.	100.00
dec_s3_stx2a_escV_aggR_uidA.fasta	4	100.00	.	100.00	.	100.00	100.00
```

这与先 `gapit screen -d ecoli_dec *.fna --output combined.tsv` 再
`gapit summary combined.tsv` 打印的矩阵完全相同。规则：

- **在终端里裸调用 `gapit summary` 仍打印帮助**（退出码 2）；挂着管道或重定向时改为
  读取 stdin。
- **`gapit summary -`** 是显式 stdin 标记——即使挂着终端也从 stdin 读取，直到 EOF。
- **`-` 与文件参数混用是用法错误**（v1 只接受 stdin 或文件之一，绝不混用）。
- **非 UTF-8 的 stdin** 以指明 `-` 的 `SUMMARY_MALFORMED` 信封退出码 5；畸形的管道行
  会指明 `-` 及其行号。

## 单元格模式

默认单元格为**存在判定**：`+` 表示该文件至少有一条命中，`-` 表示缺失。

| 模式 | 调用方式 | 单元格 |
|---|---|---|
| 存在（默认） | `gapit summary a.tsv b.tsv` | `+` / `-` |
| 一致性 | `--identity` / `-i` | 各命中的 %IDENTITY，`;` 连接（`.` = 缺失） |
| 覆盖度 | `--coverage` / `-c` | 各命中的 %COVERAGE，`;` 连接（`.` = 缺失） |
| 两者 | `-ic`（组合） | 各命中 `identity/coverage`，`;` 连接（`.` = 缺失） |

两旗标可组合（`-ic`）：每个命中渲染为 `identity/coverage`。所选指标在 stderr 提示。
`--coverage` 逐字节复现 abricate `--summary` 的经典单元格形态——parity 对比在该模式下运行。

```console
$ gapit summary tests/data/summary/sample_a.tsv tests/data/summary/sample_b.tsv
#FILE	NUM_FOUND	feature_a	feature_b
sample_a.tsv	2	+	+
sample_b.tsv	2	+	+
```

单元格原样保留原始报告字符串，不做任何重新格式化或平均。`NUM_FOUND` 统计不同基因
的个数。基因全集是所有输入的 `GENE` 值的并集，按字典序排序；某行缺失的基因显示
`.`。

## 输入解析

- **分隔符自动检测，逐文件。** 每个输入通过嗅探自己的首行判断按 TSV 还是 CSV 读
  取，所以一次调用可以混用制表符与逗号报告，无需参数：

```console
$ gapit summary tests/data/summary/sample_a.csv tests/data/summary/sample_b.csv --nopath
#FILE	NUM_FOUND	feature_a	feature_b
sample_a.csv	2	99.50;52.00	76.00
sample_b.csv	2	90.00	100.00
```

- **重复输入会被跳过。** 同一路径出现两次（按原样比较，在任何文件名化之前）会在
  stderr 警告且只处理一次：

```console
$ gapit summary tests/data/summary/sample_a.tsv tests/data/summary/sample_a.tsv
WARNING: Skipping duplicate file: tests/data/summary/sample_a.tsv
#FILE	NUM_FOUND	feature_a	feature_b
tests/data/summary/sample_a.tsv	2	99.50;52.00	76.00
```

- **畸形输入是类型化错误，不是沉默。** 文件缺失以 `INPUT_NOT_FOUND` 信封退出码 5；
  行比表头映射短以 `SUMMARY_MALFORMED` 退出码 5。
- **表头处理与 abricate 的解析器一致。** 第一个文件的首行成为列名映射。首列以 `#`
  开头的行随后按表头跳过，因此常规报告只贡献数据行。这条规则带来两个怪癖：无表头
  报告的首行兼任表头映射和数据行；键以 `#` 开头的文件会被整体跳过（该检查先于
  `--nopath` 文件名化运行）。

## 输出格式

`--format tsv`（默认）与 `--format csv` 打印上面的矩阵。`--format json` 和
`--format md` 输出 `gapit.summary/1` 文档，保留每个单元格字符串并附加机器可读的
params：

```console
$ gapit summary tests/data/summary/sample_a.tsv tests/data/summary/sample_b.tsv --format json --quiet
{
  "schema": "gapit.summary/1",
  "tool": {
    "name": "gapit",
    "version": "0.5.2"
  },
  "created_at": "2026-09-19T01:12:48Z",
  "params": {
    "metric": "%COVERAGE",
    "nopath": false
  },
  "genes": [
    "feature_a",
    "feature_b"
  ],
  "rows": [
    {
      "file": "tests/data/summary/sample_a.tsv",
      "num_found": 2,
      "cells": {
        "feature_a": [
          "99.50",
          "52.00"
        ],
        "feature_b": [
          "76.00"
        ]
      }
    },
    {
      "file": "tests/data/summary/sample_b.tsv",
      "num_found": 2,
      "cells": {
        "feature_a": [
          "90.00"
        ],
        "feature_b": [
          "100.00"
        ]
      }
    }
  ]
}
```

用 `gapit schema summary` 自省 schema；字段契约见
[./outputs.md](./outputs.md)。

## 注意事项

- **按键排序，不是按标签。** 行按键的原样排序，先于 `--nopath` 文件名化。文件散落
  在不同目录时，文件名化后的标签因此可能显得未排序：

```console
$ gapit summary /tmp/gapit-demo/sortdemo/1dir/zeta.tsv /tmp/gapit-demo/sortdemo/2dir/mid.tsv --nopath
#FILE	NUM_FOUND	feature_a	feature_b
zeta.tsv	2	91.00;99.00	50.00
mid.tsv	2	99.50;52.00	76.00
```

  顺序来自 `1dir/zeta.tsv` < `2dir/mid.tsv`，所以 `zeta.tsv` 先打印，尽管 `mid.tsv`
  按字母序更靠前。不带 `--nopath` 时完整键按显而易见的顺序打印：

```console
$ gapit summary /tmp/gapit-demo/sortdemo/1dir/zeta.tsv /tmp/gapit-demo/sortdemo/2dir/mid.tsv
#FILE	NUM_FOUND	feature_a	feature_b
/tmp/gapit-demo/sortdemo/1dir/zeta.tsv	2	91.00;99.00	50.00
/tmp/gapit-demo/sortdemo/2dir/mid.tsv	2	99.50;52.00	76.00
```

- **dutch 模式的键来自报告内部。** 只有一份输入时，行标签是报告内部找到的 `FILE`
  值，不是报告自身的文件名。
- **`--quiet` 抑制 `--identity` 提示与重复文件警告**，从不影响 stdout。
