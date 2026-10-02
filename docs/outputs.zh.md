# 输出

每条 gapit 命令遵循同一个契约：数据走 stdout，诊断走 stderr，错误是带文档化退出码
的类型化 JSON 信封。本页是各格式的参考。schema 带版本号、可被机器发现，agent 无需
阅读本页即可消费。产生这些输出的命令另见 [screen.md](./screen.md)、
[reads.md](./reads.md) 和 [summary.md](./summary.md)。

## stdout 与 stderr

- **stdout 只承载数据**：TSV/CSV 表、JSON 文档、Markdown 报告，或 `gapit db fetch`
  与 `gapit db install` 打印的单行回执。使用 `gapit screen --output PATH` 时数据改写
  入该文件，stdout 保持为空。
- **stderr 承载诊断**：筛查期间的 `Processing: <file>` 等进度行、抓取进度、警告。
  `--quiet` 静默 stderr；绝不触碰 stdout。
- 输出是确定性的：稳定的排序、固定的工具参数、数据载荷内没有墙钟时间戳（只有
  `created_at` 元数据字段）。
- **发射时机**（`gapit screen`）：tsv/csv/md 流式输出 —— 先表头或静态 frontmatter，
  随后每个文件完成的那一刻输出该文件的行/小节（输入顺序；`--jobs` 下按队头阻塞）。
  json 是单一文档，在结束时一次性写出。两种方式的拼接字节完全一致。

## TSV（默认格式）

`gapit screen` 打印 abricate 兼容的 TSV：一行表头，随后每个幸存命中一行。来自三基
因夹具数据库的真实输出：

```console
$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr /tmp/opencode/gapit-own-db/gap.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
/tmp/opencode/gapit-own-db/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
```

### 15 列

| # | 列 | 类型 | 含义 |
|---|---|---|---|
| 1 | `FILE` | string | 按命令行原样的输入路径（`--nopath` 时仅文件名） |
| 2 | `SEQUENCE` | string | 输入文件内的 contig 标识 |
| 3 | `START` | integer | contig 上比对起点，1 基 |
| 4 | `END` | integer | contig 上比对终点，始终大于 `START` |
| 5 | `STRAND` | string | `+` 或 `-`；是受试链上基因的链方向，不是 read 方向 |
| 6 | `GENE` | string | 数据库表头中的基因名 |
| 7 | `COVERAGE` | string | `sstart-send/slen`：基因全长上的比对跨度 |
| 8 | `COVERAGE_MAP` | string | 15 字符比对草图，见下文 |
| 9 | `GAPS` | string | `gapopen/gaps`：缺口开口数 / 缺口列总数 |
| 10 | `%COVERAGE` | float | `100 * (length - gaps) / slen`，显示两位小数 |
| 11 | `%IDENTITY` | float | BLAST 打印的 `pident`，从不过滤，两位小数 |
| 12 | `DATABASE` | string | 来源数据库名 |
| 13 | `ACCESSION` | string | 数据库表头中的登录号或坐标跨度 |
| 14 | `PRODUCT` | string | 产物描述；表头没有时为 `n/a` |
| 15 | `RESISTANCE` | string | 数据库表头中的耐药或功能类别 |

每个输入文件内，行按 `SEQUENCE`（字典序）再按 `START`（数值）排序；文件按参数顺序
出现。

### 两个百分比

两者都直接来自 BLAST 行，不做平滑：

- **%IDENTITY** 是 BLAST 的 `pident`，按 `%.2f` 显示。从不再过滤：`--minid` 在
  blastn 内部经 `-perc_identity` 强制执行。
- **%COVERAGE** 是 `100 * (length - gaps) / slen`，无缺口比对列数除以基因全长。
  `--mincov` 阈值与未取整浮点比较，然后值按 `%.2f` 显示。79.996% 的命中显示为
  `80.00` 但会被丢弃。这是 abricate 的确切行为，刻意保留。

### COVERAGE_MAP

minimap 用 15 个字符格勾画比对落在基因的什么位置：

- 通常 15 格。比对有缺口开口时，最后一格让给 `/` 标记，剩 14 格。
- 每格覆盖 `slen / cells` 个基因碱基。比对起止之间的格子是 `=`，其余是 `.`。有缺
  口时 `/` 落在中间格位置。
- 坐标截断为整数，所以长基因上整段比对也可能有边缘格停留在 `.`。abricate 有同样
  的怪癖，gapit 复现它而不是去"修"。

上文的示例行是 `sul1` 的带缺口比对（94 bp 基因，1 个缺口开口跨 3 列）：

```
========/======
|||||||| ||||||
01234567 89...13   cells 0..13 (14, because gapopen > 0)
                   '/' sits after cell 7, the middle position
```

14 格全是 `=`，因为比对横跨整个基因；`/` 记录这次比对带缺口。无缺口的全长比对打
印 15 个 `=`，如下面这行 `tetA`：

```text
/tmp/opencode/gapit-own-db/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
```

### 变体

- `--format csv` 把分隔符换成 `,`。
- `--noheader` 去掉 `#FILE ...` 表头行。
- `--nopath` 把 `FILE` 列文件名化。
- 不同 query 跨度上的重叠基因都报告。共享同一 `(contig, start, end)` 的命中去重，
  第一条 BLAST 行胜出。默认路径不合并区间；这与 abricate 一致。

### typed 基因数据库：判定是第二条命令

携带 `typing.json` 的基因数据库（**typed** 数据库 —— 捆绑的 `ecoli_dec` 即是示例）
的筛查输出与 untyped 完全相同：上面冻结的 15 列 abricate 表，每种格式都逐字节一致。
判定是两阶段流水线的第二阶段 —— [`gapit typing`](./typing.zh.md#两阶段判定工作流)
读取筛查结果表并渲染判定（真实输出；结果表由 `screen --output` 写出）：

```console
$ gapit screen dec_s2_pic_astA_uidA.fasta --db ecoli_dec --output dec.tsv --quiet
$ gapit typing dec.tsv --quiet
FILE	SCHEME	PHENOTYPE	CONFIDENCE	SCORE	RUNNER_UP	NOTES
dec_s2_pic_astA_uidA.fasta	gb4789_6	EAEC	high	1.0000	EHEC (0.0000)	GB 4789.6: any of aggR/pic/astA
dec_s2_pic_astA_uidA.fasta	risk_monitoring	non-DEC	low	0.0000	STEC (0.0000)	GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE); severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up
```

无 aggR 的 pic+astA 谱是头条分歧 —— GB 4789.6 判 EAEC，而风险监测 scheme 回退到
non-DEC。零命中文件不产生筛查行，因此不会出现在判定输出里（完全没有数据行的表是
`TYPING_NO_DATA` 输入错误）。

## JSON

各筛查入口在 `--format json` 下输出带版本号的 JSON 文档。每个文档以 `schema` 字段
开头，声明其契约。字段名是 snake_case，单位显式（`identity_pct`、`coverage_pct`、
`breadth_pct`）。数值命中字段四舍五入到两位小数，与 TSV 显示一致。

### gapit.report/1（contig 筛查）

真实文档，与上面的 TSV 同一次运行：

```json
{
  "schema": "gapit.report/1",
  "tool": {
    "name": "gapit",
    "version": "0.5.0"
  },
  "created_at": "2026-09-19T01:12:20Z",
  "params": {
    "db": "tinyamr",
    "minid": 80.0,
    "mincov": 80.0,
    "threads": 1
  },
  "files": [
    {
      "file": "/tmp/opencode/gapit-own-db/gap.fa",
      "hits": [
        {
          "sequence": "contig1",
          "start": 1,
          "end": 97,
          "strand": "+",
          "gene": "sul1",
          "coverage": "1-94/94",
          "coverage_map": "========/======",
          "gaps": "1/3",
          "coverage_pct": 100.0,
          "identity_pct": 96.91,
          "database": "tinyamr",
          "accession": "U12338.4:1-940",
          "product": "sulfonamide-resistant dihydropteroate synthase Sul1",
          "resistance": "SULFONAMIDE"
        }
      ]
    }
  ]
}
```

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.report/1` |
| `tool` | object | gapit 的 `{name, version}` |
| `created_at` | string | ISO-8601 UTC 时间戳，秒精度 |
| `params` | object | 生效的筛查参数 |
| `params.db` | string | 数据库名 |
| `params.minid` | number | 最小一致性阈值 |
| `params.mincov` | number | 最小覆盖度阈值 |
| `params.threads` | integer | BLAST 线程数 |
| `files` | array | 每个输入文件一项，含零命中文件 |
| `files[].file` | string | 按原样的输入路径 |
| `files[].hits` | array | 幸存命中，无命中时为空列表 |

每个 hit 与 TSV 列一一对应：

| 字段 | 类型 | 含义 |
|---|---|---|
| `sequence` | string | contig 标识 |
| `start`、`end` | integer | contig 比对跨度，1 基 |
| `strand` | string | `+` 或 `-` |
| `gene` | string | 基因名 |
| `coverage` | string | `sstart-send/slen` |
| `coverage_map` | string | 15 字符 minimap |
| `gaps` | string | `gapopen/gaps` |
| `coverage_pct` | number | 覆盖度百分比，两位小数 |
| `identity_pct` | number | 一致性百分比，两位小数 |
| `database` | string | 来源数据库名 |
| `accession` | string | 登录号或坐标跨度 |
| `product` | string | 产物描述 |
| `resistance` | string | 耐药或功能类别（字段名冻结，为 TSV 一致性保留） |

`gapit screen --merge-fragments`（跨 contig 片段合并；见
[screen.md](./screen.md)）产生的 hit 带两个累加的可选字段，其他任何行都不出现：

| 字段 | 类型 | 含义 |
|---|---|---|
| `merged` | boolean | 恒为 `true`；标记由 >= 2 个基因片段拼成的 hit |
| `fragments` | array | 每个贡献片段一项 `{contig, start, end, strand, identity_pct, coverage_pct}`，按 `(contig, start)` 排序 |

合并 hit 本身携带并集 `%COVERAGE`、按比对长度加权的均值 `%IDENTITY`、求和的
`GAPS`，以及锚点片段的 query 坐标；`sequence` 是逗号串接的 contig 列表。这些字段
对 `gapit.report/1` 是累加的（无重命名、无改型、无 schema 版本号变更）。

### gapit.typing_result/1（从筛查结果判定）

`gapit typing RESULT.tsv [RESULT2.tsv ...]` 读取一份或多份 gapit/abricate 筛查结果
表（TSV 或 CSV），从所有行共享的 `DATABASE` 列解析数据库，并对其 `typing.json` 在
每个 FILE 的基因上评估 —— 每个 `(FILE, GENE)` 折叠为其 `(%IDENTITY, %COVERAGE)`
最佳的行。真实文档（有裁剪）：

```json
{
  "schema": "gapit.typing_result/1",
  "tool": {"name": "gapit", "version": "0.5.0"},
  "created_at": "2026-10-02T10:04:55Z",
  "source": ["dec.tsv"],
  "db": "ecoli_dec",
  "files": [
    {
      "file": "dec_s2_pic_astA_uidA.fasta",
      "phenotypes": {
        "gb4789_6": {
          "phenotype": "EAEC",
          "score": 1.0,
          "confidence": "high",
          "components": [
            {"name": "aggR", "score": 0.0},
            {"name": "pic", "score": 1.0},
            {"name": "astA", "score": 1.0}
          ],
          "runner_up": {"phenotype": "EHEC", "score": 0.0},
          "notes": ["GB 4789.6: any of aggR/pic/astA"]
        },
        "risk_monitoring": {
          "phenotype": "non-DEC",
          "score": 0.0,
          "confidence": "low",
          "components": [
            {"name": "escV", "score": 0.0},
            {"name": "stx1a", "score": 0.0},
            {"name": "stx1b", "score": 0.0},
            {"name": "stx2a", "score": 0.0},
            {"name": "stx2b", "score": 0.0}
          ],
          "runner_up": {"phenotype": "STEC", "score": 0.0},
          "notes": [
            "GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE)",
            "severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up"
          ]
        }
      }
    }
  ]
}
```

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.typing_result/1` |
| `tool` | object | gapit 的 `{name, version}` |
| `created_at` | string | ISO-8601 UTC 时间戳，秒精度 |
| `source` | array | 按输入原样的结果表路径 |
| `db` | string | 所有行共同筛查的数据库 |
| `files` | array | 每个有数据行的 FILE 一项，按首次出现顺序 |
| `files[].phenotypes` | object | scheme 名 → 其判定（下述累加结构） |

每条判定携带 `phenotype`（歧义为 null）、`score`、`confidence`
（`high`/`ambiguous`/`low`）、`components[{name, score}]`，以及可选的 `runner_up`、
`ambiguous[]`（并列的一对）与 `notes[]` —— 与评估器在基因簇路径上产出的结构相同。
TSV/Markdown 投影把每条判定压平为 `FILE`、`SCHEME`、`PHENOTYPE`、`CONFIDENCE`、
`SCORE`、`RUNNER_UP`、`NOTES` 七列：歧义判定的表型渲染 `-`、候选对写入 NOTES
（runner-up 单元格同样为 `-`——候选对已说明一切）。类型化错误：行间 `DATABASE`
值混杂为 `DATABASE_MISMATCH`，无数据行的表为 `TYPING_NO_DATA`，无 `typing.json` 的
数据库为 `TYPING_NO_SCHEME`，基因簇数据库为 `TYPING_CLUSTER_DB`（其判定集成在
`gapit screen` 内）。其 JSON Schema 用 `gapit schema typing_result` 打印。

### gapit.reads/1（reads 与 assembly 筛查）

reads 模式报告的是每个基因在整个 read 集（或 assembly；`--r1` 接受 FASTQ 和
FASTA）上的覆盖情况，而不是逐命中行。基因条目按 `breadth_pct` 降序；零 mapped
reads 的基因省略。字段来自 `gapit schema reads`：

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.reads/1` |
| `tool` | object | gapit 的 `{name, version}` |
| `created_at` | string | ISO-8601 UTC 时间戳，秒精度 |
| `params` | object | 生效的 read 筛查参数 |
| `params.db` | string | 数据库名 |
| `params.read_type` | string | 生效的 minimap2 预设：`sr`、`map-ont` 或 `map-hifi`（按检测到的输入解析） |
| `params.min_breadth` | number | 存在阈值，默认 90.0 |
| `params.threads` | integer | minimap2 线程数 |
| `files` | array | 每个 read 集一项 |
| `files[].reads` | array of string | 输入 read/assembly 路径（所有 lane） |
| `files[].genes` | array | 逐基因存在判定 |

每个基因条目：

| 字段 | 类型 | 含义 |
|---|---|---|
| `gene` | string | 基因名 |
| `database` | string | 来源数据库名 |
| `accession` | string | 登录号或坐标跨度 |
| `product` | string | 产物描述 |
| `resistance` | string | 功能类别 |
| `tlen` | integer | 基因长度（碱基） |
| `breadth_pct` | number | 至少被一条主比对覆盖的基因碱基百分比 |
| `mean_depth` | number | 基因全长上的平均逐碱基深度 |
| `reads_mapped` | integer | 该基因上有主比对的去重 reads 数 |
| `present` | boolean | `breadth_pct >= min_breadth` |

### gapit.reads/2（带一致性/MAPQ 过滤的 reads 筛查）

reads/1 的可选变体，仅在 `--min-identity` 或 `--min-mapq` 非零时输出：低于阈值的
PAF 比对在聚合前被丢弃，消除只看广度的存在判定在同源基因上的家族拆分型过量判定。
两个阈值都关闭时输出保持 `gapit.reads/1`，逐字节一致。字段来自
`gapit schema reads2`，除以下各项外与 reads/1 相同：

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.reads/2` |
| `params.min_identity` | number | 生效的单比对一致性下限，0 = 关闭 |
| `params.min_mapq` | integer | 生效的单比对 MAPQ 下限，0 = 关闭 |
| `files[].genes[].mean_identity_pct` | number | 幸存比对上按比对长度加权的单比对一致性均值，保留 2 位 |

单比对一致性是 `100 * (alen - nm) / alen`（PAF block 长度与 `NM:i:` 标签；没有 NM
的行按 100 计）。同源实战示例与阈值建议见
[reads.md](./reads.md)。

### gapit.cluster/1（基因簇数据库筛查）

筛查 `kind: cluster` 数据库（minimap2 `asm20` 引擎；见 [screen.md](./screen.md)）
每次运行输出一个文档：每文件一个最佳位点判定、每个被覆盖位点及其逐基因判定，以及
在携带 `typing.json` 的数据库上的表型判定和可解释的分数分解。字段来自
`gapit schema cluster`：

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.cluster/1` |
| `tool` | object | gapit 的 `{name, version}` |
| `created_at` | string | ISO-8601 UTC 时间戳，秒精度 |
| `params` | object | `db`、`preset`（恒为 `asm20`）、`min_gene_cov`、`min_gene_id`、`min_cluster_cov`、`threads` |
| `files` | array | 每个输入文件一项 |
| `files[].file` | string | 按原样的输入路径 |
| `files[].best` | object \| null | 最佳位点判定，没有位点达到 `min_cluster_cov` 时为 null |
| `files[].loci` | array | 每个被覆盖位点按排名排列（覆盖度降序、一致性降序、覆盖 bp 降序、id 升序） |

`best` 块：

| 字段 | 类型 | 含义 |
|---|---|---|
| `best.locus` | string | 排名第一的位点 id |
| `best.label` / `best.type` | string | 源 note 中的位点标签和类型（kaptive 的 `K locus`/`K type` 语义） |
| `best.coverage_pct` / `best.identity_pct` | number | 全部记录的并集覆盖度与按匹配加权的位点一致性，两位小数 |
| `best.genes_present` / `genes_partial` / `genes_absent` | integer | 最佳位点注释基因的判定计数 |
| `best.phenotype` | string \| null | 表型判定；歧义或数据库无 typing 时为 null |
| `best.phenotype_detail` | object | 累加字段，仅在 typed 运行出现：`{score, confidence, components[{name, score}], runner_up, ambiguous[]}` |

每个 `loci[]` 条目携带 `locus`、`label`、`type`、`coverage_pct`、`identity_pct`、
`rank`、`missing`（未完全存在的基因 id），以及 `genes[]`：逐基因的 `start`/`end`
（1 基，位点坐标）、`strand`、`coverage_pct`、`identity_pct` 和 `verdict`
（`present`/`partial`/`absent`）。

`phenotype_detail.confidence` 是 `high`（胜者 ≥ cutoff 且领先 ≥
`ambiguity_margin`）、`ambiguous`（两条规则落在容差内：表型为 null，前两名列入
`ambiguous`）或 `low`（低于 cutoff：文档的 fallback 字符串）。分量名即被评分规则
自身的名字：`weighted_genes` 用基因 id，`cluster_match` 用 `coverage`/`identity`/
`key_genes`，`learned_linear` 用特征字符串（加 `bias`）。

TSV 形式是每文件一行，typed 数据库的表头是 `FILE BEST_LOCUS TYPE PHENOTYPE COVERAGE
IDENTITY PRESENT PARTIAL MISSING_IDS`（歧义或未判定时 `PHENOTYPE` 为 `-`），无
typing 的数据库是同一个表头去掉 `PHENOTYPE`。没有位点判定的文件渲染成
`- - - 0.00 0.00 0 0 -` 形状的横线。

### gapit.features/1 与 gapit.typing/1（数据库侧文档）

两份基因簇数据库产物是可用 `gapit schema` 自省的带版本文档，但绝不出现在 stdout：
`features`（`db build` 写出的 `features.json` 位点/基因特征表）和 `typing`（用
`db build --typing FILE` 安装的声明式评分规范：`weighted_genes` / `cluster_match` /
`learned_linear` 规则，外加 `cutoff`、`ambiguity_margin` 和 `fallback`）。typing
文档就是表型求值器实现的契约；见 [databases.md](./databases.md)。

### gapit.summary/1（汇总矩阵）

`gapit summary --format json` 把报告表折叠成基因矩阵。真实输出，汇总两份报告：一份
有 `tetA` 命中，一份没有：

```json
{
  "schema": "gapit.summary/1",
  "tool": {
    "name": "gapit",
    "version": "0.5.0"
  },
  "created_at": "2026-09-19T01:11:09Z",
  "params": {
    "metric": "%COVERAGE",
    "nopath": false
  },
  "genes": [
    "tetA"
  ],
  "rows": [
    {
      "file": "/tmp/opencode/gapit-docs-Ifvc44/full.tsv",
      "num_found": 1,
      "cells": {
        "tetA": [
          "100.00"
        ]
      }
    },
    {
      "file": "/tmp/opencode/gapit-docs-Ifvc44/partial.tsv",
      "num_found": 0,
      "cells": {}
    }
  ]
}
```

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.summary/1` |
| `tool` | object | gapit 的 `{name, version}` |
| `created_at` | string | ISO-8601 UTC 时间戳，秒精度 |
| `params.metric` | string | `%COVERAGE` 或 `%IDENTITY`，随 `--identity` 而定 |
| `params.nopath` | boolean | 行标签是否被文件名化 |
| `genes` | array of string | 所有输入基因名的排序并集 |
| `rows` | array | 每个输入文件一行 |
| `rows[].file` | string | 输入报告路径（行标签） |
| `rows[].num_found` | integer | 该报告中不同基因的计数 |
| `rows[].cells` | object | 基因名到原始单元格字符串列表的映射；缺失的基因省略 |

单元格保留输入报告中的精确字符串、按文件顺序，所以 %COVERAGE 和 %IDENTITY 汇总都
能无重格式化地往返。

### 稳定性策略

schema 名带 semver（`gapit.report/1`）。在一个主版本内，既有字段绝不重命名或改
型；累加式变更伴随 schema 版本号变更和 `PLAN.md` 里的说明。你构建的任何东西都可以
对照 `gapit schema` 打印的 schema 校验。

## Markdown

`--format md` 把同样的数据渲染成人类可读的 Markdown 报告：**静态** YAML frontmatter
（`schema`、`tool`、`created_at`、`db`、`minid`、`mincov`、`threads` —— 不含运行总数，
总数由 JSON 文档承载），随后每个输入文件一节，用管道表承载与 TSV 相同的 15 列。
frontmatter 在第一个文件之前即可知，因此 md 与 tsv 行一样逐文件流式输出。reads 模式
和汇总模式有类似的 Markdown 形式。frontmatter 给解析器一个稳定的表头；表格在终端
或编辑器里读起来很自然。

## 错误

任何失败都会向 stderr 打印一行 JSON：

```text
{"schema": "gapit.error/1", "code": "<CODE>", "message": "...", "context": {...}}
```

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.error/1` |
| `code` | string | 机器稳定的错误码，snake_case |
| `message` | string | 人类可读的描述 |
| `context` | object | 字符串值的细节（文件路径、摘要、数据库名） |

退出码告诉你失败类别：

| 退出码 | 含义 |
|---|---|
| 0 | 成功 |
| 1 | 意外错误（信封代码 `UNEXPECTED`） |
| 2 | 用法错误 |
| 3 | 缺少依赖 |
| 4 | 数据库错误 |
| 5 | 输入错误 |

筛查不存在的文件退出码 5，stderr 真实输出：

```console
$ gapit screen /nonexistent/contigs.fa --db ncbi
{"schema":"gapit.error/1","code":"INPUT_NOT_FOUND","message":"input file not found or unreadable: /nonexistent/contigs.fa","context":{"file":"/nonexistent/contigs.fa"}}
$ echo $?
5
```

再举两例，各自由所示命令产生：

```console
$ gapit db fetch nosuchdb            # exit 2, usage
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"unknown database: nosuchdb (available: argannot, bacmet2, card, ecoh, ecoli_vf, megares, ncbi, plasmidfinder, resfinder, upec_expec_vf, vfdb, victors)","context":{"db":"nosuchdb"}}

$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr contigs.fa   # exit 4, db error, before indexing
{"schema":"gapit.error/1","code":"DATABASE_NOT_INDEXED","message":"Database /tmp/opencode/gapit-own-db/db/tinyamr/sequences is not indexed, please try: gapit setupdb","context":{"db":"/tmp/opencode/gapit-own-db/db/tinyamr/sequences"}}
```

缺外部二进制（退出码 3）时会报告是哪一个：

```console
{"schema":"gapit.error/1","code":"MISSING_DEPENDENCY","message":"required binary not found on PATH: blastdbcmd","context":{"binary":"blastdbcmd"}}
```

失败时 stdout 保持为空：没有残缺的表，也没有残缺的 JSON。

## 自描述

agent 只凭二进制就能发现整个契约。

```console
$ gapit --version --json
{"schema":"gapit.version/1","name":"gapit","version":"0.5.0"}
```

`gapit schema <name>` 打印每个文档的 JSON Schema。十个名字：`report`、
`typing_result`、`reads`、`reads2`、`cluster`、`summary`、`error`、`version`，外加
数据库侧文档 `features`（gapit.features/1，基因簇数据库的特征表）和 `typing`
（gapit.typing/1，基因簇数据
库的声明式评分规范）。`gapit schema report` 的裁剪片段：

```json
{
  "description": "gapit.report/1 \u2014 the canonical machine-readable screening output.",
  "properties": {
    "schema": {
      "const": "gapit.report/1",
      "default": "gapit.report/1",
      "title": "Schema",
      "type": "string"
    },
    "created_at": {
      "title": "Created At",
      "type": "string"
    },
    "params": {
      "$ref": "#/$defs/ParamsDocument"
    },
    "files": {
      "items": {
        "$ref": "#/$defs/FileDocument"
      },
      "title": "Files",
      "type": "array"
    }
  },
  "required": [
    "created_at",
    "params",
    "files"
  ],
  "title": "ReportDocument",
  "type": "object"
}
```

消费文档之前，把任何 JSON Schema 校验器指向这个输出即可校验。
[agent 指南](./agents.md)展示了完整的自省工作流。
