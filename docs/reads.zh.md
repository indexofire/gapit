# reads 与 assembly 筛查（FASTQ / FASTA）

`gapit screen --r1/--r2` 用 minimap2 把原始 FASTQ reads 或整条 assembly FASTA 比对到
基因数据库，按比对广度判定基因存在。abricate 完全不能筛查 reads；这个模式是 gapit
扩展 —— 存在阈值（`--min-breadth`）与 contig 流水线的一致性/覆盖度下限不同，但输出
格式默认值（tsv）与全局的人类默认值一致。

无需 BLAST 索引：minimap2 在内存中为数据库的 `sequences` FASTA 建索引。数据库：
[./databases.md](./databases.md)。contig 模式：[./screen.md](./screen.md)。

## 调用方式

```text
gapit screen --r1 R1[,R1b,...] [--r2 R2[,R2b,...]] --db NAME [--read-type sr|map-ont|map-hifi]
```

- `--r1` 接受 FASTQ reads 或 FASTA assembly（见下文）；逗号分隔，每条 lane 一项。
- 第 i 条 lane 把 `r1[i]` 与 `r2[i]` 配对，所以 `--r2` 的数量必须等于 `--r1`。
- 所有 lane 聚合成一个样本：逐基因指标汇总每条 lane 的比对，一个基因可以靠低于阈值
  的各 lane 的并集越过广度阈值。
- reads 模式与位置参数 contig 文件互斥 —— 唯一的例外是下文的全 FASTQ 通配符，它以
  位置参数方式进入 reads 模式。只有 `--r2` 没有 `--r1` 是用法错误（退出码 2）。支持
  gzip 输入。
- **输入检测。** 每个 `--r1`/`--r2` 文件在验证时按内容检测：首个非空白字节是 `>` 即
  FASTA，`@` 即 FASTQ（gzip 包装的文件透过解压器窥探）。其他情况是输入错误，退出码
  5，代码 `INVALID_READS_FORMAT`。同一个 `--r1` 列表里混用 FASTA 和 FASTQ、FASTA 配
  `--r2`、给 FASTA 显式指定 `sr`/`map-hifi` 预设，都是用法错误（退出码 2）。

## 用通配符筛查一批 FASTQ 文件（自动配对样本）

`--r1`/`--r2` 需要你自己点名各条 lane。批量任务里更省事的做法是把 shell 的通配符
展开直接交给 `screen`，由 gapit 推断样本：

```console
$ gapit screen -d ecoli_dec *.gz
```

当**每个**位置参数文件都是 FASTQ —— 扩展名是 `.fastq`/`.fq`（可带 `.gz`），或扩展名
含糊时内容嗅探为 FASTQ —— 命令即无需 `--r1`/`--r2` 进入 reads 模式。文件按文件名
（大小写不敏感地剥掉 reads 扩展名之后）分组为样本：

| mate 标记（最长后缀优先） | 示例一对 | 样本键 |
|---|---|---|
| `_R1_001` / `_R2_001`（bcl2fastq） | `run_S1_L001_R1_001.fastq.gz` + `run_S1_L001_R2_001.fastq.gz` | `run_S1_L001` —— `_L00x` lane 标记保留在样本键里 |
| `_R1` / `_R2` | `s2_R1.fq.gz` + `s2_R2.fq.gz` | `s2` |
| `_1` / `_2` | `s1_1.fq.gz` + `s1_2.fq.gz` | `s1` |
| `.1` / `.2` | `a.1.fastq` + `a.2.fastq` | `a` |

标记匹配大小写不敏感。没有可识别标记的文件自成一个单端样本，样本键为其剥掉扩展名
的词干。文档**每个样本一个 `files[]` 条目**（条目的 `reads` 字段是样本键；lane 文件
列在 stderr），逐基因指标取该样本各 lane 的并集：

```json
"files": [
  {"reads": ["s1"], "genes": [{"gene": "tetX", "...": "..."}]},
  {"reads": ["s2"], "genes": [{"gene": "tetX", "...": "..."}]}
]
```

- **配不上对的文件警告而不报错。** 标记了 mate 但 glob 里缺失另一半的文件以单端筛
  查，stderr 提示一行 `WARNING: no mate found for X — screening single-end`（用
  `--quiet` 静默）；无标记文件同样处理。
- **完全同名的文件合并**为同一样本的额外 lane（分 lane 目录的多 lane 测序），按排序
  顺序配对；多出来的文件退回单端并警告。
- 样本按字典序排序，输出与 glob 顺序无关、完全确定。
- reads 模式的用法守卫原样生效：拒绝 `--minid`/`--mincov`、`--fofn`、`--noheader`、
  `--nopath`、`--merge-fragments`；而 `--min-identity`/`--min-mapq` 可用（在这里它们
  就是 reads 模式参数，不再是错误）。`--jobs` 只在这条通配符路径上合法
  （`--r1`/`--r2` 单样本调用拒绝它）——见下文
  [并行与流式输出](#并行与流式输出通配符批量)。
- `--read-type` 与 `--r1`/`--r2` 路径完全同规则解析：内容检测为 FASTQ 选 `sr`，除非
  你显式指定预设。
- FASTA 与 FASTQ 位置参数混用是用法错误（退出码 2），错误信息点名 reads 文件：
  `mixed assembly and reads inputs; screen them separately: ...`。显式
  `--aligner blastn` 或 `--aligner minimap2` 保持各自的既有路由（引擎的冻结守卫生
  效）；通配符自动检测只发生在默认引擎下。
- bzip2 压缩的 FASTQ 文件名不按扩展名分类（reads 引擎只读明文与 gzip）；这类文件
  保持今天走的 contig 流水线。

### 并行与流式输出（通配符批量）

`--jobs N`/`-j N` 并发筛查 N 个样本，契约与 contig 路径完全一致：每个 worker 对
共享数据库跑自己的 minimap2，`--jobs × --threads` 超过 cpu 数时 stderr 提示一行超
订警告；某样本失败时按样本顺序在其位置抛出（退出码与错误封皮和串行运行完全相
同）。结果始终按样本顺序输出 —— `-j 8` 与 `-j 1` 的 stdout 逐字节一致。

输出按完成的样本流式写出，`--jobs` 下按样本顺序队头阻塞发射：

- **tsv**（默认格式）流式输出：先发 `#SAMPLE` 表头行，随后每完成一个样本即落出它
  的基因行（每基因一行，样本名在第一列；无基因的样本不产生行）。`--format csv` 是
  同一张表的逗号拼写。tsv 是 gapit 所有界面上的人类默认值 —— 单样本
  `--r1`/`--r2` 与 `--aligner minimap2` 路径同样默认它；要给 agent 单一文档时改传
  `--format json`。
- **md** 先发静态 frontmatter（schema、tool、`created_at`、db、阈值 —— 本路径没有
  运行总计；总计在 JSON 文档里），随后每完成一个样本即发出它的 `## <样本>` 小节。
- **json**（本路径的可选项）保持为单个文档、最后一次性写出 —— 与 contig 路径的
  json 规则相同。
- `--output PATH` 增量写入流式分块（每块即刻 flush；stdout 不再携带数据）。若第
  *k* 个样本批量中途失败，已流式写出的样本 1..*k-1* 会保留在文件里 —— json 变体
  什么都不写，因为它的单文档只在结尾才存在。

单样本的 `--r1`/`--r2` 路径在两点上不同：不接受 `--jobs`，其 Markdown
frontmatter 保留 `files:`/`genes_found:` 运行总计（它的 tsv 默认渲染同一张表，
单个样本一个分块）。

## 选项

reads 模式走同一条 `gapit screen` 命令；以下是适用的参数（转写自 `gapit screen
--help`，gapit 0.5.2）。这里未列出的 contig 模式参数（`--minid`、`--mincov`、
`--fofn`、`--noheader`、`--nopath`）不适用；`--jobs` 只适用于通配符路径（见上文）。

| 参数 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `--r1` | str | 必填 | reads 或 assembly FASTA 文件，逗号分隔，每条 lane 一个。 |
| `--r2` | str | 无 | 逗号分隔的 mate FASTQ 文件；数量必须与 `--r1` 一致。 |
| `--read-type` | sr\|map-ont\|map-hifi | FASTQ 用 `sr`，FASTA 用 `map-ont` | minimap2 预设；省略时按检测到的输入解析。 |
| `--min-breadth` | float | `90.0` | 判定存在的最小广度百分比。 |
| `--min-identity` | float | `0.0`（关闭） | 单条比对的最小一致性百分比（0-100）；任意非零值开启 gapit.reads/2 过滤（见下文）。 |
| `--min-mapq` | int | `0`（关闭） | 单条比对的最小 MAPQ；任意非零值开启 gapit.reads/2 过滤。 |
| `--db` | str | 必填 | 用于筛查的数据库（数据目录的子目录）。没有默认值 —— 请显式选择（`gapit db list`）。 |
| `--datadir` | path | `$GAPIT_DATADIR`，然后 `~/.local/share/gapit/db` | 数据库目录。 |
| `--threads` | int | `1` | minimap2 工作线程数。 |
| `--jobs` | int | `1` | 仅通配符路径：并发筛查 N 个样本（`--r1`/`--r2` 拒绝它）。输出顺序始终是样本顺序。 |
| `--quiet` | 开关 | 关闭 | 静默 stderr 诊断（包括 assembly-FASTA 提示）。 |
| `--debug` | 开关 | 关闭 | 详细的 stderr 诊断；回显 minimap2 命令行。 |
| `--format` | tsv\|csv\|json\|md | `tsv` | 输出格式。所有路径的 reads 默认输出流式 tsv 表格（csv = 逗号拼写）；`--format json` 是给 agent 的单个 `gapit.reads/1`(+`/2`) 文档，`--format md` 是 Markdown 形式。 |

`--minid` 和 `--mincov` 不适用于 reads 模式。存在与否只由广度决定。

## 选择 `--read-type`

| 测序平台 | 预设 |
|---|---|
| Illumina（短 reads，双端或单端） | `sr` |
| Oxford Nanopore | `map-ont` |
| PacBio HiFi | `map-hifi` |

省略 `--read-type` 时，预设按检测到的输入解析：FASTQ 用 `sr`（历史默认值），
assembly FASTA 用 `map-ont`，并在 stderr 用一行提示宣布
（`assembly FASTA detected; using map-ont`）。与输入矛盾的显式预设是用法错误：
assembly FASTA 要求 `map-ont`（退出码 2）。FASTQ 接受任何显式预设；测长 reads 时请
设置它，否则比对质量会受损。

## 指标

对每个基因（数据库序列，长度 `tlen`），统计其上 mapped reads 的主比对：

| 指标 | 含义 |
|---|---|
| `breadth_pct` | `100 * covered_bases / tlen`，其中覆盖碱基是全部主比对受试跨度的并集。 |
| `mean_depth` | `sum(per-base depth) / tlen`，在全基因长度上取平均。 |
| `reads_mapped` | 在该基因上至少有一条主比对的去重 read 名数量。 |

`breadth_pct >= --min-breadth`（默认 90.0，srst2 风格）时基因判定为 `present`。
零 mapped reads 的基因完全不进输出。

调参建议：默认 90% 广度是刻意从严的。用 `sr` 筛查极短的参考（大约 100 nt 以下）很
难达到，因为 minimap2 会在每个比对末端软裁剪几个碱基，所以小型自定义数据库应把
`--min-breadth` 调低。噪声大的长 reads 或许也需要调低。逐 read 的一致性与 MAPQ 过
滤是可选项，见 `--min-identity` / `--min-mapq`（下一节）。

## 输出

默认 reads 表只列出**present 基因**（`PRESENT` 列恒为 `yes`）；`--all-genes` / `-A`
加入缺席判定（breadth 低于阈值）。表尾列为 `PRODUCT`（contig 表保留 abricate 的
`RESISTANCE` 列；reads 表不携带该列）。


默认是流式 **tsv** 表格（tsv 是 gapit 所有界面上的人类默认值；`--format csv` 是它的
逗号拼写）：

```console
$ gapit screen --r1 tetx_full.fq --db tinyreads
Screening reads: tetx_full.fq
Detected 1 present genes in tetx_full.fq
#SAMPLE	GENE	BREADTH%	DEPTH	READS	PRESENT	DATABASE	ACCESSION	PRODUCT
tetx_full.fq	tetX	97.70	2.09	12	yes	tinyreads	SYN-001	extended resistance determinant tetX	TETRACYCLINE
```

（`Screening reads:`/`Detected` 行走 stderr；stdout 是表格。）

`--format json` 生成版本化的 agent 文档，schema `gapit.reads/1`：

- `files[]` 与输入一一对应：每项列出它筛查的 reads 和检出的基因。
- 基因条目按 `breadth_pct` 降序、再按基因名排序。
- 用 `gapit schema reads` 自省 schema；细节见
  [./outputs.md](./outputs.md)。

`--format md` 生成带 YAML frontmatter 的 Markdown 形式。

## 按一致性和 MAPQ 过滤比对（gapit.reads/2）

reads 模式的存在判定默认只看广度，这会对同源基因家族过量判定：新等位基因的 reads
以低一致性主比对堆到每一条相似的数据库条目上，广度不断累积，直到错误的基因越过阈
值。KP 基准测试让这个失败变得具体：blastn contig 筛查确认了 18 个基因，reads 流水
线在 Illumina `sr` 下只判定出 11 个（拆分比对让家族成员双双挨饿），在 ONT 下判定出
140 个（噪声 reads 几乎全部过量判定）。`--min-identity` 和 `--min-mapq` 是可选的修
正手段：

- `--min-identity FLOAT`（0-100，默认 0 = 关闭）：只保留 `identity >= 阈值` 的比对。
  单条比对的一致性是 `100 * (alen - nm) / alen`，基于 PAF block 长度和 `NM:i:` 错配
  数；没有 NM 的行按 100% 计（无法评估）。
- `--min-mapq INT`（默认 0 = 关闭）：只保留 PAF MAPQ >= 阈值的比对，对付把广度摊到
  重复基因拷贝上的多比对 reads。
- 两个过滤器都在 primary-only 规则之后、聚合之前运行；广度、深度和 `reads_mapped`
  只在幸存比对上计算。丢失全部比对的基因完全不进输出。
- 两个阈值都关闭时一切照旧：运行保持 `gapit.reads/1`，逐字节一致，minimap2 调用也
  不变。任何非零值会把输出切换成 **`gapit.reads/2`** 文档：形状相同，外加
  `params.min_identity` / `params.min_mapq` 和逐基因的 `mean_identity_pct`（幸存比
  对按比对长度加权的一致性均值）。用 `gapit schema reads2` 自省。`--format md` 下
  frontmatter 多出这两个阈值，每个基因行多出 `Identity%` 列。
- 一个细节：`/2` 调用会给 minimap2 传 `--cs`（唯一会发出 `NM:i:` 的 PAF 侧参数），
  两种几何下 minimap2 发出的跨度可能略有差异：同一条 reads 的未过滤 `/1` 与过滤
  `/2` 运行之间，广度会有小幅出入。
- 这些参数仅限 reads 引擎：与 blastn contig 流水线组合是用法错误（退出码 2）。

经验法则：`--min-identity 95` 近似等位基因级严格度（真基因的比对在 ~98-100% 存活，
同源基因 ~85-92% 被丢弃）；原始 ONT reads 的真比对更嘈杂，可降到 90 左右。数据库里
有重复或倍增的基因拷贝、广度在拷贝间分摊时，用 `--min-mapq 20` 或更高。

### 实战示例：同源夹具

仓库夹具确定性地复现了这个失败与修正
（`tests/data/reads2_db/homologs` + `tests/data/reads2`）：数据库有一条 600 nt 的
`geneA` 加一条 400 nt 的部分同源 `geneB`（一致性约 85%，样本中不存在）；样本携带
`geneA` 和一条新的 B 类等位基因，其 reads 以约 90% 一致性落在 `geneB` 上。ONT 风格
reads，未过滤，两个基因都判定为存在：

```console
$ mkdir -p /tmp/gapit-demo/readdb
$ cp -r tests/data/reads2_db/homologs /tmp/gapit-demo/readdb/
$ export GAPIT_DATADIR=/tmp/gapit-demo/readdb
$ cd tests/data/reads2
$ gapit screen --r1 ont_homologs.fq --db homologs --read-type map-ont --format json --quiet
{
  "schema": "gapit.reads/1",
  ...
  "files": [
    {
      "reads": [
        "ont_homologs.fq"
      ],
      "genes": [
        {
          "gene": "geneB",
          ...
          "tlen": 400,
          "breadth_pct": 98.75,
          "mean_depth": 14.03,
          "reads_mapped": 15,
          "present": true
        },
        {
          "gene": "geneA",
          ...
          "tlen": 600,
          "breadth_pct": 98.5,
          "mean_depth": 14.43,
          "reads_mapped": 15,
          "present": true
        }
      ]
    }
  ]
}
```

`geneB` 是误判，样本并不携带它。加 `--min-identity 95` 丢掉 `geneB` 上所有约
89% 一致性的比对（该基因消失；零 read 基因不输出），保留真基因约 98% 的比对：

```console
$ gapit screen --r1 ont_homologs.fq --db homologs --read-type map-ont --min-identity 95 --format json --quiet
{
  "schema": "gapit.reads/2",
  ...
  "params": {
    "db": "homologs",
    "read_type": "map-ont",
    "min_breadth": 90.0,
    "threads": 1,
    "min_identity": 95.0,
    "min_mapq": 0
  },
  "files": [
    {
      "reads": [
        "ont_homologs.fq"
      ],
      "genes": [
        {
          "gene": "geneA",
          ...
          "breadth_pct": 100.0,
          "mean_depth": 14.98,
          "reads_mapped": 15,
          "present": true,
          "mean_identity_pct": 98.04
        }
      ]
    }
  ]
}
```

短 reads 夹具行为相同，只需一次校准：`sr` 的软裁剪把 `geneB` 未过滤广度压在约
89.8%，所以未过滤一组用 `--min-breadth 80` 跑（`geneA` 94.83% / `geneB` 89.75%，
两者存在；过滤后：`geneA` 96.5%、一致性 100.0，`geneB` 消失）。过滤运行的
Markdown 形式：

```console
$ gapit screen --r1 ont_homologs.fq --db homologs --read-type map-ont --min-identity 95 --format md --quiet
---
schema: gapit.reads/2
tool: gapit 0.5.2
created_at: 2026-09-20T14:23:33Z
db: homologs
read_type: map-ont
min_breadth: 90.0
threads: 1
min_identity: 95.0
min_mapq: 0
files: 1
genes_found: 1
---

# gapit read screening report

## `ont_homologs.fq`

| Gene | Breadth% | Depth | Reads | Present | Database | Accession | Product | Resistance | Identity% |
|---|---|---|---|---|---|---|---|---|---|
| geneA | 100.00 | 14.98 | 15 | yes | homologs | SYN-A | true allele carried by the sample | TETRACYCLINE | 98.04 |
```

## 按基因设置一致性下限（gapit.floors/1，数据库侧）

`--min-identity` 是每次运行时由调用方做的决定。有时这个决定应该属于数据库本身：
一张精选面板知道自己的哪些基因身处同源家族之中、真等位基因必须达到多少一致性。
这类知识以可选的 `floors.json` 边车文件随数据库发布（与 `sequences` 同目录），
schema 为 `gapit.floors/1`（用 `gapit schema floors` 内省）：

```json
{"schema": "gapit.floors/1", "default": null, "genes": {"pic": 90.0}}
```

- `genes` 把基因名映射到其最小比对一致性（%）；`default` 作用于未列出的基因
  （`null` = 无下限）。没有该文件就完全没有下限。
- 在 reads 模式下，落在带下限基因上、逐比对一致性（与 `--min-identity` 相同的
  `100 * (alen - nm) / alen` 规则）低于下限的比对会在聚合**之前**被丢弃：广度、
  深度与 `reads_mapped` 都从幸存行重新计算——以 86% 一致性堆上基因、广度 97% 的
  同源序列因此塌缩到存在阈值之下，不再误报。
- 下限是**数据库驱动的可选开启**：CLI 契约没有任何变化，不带边车的数据库输出与
  之前逐字节一致——既有数据库和已发布的输出 schema（`gapit.reads/1`、
  `gapit.reads/2`）都不受影响。blastn contig 管线从不读取 floors；它的
  `--minid` 本就按一致性过滤。
- floors 与 `--min-identity` 是合取关系：两者同时生效时，一条比对必须同时通过
  基因下限和全局过滤（逐行取更严者）。
- 安装方式：`gapit db build NAME genes.fa --floors floors.json`——先校验
  （`genes` 中的每个基因必须存在于 FASTA、取值在 [0, 100] 内），再原样拷入数据库
  目录成为 `floors.json`。把 floors.json 直接放进既有数据库目录同样有效；筛查时
  读取。详见
  [databases.zh.md](./databases.zh.md#floorsjson--按基因设置一致性下限)。

### 实战示例：`pic` 的 SPATE 假阳性（内置 `ecoli_dec`）

分离株 26ECO0084 携带一个与面板 `pic` 基因约 86.6% 一致的 SPATE 家族同源序列。
reads 模式下，该同源序列的比对仍覆盖了基因长度的 97.59%，广度判存在因此报出了
`pic`——而 blastn contig 管线在默认 80% 一致性下限处正确地把它过滤掉了。于是内置
的 `ecoli_dec` 携带了仅一条下限的 `floors.json`：

```json
{"schema": "gapit.floors/1", "default": null, "genes": {"pic": 90.0}}
```

真 `pic` 等位基因一致性 ≥95%，所以 90 的下限只移除 SPATE 同源序列的低一致性比对
——源自同源序列的 reads 里 `pic` 翻转为不存在，而真等位基因的判定不变。`astA`
刻意不设下限：其真等位基因一致性在 90.4–91% 之间，90 的下限之下毫无余量。

## 实战示例

仓库测试夹具还包含一个微型 reads 数据库（`tinyreads`，一条 522 nt 的 `tetX` 基因加
一条部分覆盖的 `sulY`）和几个小 FASTQ 文件。把数据库复制到临时数据目录，在 reads
目录下运行：

```console
$ mkdir -p /tmp/gapit-demo/readdb
$ cp -r tests/data/reads_db/tinyreads /tmp/gapit-demo/readdb/
$ export GAPIT_DATADIR=/tmp/gapit-demo/readdb
$ cd tests/data/reads
$ gapit screen --r1 tetx_full.fq --db tinyreads
Screening reads: tetx_full.fq
Detected 1 present genes in tetx_full.fq
#SAMPLE	GENE	BREADTH%	DEPTH	READS	PRESENT	DATABASE	ACCESSION	PRODUCT
tetx_full.fq	tetX	97.70	2.09	12	yes	tinyreads	SYN-001	extended resistance determinant tetX	TETRACYCLINE
```

（`Screening reads:` 和 `Detected N present genes` 行走 stderr；stdout 是表格。
`--format json` 换成版本化文档：

```console
$ gapit screen --r1 tetx_full.fq --db tinyreads --format json
{
  "schema": "gapit.reads/1",
  "tool": {
    "name": "gapit",
    "version": "0.5.2"
  },
  "created_at": "2026-09-19T01:12:25Z",
  "params": {
    "db": "tinyreads",
    "read_type": "sr",
    "min_breadth": 90.0,
    "threads": 1
  },
  "files": [
    {
      "reads": [
        "tetx_full.fq"
      ],
      "genes": [
        {
          "gene": "tetX",
          "database": "tinyreads",
          "accession": "SYN-001",
          "product": "extended resistance determinant tetX",
          "resistance": "TETRACYCLINE",
          "tlen": 522,
          "breadth_pct": 97.7,
          "mean_depth": 2.09,
          "reads_mapped": 12,
          "present": true
        }
      ]
    }
  ]
}
```

### 双端与多条 lane

一条双端 lane：传入两条 mate。两条单端 lane 合成一个样本：用逗号串接。两者在这里
聚合成的相同结果：

```console
$ gapit screen --r1 tetx_R1.fq --r2 tetx_R2.fq --db tinyreads --format json --quiet
{
  "schema": "gapit.reads/1",
  ...
  "files": [
    {
      "reads": [
        "tetx_R1.fq",
        "tetx_R2.fq"
      ],
      "genes": [
        {
          "gene": "tetX",
          ...
          "breadth_pct": 97.7,
          "mean_depth": 2.09,
          "reads_mapped": 12,
          "present": true
        }
      ]
    }
  ]
}
$ gapit screen --r1 tetx_lane1.fq,tetx_lane2.fq --db tinyreads --format json --quiet
{
  "schema": "gapit.reads/1",
  ...
  "files": [
    {
      "reads": [
        "tetx_lane1.fq",
        "tetx_lane2.fq"
      ],
      "genes": [
        {
          "gene": "tetX",
          ...
          "breadth_pct": 97.7,
          "mean_depth": 2.09,
          "reads_mapped": 12,
          "present": true
        }
      ]
    }
  ]
}
```

（`...` 表示与上方第一个 JSON 块相同的输出。）

### Markdown 输出

```console
$ gapit screen --r1 tetx_full.fq --db tinyreads --format md
---
schema: gapit.reads/1
tool: gapit 0.5.2
created_at: 2026-09-19T01:12:25Z
db: tinyreads
read_type: sr
min_breadth: 90.0
threads: 1
files: 1
genes_found: 1
---

# gapit read screening report

## `tetx_full.fq`

| Gene | Breadth% | Depth | Reads | Present | Database | Accession | Product | Resistance |
|---|---|---|---|---|---|---|---|---|
| tetX | 97.70 | 2.09 | 12 | yes | tinyreads | SYN-001 | extended resistance determinant tetX | TETRACYCLINE |
```

### 调节 `--min-breadth`

`suly_partial.fq` 覆盖 261 nt `sulY` 基因的 63.98%。默认阈值下基因以 `present:
false` 报告；降低阈值即可翻转判定，指标不受影响：

```console
$ gapit screen --r1 suly_partial.fq --db tinyreads --format json --quiet
{
  "schema": "gapit.reads/1",
  ...
  "genes": [
    {
      "gene": "sulY",
      ...
      "tlen": 261,
      "breadth_pct": 63.98,
      "mean_depth": 1.0,
      "reads_mapped": 3,
      "present": false
    }
  ]
}
$ gapit screen --r1 suly_partial.fq --db tinyreads --min-breadth 50 --format json --quiet
{
  "schema": "gapit.reads/1",
  ...
  "genes": [
    {
      "gene": "sulY",
      ...
      "tlen": 261,
      "breadth_pct": 63.98,
      "mean_depth": 1.0,
      "reads_mapped": 3,
      "present": true
    }
  ]
}
```

### 默认表格

默认运行就是这张表（tsv 是 gapit 所有界面上的人类默认值；`--format csv` 是逗号拼
写）：

```console
$ gapit screen --r1 tetx_full.fq --db tinyreads
#SAMPLE	GENE	BREADTH%	DEPTH	READS	PRESENT	DATABASE	ACCESSION	PRODUCT
tetx_full.fq	tetX	97.70	2.09	12	yes	tinyreads	SYN-001	extended resistance determinant tetX	TETRACYCLINE
```

样本名在第一列；agent 需要版本化文档时传 `--format json`。

## 筛查 assembly（快速存在性普查）

主要写法是 `gapit screen --aligner minimap2 assembly.fa`：位置文件被送进 minimap2
引擎，且每个输入都必须是 FASTA 内容，gzip 或 plain 均可；FASTQ 文件放这里是用法
错误（退出码 2）。`--r1 assembly.fa` 写法等价。minimap2 原生接受 FASTA query，gapit
与处理 FASTQ 时一样原样传文件；变化只在预设：内容检测强制 `map-ont`（比如一条连续
的 522 nt contig 在 `sr` 下只能比对到基因的约 14%，短 read 软裁剪会毁掉长 query 的
比对），gapit 会在 stderr 上说明。每条 contig 相当于一条长 read：`reads_mapped` 统
计 contig 数，`mean_depth` 在被覆盖区间附近徘徊，`present` 仍表示
`breadth_pct >= --min-breadth`。文档保持 `gapit.reads/1`；`params.read_type` 报告解
析出的预设。默认输出是流式 tsv 表格（单样本一个分块）；`--format json` 给出下文展
示的文档。

用 tinyreads 夹具冒充 assembly（任何多 contig FASTA 行为相同）：

```console
$ mkdir -p /tmp/gapit-demo/readdb
$ cp -r tests/data/reads_db/tinyreads /tmp/gapit-demo/readdb/
$ export GAPIT_DATADIR=/tmp/gapit-demo/readdb
$ cp tests/data/reads_db/tinyreads/sequences /tmp/gapit-demo/assembly.fa
$ gapit screen --aligner minimap2 /tmp/gapit-demo/assembly.fa --db tinyreads --format json
assembly FASTA detected; using map-ont
Screening reads: /tmp/gapit-demo/assembly.fa
Detected 2 present genes in /tmp/gapit-demo/assembly.fa
{
  "schema": "gapit.reads/1",
  "tool": {
    "name": "gapit",
    "version": "0.5.2"
  },
  "created_at": "2026-09-19T14:23:00Z",
  "params": {
    "db": "tinyreads",
    "read_type": "map-ont",
    "min_breadth": 90.0,
    "threads": 1
  },
  "files": [
    {
      "reads": [
        "/tmp/gapit-demo/assembly.fa"
      ],
      "genes": [
        {
          "gene": "tetX",
          "database": "tinyreads",
          "accession": "SYN-001",
          "product": "extended resistance determinant tetX",
          "resistance": "TETRACYCLINE",
          "tlen": 522,
          "breadth_pct": 97.89,
          "mean_depth": 0.98,
          "reads_mapped": 1,
          "present": true
        },
        {
          "gene": "sulY",
          "database": "tinyreads",
          "accession": "SYN-002",
          "product": "partial coverage test determinant sulY",
          "resistance": "SULFONAMIDE",
          "tlen": 261,
          "breadth_pct": 95.79,
          "mean_depth": 0.96,
          "reads_mapped": 1,
          "present": true
        }
      ]
    }
  ]
}
```

（`assembly FASTA detected`、`Screening reads:`、`Detected ...` 行走 stderr；
`--quiet` 把提示和杂音一并静默，stdout 逐字节不变。）

### 两段式模式：先普查，后确认

minimap2 一段完全跳过 BLAST 索引，把 assembly 送进 minimap2 引擎比 BLAST contig 流
水线快大约一个数量级，代价是等位基因级精度。这个取舍适合对大量样本做两段式流程：

1. **普查**：用 minimap2 加刻意放宽的广度下限筛查每个样本。
2. **确认**：只对阳性（或任何有命中的样本）跑常规 contig 流水线，按 abricate 一致
   性应用一致性与覆盖度下限。

真实数字，一株 K. pneumoniae RefSeq assembly（GCF_000240185.1，5.3 Mb，`--db ncbi`，
单线程，gapit 0.5.2）：

```console
$ # Stage 1: survey, ~0.9 s
$ gapit screen --aligner minimap2 kpneu_mgh78578.fna.gz --db ncbi --min-breadth 50 --format json --quiet
{
  "schema": "gapit.reads/1",
  ...
  "params": {
    "db": "ncbi",
    "read_type": "map-ont",
    "min_breadth": 50.0,
    "threads": 1
  },
  "files": [
    {
      "reads": ["kpneu_mgh78578.fna.gz"],
      "genes": [
        { "gene": "aph(3'')-Ib", "breadth_pct": 99.75, "present": true, ... },
        { "gene": "dfrA12", "breadth_pct": 99.6, "present": true, ... },
        { "gene": "tet(G)", "breadth_pct": 99.57, "present": true, ... },
        { "gene": "blaKPC-2", "breadth_pct": 99.55, "present": true, ... },
        ... 11 more present genes (floR2 dfrA50 blaCTX-M-14 rmtB1 aadA2 sul2 blaTEM-1 aph(6)-Id fosA6 blaSHV-155 aac(3)-IId) ...
        { "gene": "sul1", "breadth_pct": 62.74, "present": true, ... },
        { "gene": "tmexD3", "breadth_pct": 60.22, "present": true, ... },
        { "gene": "tmexD2", "breadth_pct": 56.2, "present": true, ... }
      ]
    }
  ]
}
$ # Stage 2: confirm the positives with the contig pipeline, ~13.6 s
$ gapit screen kpneu_mgh78578.fna.gz --db ncbi --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
kpneu_mgh78578.fna.gz	NC_016838.1	31843	32718	-	blaCTX-M-14	1-876/876	===============	0/0	100.00	100.00	ncbi	NG_048929.1	extended-spectrum class A beta-lactamase CTX-M-14	CEPHALOSPORIN
kpneu_mgh78578.fna.gz	NC_016838.1	108291	108766	-	dfrA50	2-477/477	===============	0/0	99.79	98.74	ncbi	NG_242637.1	trimethoprim-resistant dihydrofolate reductase DfrA50	TRIMETHOPRIM
... 16 more rows: aac(3)-IId aadA2 aph(3'')-Ib aph(6)-Id blaKPC-2 blaSHV-158 blaTEM-1 dfrA12 floR2 fosA6 rmtB1 sul2 tet(G) ...
```

（`...` 省略字段和行；普查 JSON 裁剪到下文讨论的基因。）

在这个基因组上普查耗时 0.9 秒，BLAST 13.6 秒（约 17 倍），18 个确认基因中 14 个同
时出现在两个名单里。差异正是重点：

- **家族级分辨率，不是等位基因级。** 普查把 blaSHV 位点标成 `blaSHV-155`；BLAST 确
  认位点但等位基因判定为 `blaSHV-158`。minimap2 的 primary-only 分配把共享 contig
  交给近缘家族成员中的一个，所以把普查基因名当作家族级线索，等位基因命名交给确认
  阶段。
- **没有一致性下限。** reads 模式的存在判定只看广度。`sul1`（62.7% 广度）和
  `tmexD2`/`tmexD3` 对（56-60%）过了放宽的 50% 普查下限，但它们是部分或分歧位点，
  确认阶段的 80/80 一致性/覆盖度阈值会拒绝。放宽普查下限是刻意用精度换召回；精确
  过滤属于第二阶段——唯一的例外是数据库带有
  [按基因下限](#按基因设置一致性下限gapitfloors1数据库侧) 的基因，普查同样尊重它。

需要一次拿到精确等位基因时，直接用 contig 流水线；需要快速从大量 assembly 得到基
因家族层面的答案时，用普查。

## 注意事项

- **只统计主比对。** 每条 read 对每个基因至多贡献一条比对（`tp:A:P` 记录；没有
  `tp` 标签的记录也算）。次级比对从不抬高广度或深度。
- **同源基因家族。** 近缘等位基因争夺同一批 reads。primary-only 分配会把共享 reads
  错分给近缘家族成员，所以单个成员的广度可能偏低而整个家族覆盖良好，新等位基因的
  reads 也可能把缺席的同源基因推过广度下限。`--min-identity`（上一节）能去掉这类失
  败的低一致性一侧。v1 没有 SNP 级等位基因判定。
- **极短的基因。** `sr` 预设下，比对末端的软裁剪限制了可达广度；100 nt 以下的基因
  可能永远到不了 90%。reads 模式面向正常长度的基因（数百 nt 及以上）；微型数据库请
  调低 `--min-breadth`。
- **Gzip 输入。** FASTQ 和 FASTA 都支持 `.gz`；输入检测透过 gzip 包装窥探，minimap2
  原生解压。
