# 快速开始

一次完整的首个会话，离线、可复现，使用测试套件自带的微型夹具数据库。你将筛查两个
contig 文件，以 TSV 和 JSON 两种形式读结果，再把两份报告汇总成一个矩阵。

一切都在仓库根目录运行，并让 gapit 环境在 PATH 上
（`export PATH="$PWD/.pixi/envs/default/bin:$PATH"`，或给每条命令加 `pixi run`
前缀）。

## 1. 准备一个用完即弃的数据目录

夹具数据库 `tests/data/db/tinyamr` 里有三个短 AMR 基因：`tetA`、`blaTEM-1` 和
`sul1`。把它复制进一个临时数据目录并构建 BLAST 索引，这正是测试套件为 MCP 夹具
做的事：

```bash
mkdir -p /tmp/gapit-quickstart/db
cp -r tests/data/db/tinyamr /tmp/gapit-quickstart/db/
makeblastdb -in /tmp/gapit-quickstart/db/tinyamr/sequences \
  -title tinyamr -dbtype nucl -logfile /dev/null
export GAPIT_DATADIR=/tmp/gapit-quickstart/db
```

`gapit db list` 展示数据库目录和安装状态（下面有删节；`--json` 返回
`gapit.dblist/1` 文档）。夹具 `tinyamr` 是一个普通的自定义数据库，不在目录里，所以
不会出现在列表里，下一节的筛查会证明它可用：

```console
$ gapit db list
NAME	PROVIDER	STATUS	DBTYPE	DESCRIPTION
argannot	IHU Méditerranée-Infection	available	nucl	ARG-ANNOT acquired resistance genes
ncbi	NCBI	available	nucl	NCBI AMRFinderPlus (reference finder) curated AMR
```

（共十二个目录数据库，在这个临时数据目录里全部为 `available`。）

## 2. 筛查一个 contig 文件

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
```

TSV 与 abricate 兼容，可以直接丢进任何已在解析 abricate 的流水线。`Processing:`
进度行走 stderr；stdout 只承载数据。各列含义：

| 列 | 含义 |
|---|---|
| `#FILE` | 输入文件，按命令行上给定的路径 |
| `SEQUENCE` | 命中所在的 contig |
| `START`、`END` | contig 上的命中坐标（1 基） |
| `STRAND` | 基因方向，`+` 或 `-` |
| `GENE` | 数据库头中的基因名 |
| `COVERAGE` | 基因上的比对跨度，`range/length` |
| `COVERAGE_MAP` | 基因的定宽覆盖度草图；`=` 表示覆盖，`/` 标记一段缺口 |
| `GAPS` | 比对中的缺口开口数/缺口碱基数 |
| `%COVERAGE` | 基因被覆盖的比例（内部不取整，显示两位小数） |
| `%IDENTITY` | BLAST 一致性百分比，显示两位小数 |
| `DATABASE` | 筛查所用数据库 |
| `ACCESSION` | 数据库头中的上游登录号 |
| `PRODUCT` | 基因产物描述 |
| `RESISTANCE` | 数据库头中的耐药类别 |

默认值沿用 abricate：命中需通过 `--minid 80` 和 `--mincov 80`（一致性百分比、基因
覆盖度百分比）。所有参数与过滤规则见[筛查](./screen.md)。

## 3. 同一结果的 JSON 形式

```bash
gapit screen tests/data/contigs/full.fa --db tinyamr --format json
```

```json
{
  "schema": "gapit.report/1",
  "tool": {
    "name": "gapit",
    "version": "0.5.1"
  },
  "created_at": "2026-09-19T01:11:06Z",
  "params": {
    "db": "tinyamr",
    "minid": 80.0,
    "mincov": 80.0,
    "threads": 1
  },
  "files": [
    {
      "file": "tests/data/contigs/full.fa",
      "hits": [
        {
          "sequence": "contig1",
          "start": 1,
          "end": 79,
          "strand": "+",
          "gene": "tetA",
          "coverage": "1-79/79",
          "coverage_map": "===============",
          "gaps": "0/0",
          "coverage_pct": 100.0,
          "identity_pct": 100.0,
          "database": "tinyamr",
          "accession": "NC_000913.3:100-900",
          "product": "tetracycline efflux pump TetA",
          "resistance": "TETRACYCLINE"
        }
      ]
    }
  ]
}
```

文档携带 `"schema": "gapit.report/1"`，消费方可以锁定契约；`gapit schema report`
打印它的 JSON Schema。字段级细节见[输出](./outputs.md)。

## 4. 把报告汇总成矩阵

保存两份报告，再折叠成一个"基因 × 文件"矩阵。第二个夹具文件 `gap.fa` 带有
`sul1`，比对中有缺口：

```bash
gapit screen tests/data/contigs/full.fa --db tinyamr > full.tsv
gapit screen tests/data/contigs/gap.fa --db tinyamr > gap.tsv
gapit summary full.tsv gap.tsv
```

```console
#FILE	NUM_FOUND	sul1	tetA
full.tsv	1	.	100.00
gap.tsv	1	100.00	.
```

每行对应一份报告文件；列是检出的基因并集。单元格存 %COVERAGE（缺失时为 `.`），
`NUM_FOUND` 统计不同基因的个数。加 `--format json|md` 得到机器或人类可读的矩阵；
细节见[汇总报告](./summary.md)。

## 5. 试试基因簇引擎

基因**簇**数据库是第二种数据库：记录不再是单个基因，而是整个位点（Kaptive 风格的
抗原位点、荚膜基因簇），从 GenBank/GFF 输入构建。夹具
`tests/data/cluster/screening.gbk` 含两个合成位点；`--typing` 会同时安装一份
`gapit.typing/1` 表型评分规范：

```bash
gapit db build tinykps tests/data/cluster/screening.gbk \
  --datadir /tmp/gapit-quickstart/db \
  --typing tests/data/cluster/typing_screen.json
# 从数据库里取出一个位点，作为单位点 query 装配体
awk '/^>locusA$/{p=1} /^>/{if($0!~/>locusA$/)p=0} p' \
  /tmp/gapit-quickstart/db/tinykps/sequences > /tmp/gapit-quickstart/locusA.fa
gapit screen /tmp/gapit-quickstart/locusA.fa --db tinykps
```

筛查基因簇数据库会分发到 minimap2 基因簇引擎；TSV 的形状变为每文件一个最佳位点
判定，typing 规范再把它转成表型：

```console
Processing: /tmp/gapit-quickstart/locusA.fa
Best locus in /tmp/gapit-quickstart/locusA.fa: locusA
FILE	BEST_LOCUS	TYPE	PHENOTYPE	COVERAGE	IDENTITY	PRESENT	PARTIAL	MISSING_IDS
/tmp/gapit-quickstart/locusA.fa	locusA	KL101	K101	100.00	100.00	3	0	-
```

`PHENOTYPE` 为 `-` 表示没有表型越过规范的 cutoff，或两个表型打平在歧义容差之内
（把两个位点一起筛查，K101/K102 打平正是这种情况）。七个 kaptive 提供商
（`gapit db fetch kpsc_k`、`kpsc_o`、`kosc_k`、`kosc_o`、`ab_k`、`ab_o`、
`ecoli_kps`）会以同样的方式安装真实的
Kaptive 位点数据库。完整细节见[筛查](./screen.md)与[数据库](./databases.md)。

## 下一步

- [筛查](./screen.md)：contig 模式的全部参数、阈值、多输入、`--jobs`
- [reads（FASTQ）筛查](./reads.md)：FASTQ 输入走 minimap2
- [数据库](./databases.md)：安装真实数据库（`gapit db fetch`）、提供商、数据目录、基因簇数据库与 typing 规范
- [输出](./outputs.md)：格式、schema、错误信封、退出码
- [MCP 服务器](./mcp.md)与[面向 agent 的指南](./agents.md)：从 agent 驱动 gapit
