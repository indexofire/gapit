# 数据库

gapit 数据库是一个参考基因目录，你用 contigs 和 reads 去筛查它。本页介绍数据库存放
在哪里、如何获取、gapit 在磁盘上写了什么，以及如何带上自己的数据库。筛查如何使用
它们见 [screen.md](./screen.md)。

## 数据库存放在哪里

每个数据库位于一个**数据目录**（datadir）下，一个数据库一个子目录。gapit 按以下顺
序解析数据目录：

1. 单次调用上的 `--datadir`
2. `$GAPIT_DATADIR`
3. `~/.local/share/gapit/db`

解析出的路径必须存在，否则命令以类型化错误失败并退出码 4：

```console
$ gapit db list --datadir /no/such/dir
{"schema":"gapit.error/1","code":"DATADIR_NOT_FOUND","message":"datadir does not exist: /no/such/dir","context":{"datadir":"/no/such/dir"}}
```

## 查看已安装的数据库

`gapit db list` 列出所有已知提供商及其安装状态：

```console
$ gapit db list
PROVIDER	STATUS	DBTYPE	DESCRIPTION
argannot	installed (2224)	nucl	ARG-ANNOT acquired resistance genes
bacmet2	installed (746)	prot	BacMet2 experimentally confirmed biocide/resistance genes (protein)
card	installed (6059)	nucl	CARD protein homolog resistance models
ecoh	installed (597)	nucl	E. coli O and H antigens (srst2 EcOH)
ecoli_vf	installed (2701)	nucl	E. coli virulence factors (phac-nml)
megares	installed (7425)	nucl	MEGARes antimicrobial resistance genes
ncbi	installed (8373)	nucl	NCBI AMRFinderPlus (reference finder) curated AMR
plasmidfinder	installed (488)	nucl	CGE PlasmidFinder replicons
resfinder	installed (3206)	nucl	CGE ResFinder acquired resistance genes
upec_expec_vf	installed (77)	nucl	UPEC/ExPEC virulence genes (FordeGenomics)
vfdb	installed (4769)	nucl	VFDB virulence factors (set A, nucleotide)
victors	installed (4402)	nucl	Victors virulence factors
```

当 `<datadir>/<name>/gapit-manifest.json` 存在时 STATUS 显示 `installed (N)`，N 是
记录数，否则显示 `available`。面向 agent，`--json` 输出 `gapit.dblist/1` 文档
（十二个提供商中的前三个，输出有删节）：

```console
$ gapit db list --json
{
  "schema": "gapit.dblist/1",
  "providers": [
    {
      "name": "argannot",
      "description": "ARG-ANNOT acquired resistance genes",
      "dbtype": "nucl",
      "installed": true,
      "records": 2224
    },
    {
      "name": "bacmet2",
      "description": "BacMet2 experimentally confirmed biocide/resistance genes (protein)",
      "dbtype": "prot",
      "installed": true,
      "records": 746
    },
    {
      "name": "card",
      "description": "CARD protein homolog resistance models",
      "dbtype": "nucl",
      "installed": true,
      "records": 6059
    }
  ]
}
```

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.dblist/1` |
| `providers` | array | 每个提供商一项 |
| `providers[].name` | string | 提供商名称，即传给 `gapit db fetch` 的名字 |
| `providers[].description` | string | 内容简述 |
| `providers[].dbtype` | string | `nucl`（用 blastn 筛查）或 `prot`（用 blastx 筛查） |
| `providers[].installed` | boolean | 数据目录里存在 manifest 时为 true |
| `providers[].records` | integer | 记录数，数据库未安装时省略 |

`gapit db list` 是唯一的列表入口：即上面的提供商目录，`--json` 返回
`gapit.dblist/1` 文档。（早期独立的列表命令及其 schema 已移除；这个入口刻意放弃了
abricate `--list` 的逐字节一致性。）见 [outputs.md](./outputs.md)。

## 检查数据库新鲜度

`gapit db outdated` 报告每个已安装数据库的年龄，并标记两种更新情况。对装满的数据
目录（输出裁剪到十二行中的三行）：

```console
$ gapit db outdated --days 30
NAME	FETCHED_AT	AGE_DAYS	STATUS
argannot	2026-09-17T23:22:06Z	2.58	ok
bacmet2	2026-09-17T23:13:51Z	2.58	ok
card	2026-09-17T23:14:25Z	2.58	ok
...
```

安装已久且被内置快照超越的数据库会同时报告两种标记：

```console
$ gapit db outdated --datadir /tmp/opencode/gapit-outdated-demo
NAME	FETCHED_AT	AGE_DAYS	STATUS
card	2020-01-01T00:00:00Z	2454.55	stale+snapshot-update
```

| 状态 | 含义 |
|---|---|
| `ok` | 足够新且没有更新的内置快照 |
| `stale` | `age_days` 超过 `--days`（默认 90；`--days 0` 把一切都标为 stale） |
| `snapshot-update` | 提供商的内置快照比已安装副本新（card、vfdb） |
| `stale+snapshot-update` | 以上两者同时成立 |

过期只是一份报告，绝不是错误状态：无论多旧命令都退出 0。退出码 4
（`DATADIR_NOT_FOUND`、`DATADIR_EMPTY`）对应缺失的数据目录或没有任何已安装数据库的
目录；无法解析的 `fetched_at` 是 `MANIFEST_MALFORMED`（退出码 5）。面向 agent，
`--json` 输出 `gapit.dboutdated/1` 文档（形如 `gapit.dblist/1` 的 CLI 列表，未注册
到 `gapit schema`）：

```console
$ gapit db outdated --days 3 --json
{
  "schema": "gapit.dboutdated/1",
  "databases": [
    {
      "db": "argannot",
      "fetched_at": "2026-09-17T23:22:06Z",
      "age_days": 2.58,
      "status": "ok",
      "upstream_version": ""
    },
    {
      "db": "card",
      "fetched_at": "2026-09-17T23:14:25Z",
      "age_days": 2.58,
      "status": "ok",
      "upstream_version": ""
    }
  ]
}
```

## 跨数据库检索记录

`gapit db search TERM` 在每个已安装数据库的 `records.jsonl` 事实存储中查找基因：
不做 BLAST，只是快速扫描。默认按大小写不敏感的子串匹配；`--exact` 切换到全字段相
等：

```console
$ gapit db search ctx-m --limit 3
argannot	(Bla)blaCTX-M-1	X92506:63-938		(Bla)blaCTX-M-1	876
argannot	(Bla)blaCTX-M-10	AF255298:1-873		(Bla)blaCTX-M-10	873
argannot	(Bla)blaCTX-M-100	FR682582:1-876		(Bla)blaCTX-M-100	876
$ gapit db search "blaCTX-M-1" --field gene --exact
ncbi	blaCTX-M-1	NG_048897.1	CEPHALOSPORIN	extended-spectrum class A beta-lactamase CTX-M-1	876
$ gapit db search virulence --db vfdb --field function --limit 2
vfdb	AAA92657	AAA92657	virulence	(AAA92657) unknown protein [TraJ (VF0241) - Invasion (VFC0083)] [Escherichia coli]	606
vfdb	AAC38364	AAC38364	virulence	(AAC38364) Orf1 [Ler (VF0189) - Regulation (VFC0301)] [Escherichia coli O127:H6 str. E2348/69]	390
```

行格式为 `DB\tGENE\tACCESSION\tFUNCTION\tPRODUCT\tLENGTH`，按数据库再按文件的顺序流
式输出；`FUNCTION` 把记录的功能类别用 `;` 串接。`--field` 选择
`gene|accession|function|product|any`（默认 `any` 搜索全部，功能类别逐个匹配）。
`--limit N` 截断输出（默认 100；`0` = 不限），截断时在 stderr 提示（`--quiet` 静默
）；零命中退出 0 且 stdout 为空。`--json` 每个命中打印一行 JSON 对象，字段同名转
snake_case（`function` 为数组）：

```console
$ gapit db search "tet(M)" --field gene --exact --json | head -1
{"db":"card","gene":"tet(M)","accession":"AB039845.1:25-1945","function":["tetracycline"],"product":"Tet(M) is a ribosomal protection protein that confers tetracycline resistance. It is found on transposable DNA elements and its horizontal transfer between bacterial species has been documented.","length":1920}
```

全库扫描时，`records.jsonl` 缺失的已安装数据库会被跳过并在 stderr 警告；但显式指
定它（`--db NAME`）会以退出码 4 `DB_INCOMPLETE` 失败；未知的 `--db NAME` 是用法错
误（退出码 2），并列出已安装的数据库。

## 提供商

gapit 自带十九个提供商。`card` 和 `vfdb` 还以内置快照的形式随包发布，零网络安装。
七个 kaptive 提供商是**基因簇**数据库（抓取时下载并构建；GPL 内容绝不内置，见下文
Kaptive 提供商一节）。

| 名称 | 内容 | dbtype |
|---|---|---|
| `ncbi` | NCBI AMRFinderPlus (reference finder) curated AMR（默认数据库） | nucl |
| `card` | CARD protein homolog resistance models | nucl |
| `resfinder` | CGE ResFinder acquired resistance genes | nucl |
| `argannot` | ARG-ANNOT acquired resistance genes | nucl |
| `plasmidfinder` | CGE PlasmidFinder replicons | nucl |
| `megares` | MEGARes antimicrobial resistance genes | nucl |
| `ecoh` | E. coli O and H antigens (srst2 EcOH) | nucl |
| `vfdb` | VFDB virulence factors (set A, nucleotide) | nucl |
| `ecoli_vf` | E. coli virulence factors (phac-nml) | nucl |
| `bacmet2` | BacMet2 experimentally confirmed biocide/resistance genes (protein) | prot |
| `victors` | Victors virulence factors | nucl |
| `upec_expec_vf` | UPEC/ExPEC virulence genes (FordeGenomics) | nucl |
| `kpsc_k` | K. pneumoniae species complex K locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `kpsc_o` | K. pneumoniae species complex O locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `kosc_k` | K. oxytoca species complex K locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `kosc_o` | K. oxytoca species complex O locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `ab_k` | A. baumannii K locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `ab_o` | A. baumannii OC locus — 官方关键字 `ab_o`（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `ecoli_kps` | E. coli group 2+3 capsular polysaccharide loci（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |

蛋白质数据库（`bacmet2`）经 `blastx` 筛查；核苷酸数据库经 `blastn`。基因簇数据库
经 minimap2 基因簇引擎筛查（[screen.md](./screen.md)）。

### Kaptive 提供商（GPL，抓取时下载）

七个 kaptive 提供商封装了 [Kaptive](https://github.com/klebgenomics/Kaptive) 的参考
数据库（Wyres et al., J Clin Microbiol 2020：使用这些数据库的结果请引用
Kaptive）。提供商名称就是 Kaptive v3 数据库文档中的**官方安装关键字**
（[Available databases](https://klebgenomics.github.io/Kaptive/db/overview.html#available-databases)），
每个提供商直接从持续维护的按物种上游仓库的 `main` 头部抓取原始 GenBank 文件 ——
`klebgenomics/KpSC_surface_antigen_loci`（`kpsc_k`、`kpsc_o`）、
`klebgenomics/KoSC-surface-antigen-loci`（`kosc_k`、`kosc_o`）、
`johannajkenyon/Abaumannii_surface_polysaccharide_loci`（`ab_k`、`ab_o`）和
`rgladstone/EC-K-typing`（`ecoli_kps`）—— 因此抓到的数据库跟随上游的最新维护，
而不是冻结的发布版。这些仓库还带有包含上游一致性阈值的 `.toml` 元数据；gapit
在 v1 不抓取它（未来 `typing.json` 的来源）。

数据库内容是 **GPL-3.0**，而 gapit 是 MIT，所以任何 kaptive 内容都不内置：
`gapit db fetch kpsc_k`（以及其余六个关键字）在抓取时下载 GenBank 文件，并
经基因簇流水线构建 `kind: cluster` 数据库。manifest 记录上游 URL、
`GPL-3.0 (database content)` 许可证和引用说明。它们不带 typing 模型（表型判定保持
null）；你得到 kaptive 风格的最佳位点判定，之后可以用 `gapit db build --typing`
重建来安装自己的 `typing.json` 语义。

## 抓取数据库

### 内置快照，零网络

不带名字的 `gapit db fetch` 会从包内快照安装默认集合：先 `card` 后 `vfdb`。完全不
碰网络：快照档案携带 `records.jsonl` 和 manifest，gapit 在本地重建 `sequences` 和
BLAST 索引。重建是确定且快速的，并且与你机器上实际安装的 BLAST 版本匹配。

内置数据不会腐坏：一个定时工作流
（`.github/workflows/snapshot-refresh.yml`）每月从上游重新抓取 card 和 vfdb，记录
有变化时就开一个 pull request。那个 PR 就是评审关卡：数据更新由人签字后才能合并。
GitHub Actions 之外，`gapit db fetch <name> --from-source` 仍是手动走上游的路径。

### 抓取指定提供商

`gapit db fetch <name>` 对单个提供商跑完整流水线。下面是从内置快照安装 `card`
（stderr 进度行，随后 stdout 一行 JSON 回执）：

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo card
gapit: installed card from bundled snapshot card.tar.gz
gapit: generated /tmp/opencode/gapit-dbs-demo/card/sequences
gapit: self-check passed for card
gapit: BLAST index built (nucl)
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/opencode/gapit-dbs-demo/card"}
```

| 回执字段 | 类型 | 含义 |
|---|---|---|
| `db` | string | 提供商名称 |
| `records` | integer | 写入 `records.jsonl` 的记录数 |
| `dbtype` | string | `nucl` 或 `prot` |
| `destination` | string | 安装后的数据库目录 |

没有内置快照的提供商在抓取时从上游下载，需要网络。未知名字在解析数据目录之前就会
被拒绝：

```console
$ gapit db fetch nosuchdb
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"unknown provider: nosuchdb (available: argannot, bacmet2, card, ecoh, ecoli_vf, megares, ncbi, plasmidfinder, resfinder, upec_expec_vf, vfdb, victors)","context":{"provider":"nosuchdb"}}
```

### 重新抓取与强制走上游

已有的数据库目录永远不会被静默覆盖。不带参数重新抓取会以退出码 4 失败：

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo card
{"schema":"gapit.error/1","code":"DB_ALREADY_EXISTS","message":"won't overwrite existing database card (use --force)","context":{"db":"card"}}
```

`--force` 删除后原地重建：

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo --force card
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/opencode/gapit-dbs-demo/card"}
```

`--from-source` 跳过内置快照，即使快照存在也强制走上游下载。反过来，提供商没有快
照档案时，抓取会自动落到网络路径。

## 安装本地文件

`gapit db install SOURCE --sha256 HASH --output TARGET` 是带校验的本地文件安装。它
对提供商、档案、序列内容一无所知：把 SOURCE 流式过 SHA256，与 `--sha256` 比对摘
要，通过后才原子替换 TARGET。校验失败时已有的 TARGET 不受影响。

```console
$ sha256sum src/gapit/data/snapshots/card.tar.gz
65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025  src/gapit/data/snapshots/card.tar.gz
$ gapit db install src/gapit/data/snapshots/card.tar.gz \
    --sha256 65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025 \
    --output /tmp/opencode/card-copy.tar.gz
{"destination":"/tmp/opencode/card-copy.tar.gz","sha256":"65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025"}
```

摘要错误以退出码 5 中止：

```console
$ gapit db install src/gapit/data/snapshots/card.tar.gz \
    --sha256 0000000000000000000000000000000000000000000000000000000000000000 \
    --output /tmp/opencode/card-bad.tar.gz
{"schema":"gapit.error/1","code":"CHECKSUM_MISMATCH","message":"SHA256 mismatch for src/gapit/data/snapshots/card.tar.gz","context":{"source":"src/gapit/data/snapshots/card.tar.gz","expected":"0000000000000000000000000000000000000000000000000000000000000000","actual":"65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025"}}
```

## gapit 构建的数据库里有什么

```
<datadir>/<name>/
  records.jsonl         truth source: one Record JSON object per line
  sequences             generated FASTA projection (gapit/v1 tagged headers)
  sequences.n*|p*       BLAST index built from sequences
  gapit-manifest.json   provenance sidecar, written last
```

### records.jsonl

这是可编辑的事实来源。每个筛查产物都由它生成，所以编辑记录后用 `--force` 重新抓
取会再生一切下游产物。每行一个 JSON 对象；card 的第一条记录，原样如下：

```json
{"db":"card","gene":"23S_rRNA_(adenine(2058)-N(6))-methyltransferase_Erm(A)","sequence":"ATGAAACAGAAAAACCCGAAAAATACGCAAAATTTCATTACATCTAAAAAGCATGTAAAGGAAATATTAAAATATACGAATATCAATAAACAAGATAAAATAATAGAAATTGGGTCAGGAAAAGGACATTTTACCAAGGAACTTGTGGAAATGAGTCAACGGGTGAATGCTATAGAGATTGATGAAGGTTTATGTCATGCCACGAAAAAAGCAGTTGAACCTTTTCAGAATATAAAAGTTATTCATGAGGATATTTTGAAGTTTAGCTTTCCTAAAAATACAGACTATAAAATATTTGGTAATATTCCCTACAATATTAGTACTGATATTGTAAAAAAGATTGCTTTTGATAGTCAAGCGAAATATAGCTACCTTATTGTAGAGAGGGGATTTGCTAAAAGGTTGCAAAATACCCAACGAGCTTTAGGTTTGCTGTTAATGGTGGAAATGGATATAAAAATTCTTAAAAAAGTGCCACGAGCATATTTTCACCCTAAGCCTAATGTAGATTCTGTATTGATTGTACTTGAAAGGCATAAACCATTTATTTTAAAGAAGGACTACAAAAAGTATAGATTTTTCGTTTATAAATGGGTAAACAGGGAATATCATGTTCTTTTTACTAAAAATCAATTAAGACAGGTGCTGAAGCATGCGAATGTTACTGATCTTGATAAATTATCCAATGAACAATTTTTGTCTGTTTTCAATAGTTACAAATTATTTCAATAA","accession":"AF002716.1:210-942","function":["lincosamide","macrolide","streptogramin"],"product":"Variant of ErmA (ARO:3000347) found in Streptococcus pyogenes. Confers the MLSb phenotype.","source_id":"3005099"}
```

| 字段 | 类型 | 含义 |
|---|---|---|
| `db` | string | 来源数据库名 |
| `gene` | string | 基因名，gapit 报告所依赖的标识 |
| `sequence` | string | 归一化序列（大写；歧义碱基记 `N`/`X`） |
| `accession` | string | 登录号或坐标跨度，未发表时为空 |
| `function` | array of string | 功能类别，已排序（见下方词汇表） |
| `product` | string | 产物描述，缺失时为 `n/a` |
| `source_id` | string | 提供商原生的记录标识 |

### gapit-manifest.json

manifest 为构建背书：抓了什么、何时、从哪里、用哪些工具版本。下面是一份真实的
manifest，来自 plasmidfinder 数据库：

```json
{
  "schema": "gapit.manifest/1",
  "name": "plasmidfinder",
  "source_urls": [
    "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip"
  ],
  "fetched_at": "2026-09-17T23:11:36Z",
  "sha256": "a18dc9dab762fe1e23a90c11314ee179c37c59dad833a5f465a78dc245970c2a",
  "n_records": 488,
  "dbtype": "nucl",
  "header_format": "gapit/v1",
  "upstream_version": "",
  "tool": {
    "name": "gapit",
    "version": "0.4.0"
  },
  "makeblastdb_version": "blastn: 2.17.0+",
  "minimap2_version": "2.31-r1302"
}
```

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.manifest/1` |
| `name` | string | 数据库名 |
| `source_urls` | array of string | 内容的上游来源 URL |
| `fetched_at` | string | ISO-8601 UTC 时间戳，秒精度 |
| `sha256` | string | 构建出的 `sequences` 文件摘要，64 个小写十六进制字符 |
| `n_records` | integer | 记录数，与 `records.jsonl` 一致 |
| `dbtype` | string | `nucl` 或 `prot`，构建时显式选择 |
| `header_format` | string | gapit 构建的数据库恒为 `gapit/v1` |
| `upstream_version` | string | 源有版本标签时的上游版本名 |
| `license` | string | 数据库内容许可证，仅提供商声明时出现（kaptive：`GPL-3.0 (database content)`） |
| `note` | string | 自由文本来源说明（kaptive：源文件 + 引用），设置时才出现 |
| `tool` | object | 构建数据库的 gapit 的 `{name, version}` |
| `makeblastdb_version` | string | 构建索引的 BLAST+ 版本 |
| `minimap2_version` | string | 构建环境中的 minimap2 版本（reads 模式在内存中建索引，不生成 `.mmi`） |

`records.jsonl` 和 manifest 是文件契约。它们从不出现在 stdout，也未注册到
`gapit schema`。

## 表头格式与兼容性

gapit 读取两种表头格式，逐记录检测：

- **abricate 遗留格式**：`>DB~~~GENE~~~ACCESSION~~~RESISTANCE PRODUCT`。字段可以缺
  失；缺失的尾部字段解码为空字符串。
- **gapit 原生 `gapit/v1`**：`>gapit|db=<v>|gene=<v>|acc=<v>|func=<v> PRODUCT`，值
  经百分号编码，基因名和产物因此可以逐字节穿越 BLAST 和 minimap2，空格和竖线也在
  内。

两者解码成相同的四个身份字段，TSV/JSON 输出形状也完全一致。原生数据库的
`RESISTANCE` 列承载 `func` 值。

> **警告：** abricate 读不了 gapit 构建的数据库。`gapit/v1` 表头格式不是 `~~~` 约
> 定。gapit 读 abricate 数据目录是刻意设计的单向通道；下游还有工具要跑 abricate
> 时，请保留一份 legacy 副本。

### `func` 词汇表

每个原生提供商从固定词汇表中填 `function`，这些值一路流向 TSV 的 `RESISTANCE` 列
和 JSON 的 `resistance` 字段：

| 提供商 | `func` 值 |
|---|---|
| `ncbi`、`resfinder`、`argannot`、`card` | 源的抗生素类别 |
| `megares` | MEGARes 的 class 字段 |
| `ecoh` | `H-antigen`、`O-antigen`，或由等位基因派生的抗原回退值 |
| `vfdb`、`ecoli_vf`、`victors`、`upec_expec_vf` | `virulence` |
| `plasmidfinder` | `replicon` |
| `bacmet2` | `biocide` |

## 带上自己的数据库

`gapit db build NAME FASTA` 把任意参考基因 FASTA 一条命令变成完整构建的 gapit 原生
数据库：`records.jsonl`、带 `gapit/v1` 表头的 `sequences` 投影、BLAST 索引和
manifest。下面是一个两基因的合成 FASTA 加一个元数据 TSV（字段见下），在临时数据目
录上运行：

```console
$ gapit db build tinyamr my_genes.fa --datadir ./db --tsv my_meta.tsv
gapit: generated /tmp/opencode/gapit-build-demo/db/tinyamr/sequences
gapit: self-check passed for tinyamr
gapit: BLAST index built (nucl)
{"db":"tinyamr","records":2,"dbtype":"nucl","destination":"/tmp/opencode/gapit-build-demo/db/tinyamr"}
```

stderr 承载逐步进度，stdout 是一行 JSON 回执（字段与 `db fetch` 相同，见上面的
表）。数据库立即可筛查：

```console
$ gapit screen --datadir ./db --db tinyamr contig.fa
Processing: contig.fa
Found 1 genes in contig.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
contig.fa	contig1	1	240	+	syn_betalac	1-240/240	===============	0/0	100.00	100.00	tinyamr	SYN-0001	synthetic class A beta-lactamase	ampicillin;cephalosporin
```

`ACCESSION` 和 `RESISTANCE` 来自 TSV 合并；`PRODUCT` 来自 FASTA 描述。
`--dbtype nucl|prot` 强制分子类型；默认由 abricate 的分子类型启发式从序列本身判
断。输入可以是 plain、`.gz` 或 `.bz2`。

### 表头检测

表头类型逐记录检测，混合文件也能工作：

| FASTA 表头 | Gene | Accession | Function | Product |
|---|---|---|---|---|
| `>syn_betalac synthetic class A beta-lactamase` | `syn_betalac` | — | — | `synthetic class A beta-lactamase` |
| `>olddb~~~sul1~~~U12338.4:1-940~~~SULFONAMIDE sulfonamide resistance` | `sul1` | `U12338.4:1-940` | `SULFONAMIDE` | `sulfonamide resistance` |
| `>gapit\|db=old\|gene=sul2\|acc=X\|func=streptomycin aminoglycoside` | `sul2` | `X` | `streptomycin` | `aminoglycoside` |

没有任何描述文本的普通表头回退到 `--description TEXT`，再回退到基因名。每条记录的
`db` 字段始终是 NAME（正在构建的数据库），`source_id` 原样保留原始 id 标记。畸形的
`gapit|` 表头使构建失败并报 `HEADER_MALFORMED`（退出码 4）；结构非法的 FASTA 报
`INVALID_FASTA`（退出码 5）。

### 从 TSV 合并元数据

`--tsv FILE` 为每个基因增补或覆盖 accession 和功能类别。规则：

- 表头行必须存在且必须含 `gene` 列；`accession` 和 `function` 按文件可选（缺列就
  是不合并），多余的列被忽略。缺 `gene` 列以退出码 5 `METADATA_MALFORMED` 失败。
- 行按基因取键：重复时**第一行胜出**，stderr 警告（`--quiet` 静默）。
- `function` 列多类别时用 `;` 分隔，例如 `ampicillin;cephalosporin`。
- 只出现在 TSV 里的基因产生 stderr 警告并被跳过：FASTA 是存在与否的唯一事实。

上面示例用到的 TSV，全文：

```console
$ cat my_meta.tsv
gene	accession	function
syn_betalac	SYN-0001	ampicillin;cephalosporin
```

### 重建

与 `db fetch` 一样，已有的数据库绝不被静默覆盖：

```console
$ gapit db build tinyamr my_genes.fa --datadir ./db
{"schema":"gapit.error/1","code":"DB_ALREADY_EXISTS","message":"won't overwrite existing database tinyamr (use --force)","context":{"db":"tinyamr"}}
$ gapit db build tinyamr my_genes.fa --datadir ./db --tsv my_meta.tsv --force
gapit: generated /tmp/opencode/gapit-build-demo/db/tinyamr/sequences
gapit: self-check passed for tinyamr
gapit: BLAST index built (nucl)
{"db":"tinyamr","records":2,"dbtype":"nucl","destination":"/tmp/opencode/gapit-build-demo/db/tinyamr"}
```

`records.jsonl` 是可编辑的事实来源：编辑它（或 FASTA）后用 `--force` 重建，每个下
游产物都会再生。

### 手动（abricate 风格）路径

推荐走 `db build`，但任何 abricate 格式的数据目录也能原样使用：把 `--datadir` 指向
含 `<name>/sequences`（带 `~~~` 表头的核苷酸 FASTA）的目录，跑一次 `gapit setupdb`
建 BLAST 索引。分子类型启发式沿用 abricate，决定 nucl 还是 prot。

用一个复制进临时数据目录的三基因微型数据库端到端演示：

```console
$ ls /tmp/opencode/gapit-own-db/db/tinyamr
sequences
$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr /tmp/opencode/gapit-own-db/gap.fa
Processing: /tmp/opencode/gapit-own-db/gap.fa
{"schema":"gapit.error/1","code":"DATABASE_NOT_INDEXED","message":"Database /tmp/opencode/gapit-own-db/db/tinyamr/sequences is not indexed, please try: gapit setupdb","context":{"db":"/tmp/opencode/gapit-own-db/db/tinyamr/sequences"}}
$ gapit setupdb --datadir /tmp/opencode/gapit-own-db/db
Indexed tinyamr (3 sequences, nucl)
$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr /tmp/opencode/gapit-own-db/gap.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
/tmp/opencode/gapit-own-db/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
```

`gapit setupdb` 为每个含有可读 `sequences` 文件的子目录建索引，当它们全部具备
`sequences.nin` 或 `sequences.pin` 时退出 0。

### reads 模式与遗留数据目录

reads 筛查（见 [reads.md](./reads.md)）中，minimap2 总是在内存中为 `sequences`
FASTA 建索引；不生成 `.mmi`，旧版 gapit 留在数据目录里的 `.mmi` 也绝不会被使用。
持久化索引尝试过并被否决：默认构建的 `.mmi` 会覆盖 `-x sr` 预设的索引参数
（minimap2 警告 "-k, -w or -H overridden by prebuilt index"），在决定性的基准里错
分了近缘同源基因（CTX-M/SHV 等位基因分歧），且在当前数据库规模下比内存索引更慢。
遗留的 abricate 数据目录走同一条内存路径。

从表头格式、蛋白质数据库、元数据 TSV 到故障排查表的完整 `db build` 实战演练见
[自定义数据库](./custom-db.md)。

## 基因簇数据库（GBK/GFF）

`gapit db build NAME input.gbk`（或 `.gbff`、`.gb`、`.gff`、`.gff3`，均支持
plain 或 `.gz`/`.bz2`）构建另一种数据库：一条记录 = 一个基因位点。输入类型按后缀
检测（`--kind gene|cluster` 必须一致，否则构建报用法错误）；上一节的 FASTA 基因流
水线不受影响。

```console
$ gapit db build mycps cps_loci.gbk --datadir ./db
gapit: generated ./db/mycps/sequences
gapit: self-check passed for mycps
gapit: BLAST index built (nucl)
{"db":"mycps","records":2,"dbtype":"nucl","destination":"./db/mycps"}
```

解析（不用 Biopython）：GenBank 记录变成位点（`LOCUS` 名 = 位点 id），其 `CDS` 特
征成为基因；source 特征的 note 承载标签和类型（kaptive 的 `K locus:`/`K type:`
风格，回退到 Bakta 风格再回退到位点 id）。GFF3 输入的序列来自内嵌的 `##FASTA` 块
或 `<stem>.fa/.fna/.fasta` 边车文件。复合 `join()` CDS 定位以类型化错误拒绝：请拆
分记录或使用简单 CDS 注释。

产物与基因数据库不同：`sequences` 每个位点一条裸位点 id 记录（普通表头，abricate
读不了）；`features.json` 是 `gapit.features/1` 特征表（用 `gapit schema features`
自省）；manifest 声明 `kind: cluster`；没有 `records.jsonl`。筛查分发到 minimap2
基因簇引擎（[screen.md](./screen.md)）。

### typing.json：声明式表型评分

`--typing FILE` 把一份经过校验的 `gapit.typing/1` 文档安装进数据库成为
`typing.json`；筛查随后为每个最佳判定注释表型和可解释的分数分解（TSV 的 PHENOTYPE
列和 gapit.cluster/1 的 `phenotype_detail`，见
[outputs.md](./outputs.md)）。一个最小的 `weighted_genes` 示例：

```json
{
  "schema": "gapit.typing/1",
  "rules": [
    {
      "model": "weighted_genes",
      "phenotype": "K1",
      "weights": {"wzx": 1.0, "wzy": 3.0},
      "negative": {"rfaD": -1.0},
      "identity_floor": 95.0,
      "require_any": ["wzx", "wzy"]
    }
  ],
  "cutoff": 0.9,
  "ambiguity_margin": 0.05,
  "fallback": "unknown"
}
```

三种规则类型：`weighted_genes`（基因存在性，带一致性下限、负标记和 any-of 门）、
`cluster_match`（位点覆盖度/一致性/关键基因的加权分量，各分量带下限）、
`learned_linear`（在命名特征 `gene:<id>:<cov|ident|present>` /
`cluster:<locus>:<coverage|identity>` 上训练的 sigmoid），外加决策层应用的文档级
`cutoff`、`ambiguity_margin` 和 `fallback`。完整 schema 用 `gapit schema typing`
打印。构建时校验是严格的：畸形文档以 `TYPING_MALFORMED` 失败；规则引用了输入没有
的基因或位点时，在任何产物写出之前以 `TYPING_UNKNOWN_GENE` 失败（同一检查在筛查时
重跑，所以手改过的数据库藏不住死引用）。

对照带标签的 assembly 调校 typing 文档正是 `scripts/cluster_calibration.py` 的用途
（开发者工具）：它筛查每个带标签的样本，打印逐期望表型的分数分布、判定×期望的一致
矩阵和分歧清单，即未来 `learned_linear` 训练流程要消费的那个循环：

```console
$ python scripts/cluster_calibration.py mycps --datadir ./db labels.tsv
```
