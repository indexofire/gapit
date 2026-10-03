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

`gapit db list` 列出所有已知数据库及其安装状态：

```console
$ gapit db list
NAME	PROVIDER	STATUS	DBTYPE	DESCRIPTION
ab_k	Kaptive (Kenyon lab)	available	nucl	A. baumannii K locus (Kaptive)
ab_o	Kaptive (Kenyon lab)	available	nucl	A. baumannii OC locus — official keyword ab_o (Kaptive)
argannot	IHU Méditerranée-Infection	available	nucl	ARG-ANNOT acquired resistance genes
bacmet2	University of Gothenburg	available	prot	BacMet2 experimentally confirmed biocide/resistance genes (protein)
card	McMaster University	available	nucl	CARD protein homolog resistance models
ecoli_kps	Kaptive (Gladstone lab)	available	nucl	E. coli group 2+3 capsular polysaccharide loci (Kaptive)
ecoli_vf	PHAC-NML	available	nucl	E. coli virulence factors (phac-nml)
kosc_k	Kaptive (klebgenomics)	available	nucl	K. oxytoca species complex K locus (Kaptive)
kosc_o	Kaptive (klebgenomics)	available	nucl	K. oxytoca species complex O locus (Kaptive)
kpsc_k	Kaptive (klebgenomics)	available	nucl	K. pneumoniae species complex K locus (Kaptive)
kpsc_o	Kaptive (klebgenomics)	available	nucl	K. pneumoniae species complex O locus (Kaptive)
megares	MEG Lab	available	nucl	MEGARes antimicrobial resistance genes
plasmidfinder	DTU CGE	available	nucl	CGE PlasmidFinder replicons
vfdb	USTC (VFDB)	available	nucl	VFDB virulence factors (set A, nucleotide)
victors	University of Chicago	available	nucl	Victors virulence factors
ecoh	Holt lab (srst2)	bundled	nucl	E. coli O and H antigens (srst2 EcOH)
ecoli_dec	gapit-curated (public-domain sources)	bundled	nucl	Diarrheagenic E. coli marker panel (GB 4789.6 + risk-monitoring designation)
lm_doumith	gapit-curated (public-domain INSDC sources)	bundled	nucl	Listeria monocytogenes serogrouping (Doumith 2004)
ncbi	NCBI	bundled	nucl	NCBI AMRFinderPlus (reference finder) curated AMR
resfinder	DTU CGE	bundled	nucl	CGE ResFinder acquired resistance genes
upec_expec_vf	FordeGenomics	bundled	nucl	UPEC/ExPEC virulence genes (FordeGenomics)
```

NAME 是传给 `--db` 的数据库名；PROVIDER 是上游维护机构。当
`<datadir>/<name>/gapit-manifest.json` 存在时 STATUS 显示 `installed (N)`，N 是
记录数，否则显示 `available`。最后六行是按名称排序的 wheel 内置数据库（见下文
"内置数据库（即装即用）"一节）：物化前为 `bundled`，物化后为 `installed (N)`。
其中四个（`ecoh`、`ncbi`、`resfinder`、`upec_expec_vf`）同时是注册表提供商 ——
这类名字只渲染一次，位于内置区（绝不会重复出现一条注册表 `available` 行），
`gapit db fetch <name>` 仍是它们的上游刷新路径。目录之外装进数据目录的数据库
同样会列出：排在内置区之后、按名称排序 —— `db build` 产物（基因与 cluster
两类都是）以及无 manifest 的目录（abricate 风格或 `db install` 落盘的字节）以
PROVIDER `local` 渲染，记录数取自 manifest（无 manifest 时从 FASTA 计数），
DBTYPE 依次取自 manifest、BLAST 索引后缀、abricate 的字母启发式。
在交互式终端上，同样的行会渲染为
带样式的 rich 表格；通过管道或重定向输出时始终保持上面的纯 TSV。面向 agent，
`--json` 输出 `gapit.dblist/1` 文档（前三个条目，输出有删节）：

```console
$ gapit db list --json
{
  "schema": "gapit.dblist/1",
  "providers": [
    {
      "name": "argannot",
      "vendor": "IHU Méditerranée-Infection",
      "description": "ARG-ANNOT acquired resistance genes",
      "dbtype": "nucl",
      "installed": true,
      "records": 2224
    },
    {
      "name": "bacmet2",
      "vendor": "University of Gothenburg",
      "description": "BacMet2 experimentally confirmed biocide/resistance genes (protein)",
      "dbtype": "prot",
      "installed": true,
      "records": 746
    },
    {
      "name": "card",
      "vendor": "McMaster University",
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
| `providers` | array | 每个已知数据库一项 |
| `providers[].name` | string | 数据库名，即传给 `--db` / `gapit db fetch` 的名字 |
| `providers[].vendor` | string | 上游维护机构（NCBI、DTU CGE、Kaptive (klebgenomics) 等） |
| `providers[].description` | string | 内容简述 |
| `providers[].dbtype` | string | `nucl`（用 blastn 筛查）或 `prot`（用 blastx 筛查） |
| `providers[].installed` | boolean | 数据目录里存在 manifest 时为 true |
| `providers[].records` | integer | 记录数，数据库未安装时省略 |
| `providers[].license` | string | 上游内容许可证，仅提供商声明时出现（card、vfdb、ecoli_vf、kaptive） |
| `providers[].source` | string | 目录之外装进数据目录的数据库（`db build` / `db install`）为 `local`，wheel 内置数据库为 `bundled`（见下）；注册表条目省略该字段 |

`gapit db list` 是唯一的列表入口：即上面的数据库目录，`--json` 返回
`gapit.dblist/1` 文档。（早期独立的列表命令及其 schema 已移除；这个入口刻意放弃了
abricate `--list` 的逐字节一致性。）见 [outputs.md](./outputs.md)。

## 内置数据库（即装即用）

大多数数据库在 `db fetch` 时才从上游下载，因为许可证禁止再分发（见"提供商"一节）。
少数**内置**数据库随 gapit wheel 一起发布 —— 共六个，其内容经过逐条记录的来源
审计（公共领域来源，或 Apache-2.0 / BSD-3-Clause / MIT 宽松许可面板；GPL 与
非商业许可内容绝不随 wheel 分发）。无需下载、无需确认：内置数据库开箱即用。

| 名称 | 内容 | 快照日期 | 分型 |
|---|---|---|---|
| `ecoh` | 大肠杆菌 O/H 抗原基因（597 条记录，srst2 EcOH） | 2026-10-02 | — |
| `ecoli_dec` | 致腹泻大肠杆菌标志基因面板（17 条记录，GB 4789.6 + 风险监测判定） | 2026-10-02 | `gapit.typing/2`（`gb4789_6`、`risk_monitoring` 双方案） |
| `lm_doumith` | 单核细胞增生李斯特菌血清群分型（5 条记录，Doumith 2004 标记） | 2026-10-03 | `gapit.typing/2`（方案 `doumith_serogroup`） |
| `ncbi` | NCBI AMRFinderPlus 精选 AMR（8373 条记录） | 2026-10-02 | — |
| `resfinder` | CGE ResFinder 获得性耐药基因（3206 条记录） | 2026-10-02 | — |
| `upec_expec_vf` | UPEC/ExPEC 毒力基因（77 条记录，FordeGenomics） | 2026-10-02 | — |

**许可证来源声明。** `ecoli_dec` 的每条记录都取自公共领域的一级提交（NCBI
RefSeq/GenBank/DDBJ 收录号，与权利方面板等位相同）；2026-10 审计中发现原面板里
混入的两条 VFDB 来源记录已被移除并替换。`lm_doumith` 面板以同样方式取自公共领域
的 INSDC 记录（收录号见每条头部；反向互补标记已做链向校正，`lmo0737` 经与 EGD-e
比对验证 100% 一致）。四个提供商快照（`ecoh`、`ncbi`、`resfinder`、
`upec_expec_vf`）通过了 2026-10 对**全部**记录的内容级审计：
NCBI AMRFinderPlus 内容属公共领域（美国政府作品），CGE ResFinder 数据库为
Apache-2.0，srst2 的 EcOH 为 BSD-3-Clause，FordeGenomics 的 UPEC-ExPEC 面板为
MIT —— 四者均允许随 MIT 许可的 wheel 再分发。回归测试锁定了六个内置库的头部
不出现任何 `VF*` / `VFDB` / `ARO:` 标签。

**快照是一个时间点。** 每个快照的 `sequences` 文件与快照当日 `gapit db fetch
<name>` 从上游产出的结果逐字节一致（resfinder 由回归测试对着一份固定的上游
归档副本逐字节验证）。随着上游持续更新，wheel 副本会逐渐滞后：`gapit db fetch
<name>`（单独使用，或加 `--force` 覆盖已物化副本）会重新下载**最新**上游内容
并重建数据目录里的库 —— 这是刷新路径。`bundled.json` 记录快照日期备查。
（v1 的 `db outdated` 仍以 manifest 为准：报告的是已装内容的物化/抓取时间，
不做内置-vs-上游比较 —— 将 wheel 快照与最新上游对比属未来工作。）

**物化机制。** 第一次 `gapit screen ... --db <内置库名>` 发现数据目录里没有该库
时，会自动构建它 —— stderr 一行提示（`gapit: materializing bundled database
ecoli_dec (17 records) into <datadir>`，`--quiet` 可静默），随后走标准基因库
构建管线：`records.jsonl`、gapit/v1 头部的 `sequences`、BLAST 索引、（库自带时）
`typing.json` 副本，以及盖上 `source: "bundled"` 的 manifest。全程确定性、零
网络。仅这条路径会在数据目录缺失时自动创建它；`gapit setupdb` 在建索引的同时
物化全部内置库；manifest 一旦存在，两者重跑都是无操作。`db list` 在物化前显示
`bundled`、物化后显示 `installed (17)` —— 四个列表入口（TSV、rich 表格、
`--json`、MCP `db_list`）一致。

### 李斯特菌血清群分型（`lm_doumith`）

用于单核细胞增生李斯特菌（Listeria monocytogenes）血清群预测的分型内置库：
五个标记 —— `prs`、`lmo0737`、`lmo1118`、`ORF2819`、`ORF2110` —— 以 95/95
的序列一致性与覆盖度阈值筛检，由 `doumith_serogroup` 方案折算为 Doumith 多重
PCR 分群（Doumith 等 2004，J Clin Microbiol 42:3819），并扩展了 Huang 2011 的
4b 变体：`prs` 作为属级控制基因把关整个方案，`IIa`（1/2a 或 3a）、`IIc`
（1/2c 或 3c）、`IIb`（1/2b、3b 或 7）、`IVb`（4b、4d 或 4e）与 `IVb-v`（携带
lmo0737 的 4b 变体）—— 与 LisSero 等工具的 in-silico 语义一致。判定走两阶段
管线：`gapit screen -o table.tsv --db lm_doumith`，然后 `gapit typing
table.tsv`。

已知的局限（相关处以 notes 形式随判定输出）：

- **4b/4d/4e 在基因层面无法区分** —— 三个血清型共享所有已发表的分子标记；
  需用 cgMLST 或血清学进一步判定。
- **IVb-v 是新出现的关注克隆** —— 4b 谱面中检出 lmo0737（Huang 2011
  "unusual 4b"；ST382/ST554 克隆）时，该规则先于 IVb 声明，变体优先胜出。
- **EGD-e 判为 IIc** —— 1/2a 参考基因组携带 lmo1118 而落入 IIc；这是 Doumith
  2004 原文记载的现象，不是工具假象。
- **水平基因转移是双向陷阱** —— 谱系 II 的 lmo0737 盒会在谱系间移动，罕见的
  ORF 标记水平转移也可能造成貌似 IIb 的假谱面；在流行病学上不太可能的语境里，
  单标记 IIb 判定需谨慎对待。
- **`prs` 是属级基因，并非 Lm 特异** —— 其他李斯特菌物种也携带它；prs 阳性但
  无任何血清群标记的分离株回退为 fallback，而不是被判为"非李斯特菌"。

回退判定 —— `untypeable (4a/4c, atypical profile, or non-Lm Listeria)` ——
是这五个标记无法给出答案时的诚实输出。

## 检查数据库新鲜度

`gapit db outdated` 报告每个已安装数据库相对于过期阈值的年龄。对装满的数据
目录（输出裁剪到十二行中的三行）：

```console
$ gapit db outdated --days 30
NAME	FETCHED_AT	AGE_DAYS	STATUS
argannot	2026-09-17T23:22:06Z	2.58	ok
bacmet2	2026-09-17T23:13:51Z	2.58	ok
card	2026-09-17T23:14:25Z	2.58	ok
...
```

安装已久的数据库只会报告 `stale`：

```console
$ gapit db outdated --datadir /tmp/opencode/gapit-outdated-demo
NAME	FETCHED_AT	AGE_DAYS	STATUS
card	2020-01-01T00:00:00Z	2454.55	stale
```

| 状态 | 含义 |
|---|---|
| `ok` | 在阈值之内 |
| `stale` | `age_days` 超过 `--days`（默认 90；`--days 0` 把一切都标为 stale） |

过期只是一份报告，绝不是错误状态：无论多旧命令都退出 0。退出码 4
（`DATADIR_NOT_FOUND`、`DATADIR_EMPTY`）对应缺失的数据目录或没有任何已安装数据库的
目录；无法解析的 `fetched_at` 是 `MANIFEST_MALFORMED`（退出码 5）。内置数据库的
报告与其他库一样以 manifest 为准：已物化的内置库从物化时刻起算年龄，v1 不做
内置-vs-上游比较（见"内置数据库（即装即用）"一节 —— 属未来工作）。面向 agent，
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

gapit 自带十九个提供商。每一个——包括 `card` 和 `vfdb`——都在抓取时从上游下载：
除上文六个通过审计、宽松许可的内置库之外，包里不内置任何内容，因为多家上游许可
（CARD 的 McMaster 条款、VFDB 的 CC BY-NC、
Kaptive 的 GPL-3.0）禁止随 MIT 许可的发行版再分发。七个 kaptive 提供商是**基因簇**
数据库（抓取时经基因簇管线构建，见下文 Kaptive 提供商一节）。

| 名称 | 维护方 | 内容 | dbtype |
|---|---|---|---|
| `ncbi` | NCBI | NCBI AMRFinderPlus（reference finder）精选 AMR | nucl |
| `card` | McMaster University | CARD protein homolog resistance models | nucl |
| `resfinder` | DTU CGE | CGE ResFinder acquired resistance genes | nucl |
| `argannot` | IHU Méditerranée-Infection | ARG-ANNOT acquired resistance genes | nucl |
| `plasmidfinder` | DTU CGE | CGE PlasmidFinder replicons | nucl |
| `megares` | MEG Lab | MEGARes antimicrobial resistance genes | nucl |
| `ecoh` | Holt lab (srst2) | E. coli O and H antigens (srst2 EcOH) | nucl |
| `vfdb` | USTC (VFDB) | VFDB virulence factors (set A, nucleotide) | nucl |
| `ecoli_vf` | PHAC-NML | E. coli virulence factors (phac-nml) | nucl |
| `bacmet2` | University of Gothenburg | BacMet2 experimentally confirmed biocide/resistance genes (protein) | prot |
| `victors` | University of Chicago | Victors virulence factors | nucl |
| `upec_expec_vf` | FordeGenomics | UPEC/ExPEC virulence genes (FordeGenomics) | nucl |
| `kpsc_k` | Kaptive (klebgenomics) | K. pneumoniae species complex K locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `kpsc_o` | Kaptive (klebgenomics) | K. pneumoniae species complex O locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `kosc_k` | Kaptive (klebgenomics) | K. oxytoca species complex K locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `kosc_o` | Kaptive (klebgenomics) | K. oxytoca species complex O locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `ab_k` | Kaptive (Kenyon lab) | A. baumannii K locus（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `ab_o` | Kaptive (Kenyon lab) | A. baumannii OC locus — 官方关键字 `ab_o`（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |
| `ecoli_kps` | Kaptive (Gladstone lab) | E. coli group 2+3 capsular polysaccharide loci（Kaptive；cluster，GPL-3.0，抓取时下载） | nucl |

蛋白质数据库（`bacmet2`）经 `blastx` 筛查；核苷酸数据库经 `blastn`。基因簇数据库
经 minimap2 基因簇引擎筛查（[screen.md](./screen.md)）。

### Kaptive 提供商（GPL，抓取时下载）

七个 kaptive 提供商封装了 [Kaptive](https://github.com/klebgenomics/Kaptive) 的参考
数据库（Wyres et al., J Clin Microbiol 2020：使用这些数据库的结果请引用
Kaptive）。数据库名（NAME）就是 Kaptive v3 数据库文档中的**官方安装关键字**
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

### 默认集合

`gapit db fetch all` 安装默认集合：先 `card` 后 `vfdb`。与其他所有提供商一
样，它从上游下载、变换记录、并在本地构建 `sequences` 和 BLAST 索引——因此索引永远
与你机器上实际安装的 BLAST 版本匹配。需要网络。不带参数的 `gapit db fetch`
只打印命令帮助：多库下载必须显式指定（`all` 或单个 NAME），绝不因省略参数而隐式
触发。

### 抓取指定数据库

`gapit db fetch <name>` 对单个数据库跑完整流水线。下面是安装 `card`
（stderr 进度行，随后 stdout 一行 JSON 回执）：

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo card
gapit: downloaded 1 source file(s)
gapit: read 6059 records from card
gapit: generated /tmp/opencode/gapit-dbs-demo/card/sequences
gapit: self-check passed for card
gapit: BLAST index built (nucl)
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/opencode/gapit-dbs-demo/card"}
```

| 回执字段 | 类型 | 含义 |
|---|---|---|
| `db` | string | 数据库名（NAME） |
| `records` | integer | 写入 `records.jsonl` 的记录数 |
| `dbtype` | string | `nucl` 或 `prot` |
| `destination` | string | 安装后的数据库目录 |

每个数据库在抓取时从其上游来源下载，因此需要网络。未知名字在解析数据目录之
前就会被拒绝：

```console
$ gapit db fetch nosuchdb
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"unknown database: nosuchdb (available: argannot, bacmet2, card, ecoh, ecoli_vf, megares, ncbi, plasmidfinder, resfinder, upec_expec_vf, vfdb, victors)","context":{"db":"nosuchdb"}}
```

### 重新抓取

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

## 安装本地文件

`gapit db install SOURCE --sha256 HASH --output TARGET` 是带校验的本地文件安装。它
对提供商目录、档案、序列内容一无所知：把 SOURCE 流式过 SHA256，与 `--sha256` 比对摘
要，通过后才原子替换 TARGET。校验失败时已有的 TARGET 不受影响。

```console
$ sha256sum card.json
65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025  card.json
$ gapit db install card.json \
    --sha256 65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025 \
    --output /tmp/opencode/card-copy.json
{"destination":"/tmp/opencode/card-copy.json","sha256":"65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025"}
```

摘要错误以退出码 5 中止：

```console
$ gapit db install card.json \
    --sha256 0000000000000000000000000000000000000000000000000000000000000000 \
    --output /tmp/opencode/card-bad.json
{"schema":"gapit.error/1","code":"CHECKSUM_MISMATCH","message":"SHA256 mismatch for card.json","context":{"source":"card.json","expected":"0000000000000000000000000000000000000000000000000000000000000000","actual":"65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025"}}
```

## gapit 构建的数据库里有什么

```
<datadir>/<name>/
  records.jsonl         truth source: one Record JSON object per line
  sequences             generated FASTA projection (gapit/v1 tagged headers)
  sequences.n*|p*       BLAST index built from sequences
  gapit-manifest.json   provenance sidecar, written last
  typing.json           可选的 gapit.typing/1 或 /2 评分规范（带分型的数据库）
  floors.json           可选的 gapit.floors/1 按基因一致性下限（reads 模式）
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
    "version": "0.5.1"
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
| `license` | string | 数据库内容许可证，仅提供商声明时出现（kaptive：`GPL-3.0 (database content)`；card：McMaster 非商业条款；vfdb：`CC BY-NC 4.0`） |
| `note` | string | 自由文本来源说明（kaptive：源文件 + 引用），设置时才出现 |
| `tool` | object | 构建数据库的 gapit 的 `{name, version}` |
| `makeblastdb_version` | string | 构建索引的 BLAST+ 版本 |
| `minimap2_version` | string | 构建环境中的 minimap2 版本（reads 模式在内存中建索引，不生成 `.mmi`） |

`records.jsonl` 和 manifest 是文件契约。它们从不出现在 stdout，也未注册到
`gapit schema`。

### floors.json — 按基因设置一致性下限

一个可选边车，为 reads 模式的存在判定声明逐基因的最小比对一致性（%）（schema 为
`gapit.floors/1`，用 `gapit schema floors` 内省；筛查语义与 `pic` 实战示例见
[reads.zh.md](./reads.zh.md#按基因设置一致性下限gapitfloors1数据库侧)）：

```json
{"schema": "gapit.floors/1", "default": null, "genes": {"pic": 90.0}}
```

| 字段 | 类型 | 含义 |
|---|---|---|
| `schema` | string | 恒为 `gapit.floors/1` |
| `default` | number 或 null | 未在 `genes` 中列出的基因的下限（`null` = 无下限） |
| `genes` | object | `{基因: 最小一致性 %}`，取值在 [0, 100] 内 |

文件的存在本身就是开关——没有 manifest 字段，也没有 CLI 开关。不带边车的数据库
输出与之前逐字节一致（下限是数据库驱动的可选开启，不是 CLI 契约变更）；blastn
contig 路径从不读取它们。通过 `db build --floors` 安装（见下节），或把文件直接放
进既有数据库目录；筛查时同样会校验（内容损坏以 `FLOORS_MALFORMED` 失败，exit 4）。
内置的 `ecoli_dec` 携带 `{"genes": {"pic": 90.0}}`——SPATE 同源假阳性的修复。

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

### 按基因一致性下限（`--floors`）

`--floors FILE` 把一份 `gapit.floors/1` 文档安装为数据库的 `floors.json` 边车——
reads 模式存在判定的逐基因最小比对一致性（见
[reads.zh.md](./reads.zh.md#按基因设置一致性下限gapitfloors1数据库侧)）。文档在
任何产物落盘前先校验：`genes` 里的每个基因必须存在于 FASTA
（`FLOORS_UNKNOWN_GENE`，exit 4——给缺失基因设下限是静默失效的安全配置），取值必
须在 [0, 100] 内、结构必须可解析（`FLOORS_MALFORMED`，exit 4）。随后文件原样拷入；
文件的存在本身就是开关（manifest 字段不变）。仅限基因数据库——GBK/GFF 簇输入对
`--floors` 报用法错误。

```console
$ cat my_floors.json
{"schema": "gapit.floors/1", "default": null, "genes": {"syn_betalac": 95.0}}
$ gapit db build tinyamr my_genes.fa --datadir ./db --tsv my_meta.tsv --floors my_floors.json
{"db":"tinyamr","records":2,"dbtype":"nucl","destination":"/tmp/opencode/gapit-build-demo/db/tinyamr"}
$ cat ./db/tinyamr/floors.json
{"schema": "gapit.floors/1", "default": null, "genes": {"syn_betalac": 95.0}}
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

`--typing FILE` 把一份经过校验的 `gapit.typing/1` 或 `gapit.typing/2` 文档安装进数据库成为
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

`gapit.typing/2` 把同样的规则包进命名 scheme（`schemes:
[{"name": "pathotype", "rules": [...], "cutoff": ..., "ambiguity_margin": ...,
"fallback": ...}]`，名称限 `[a-z0-9_]` 且唯一）；`/1` 文档等价于一个匿名的
`default` scheme，因此既有的 typed 数据库输出不变。`--typing` 也可用于基因（FASTA）
构建 —— v1 仅支持 `weighted_genes` 规则，基因引用对照 FASTA 记录校验。typed 基因
数据库的筛查输出与 untyped 完全相同（纯基因检测）；其判定来自两阶段流水线 ——
`gapit screen -o result.tsv` 写出结果表，`gapit typing result.tsv` 按 scheme 渲染
判定，输出 `gapit.typing_result/1`
（见[分型方案](./typing.zh.md#两阶段判定工作流)）。

typing/2 第二阶段新增六个规则/scheme 原语（全部可选、纯附加；完整定义见
`gapit schema typing`）：`exact_set` 规则（Doumith/志贺氏菌式标记表 —— 所有
`requires` 基因存在且所有 `excludes` 基因缺席时得 1.0，两个下限默认 90/90；满足的
exact_set 参与的 1.0 平局按声明顺序裁决而非歧义边际）、`weighted_genes`/`exact_set`
上的可选 `coverage_floor`、scheme 级 `control_gene`（prs/ipaH 式门控：缺席时整个
scheme 输出 fallback 并附 "control gene absent" 注记）、`unique_group` +
`mixed_phenotype`（同组多于一个基因存在时调用混合表型并把成对基因列入
`ambiguous`）、`compose` scheme（`"{o_group}:{k_group}"` 由兄弟 scheme 的判定渲染 ——
fallback 字符串原样流入，歧义成分使组合结果为 null 并携带该成分的变体），以及逐规则
`notes`（获胜时原样浮现到判定上）。

六个经调研的鉴定方案——完整的 Doumith 李斯特菌表、ShigaTyper 语义的志贺菌/EIEC、带
等位基因探针的 meningotype 血清群 panel、副溶血弧菌 O/K 的 Kaptive 模式、霍乱弧菌
O1/O139 与 Ogawa/Inaba、以及带文档的 DEC 占位——已编码为经过校验的示例文档（合成
标记夹具，由测试套件自检），详见[分型方案](./typing.md)。

对照带标签的 assembly 调校 typing 文档正是 `scripts/cluster_calibration.py` 的用途
（开发者工具）：它筛查每个带标签的样本，打印逐期望表型的分数分布、判定×期望的一致
矩阵和分歧清单，即未来 `learned_linear` 训练流程要消费的那个循环：

```console
$ python scripts/cluster_calibration.py mycps --datadir ./db labels.tsv
```
