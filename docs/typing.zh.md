# 分型方案（typing）

分型文档（typing document）是数据库侧的评分规范，把筛查命中转换为表型判定——血清型、血清群、致病型。通过 `gapit db build --typing FILE` 安装后，该数据库即可对筛查过它的任何样本作出判定。本页面是 `gapit.typing/2` 的完整参考，外加一份实战手册（cookbook）：六个经调研的鉴定方案，全部编码为经过校验、自带测试的示例文档。

- 数据库侧机制（`--typing`、校验、manifest）见
  [数据库](./databases.md) 的 typing.json 一节。
- 判定输出文档（`gapit.typing_result/1`）见 [输出](./outputs.md)。
- 分型文档的 JSON Schema 用 `gapit schema typing` 打印。

## 两阶段判定工作流

基因类判定是两阶段流水线（权利方设计）：**`gapit screen` 检测基因，`gapit typing`
从结果中获得型别。**

```console
$ gapit screen contigs.fa --db ecoli_dec --output result.tsv --quiet
$ gapit typing result.tsv --quiet
FILE	SCHEME	PHENOTYPE	CONFIDENCE	SCORE	RUNNER_UP	NOTES	GENES
contigs.fa	gb4789_6	EHEC	high	1.0000	EAEC (1.0000)	...	escV;stx2a;uidA
contigs.fa	risk_monitoring	EHEC	high	1.0000	EAEC (1.0000)	...	escV;stx2a;uidA
```

第一阶段是纯基因检测：typed 基因数据库与 untyped 的筛查输出逐字节相同——每种格式
都是冻结的 15 列 abricate 表。第二阶段读取一份或多份 gapit/abricate 筛查结果表
（TSV 或 CSV），从所有行共享的 `DATABASE` 列解析数据库，把每个 `(FILE, GENE)` 折叠
为其 `(%IDENTITY, %COVERAGE)` 最佳的行（并列时首行获胜——与报告命中原有的逐基因
折叠一致），再对该 FILE 的基因评估数据库的 `typing.json`。

两个阶段也可以直接经管道衔接——规范一行式：

```console
$ gapit screen 1.fna --db ecoli_dec --quiet | gapit typing --quiet
```

只要 stdin 不是终端（管道或重定向），`gapit typing` 就从 stdin 读取结果表；在终端
下裸调用仍然打印帮助。显式 `-` 标记——`gapit typing -`——即使挂着终端也从 stdin
读取，直到 EOF；`-` 与文件参数混用是用法错误（v1 只取其一，stdin 或文件，不可兼得）。
管道输入的 JSON `source` 渲染为 `["-"]`。

用法：`gapit typing RESULT.tsv [RESULT2.tsv ...] [-D datadir] [-f tsv|json|md] [-o FILE]
[-q]`，或如上经管道输入（`-` 为显式 stdin）。默认 TSV 按 (FILE, scheme) 一行，列为
`FILE`、`SCHEME`、`PHENOTYPE`、
`CONFIDENCE`、`SCORE`、`RUNNER_UP`、`NOTES`、`GENES`——歧义判定的表型（与 runner-up）
渲染为 `-`，候选对写入 NOTES；`GENES` 以该 FILE 的 present 基因名（排序后、`;` 串接）
收尾每行，且在该 FILE 的每个 scheme 行上重复——即驱动该判定的命中，一目了然。
`json` 输出 `gapit.typing_result/1`（完整分数分解外加同样的逐文件 `genes` 列表；
`gapit schema typing_result`）；`md` 在标明数据库与来源表的 frontmatter 下渲染同样的
八列表。类型化错误：`DATABASE` 值混杂为 `DATABASE_MISMATCH`，无数据行的表为
`TYPING_NO_DATA`（要评估的数据库只能从数据行得知），无 `typing.json` 的数据库为
`TYPING_NO_SCHEME`，基因簇数据库为 `TYPING_CLUSTER_DB`——基因簇判定不经过该命令。

**基因簇例外。** 基因簇数据库的判定保持集成在 `gapit screen` 内：其判定单元是最佳
*位点* call，只有 minimap2 基因簇引擎能算——不存在可供重读的筛查表。筛查 typed 基因
簇数据库会直接注释 `gapit.cluster/1`（仅限 typed 的 TSV `PHENOTYPE` 列、JSON/MD 中的
`best.phenotype` + `best.phenotype_detail`），与之前完全一致。

## 一份文档，两个引擎

`gapit.typing/2` 文档是一组命名的 **scheme**；每个 scheme 是一组 **rule** 加上作用于
这些规则的判定阈值：

```json
{
  "schema": "gapit.typing/2",
  "schemes": [
    {
      "name": "serogroup",
      "rules": [ /* 每个候选表型一条 rule */ ],
      "cutoff": 0.9,
      "ambiguity_margin": 0.05,
      "fallback": "NG"
    }
  ]
}
```

由哪个引擎评估文档取决于数据库类型：

| 数据库类型 | 输入 | 引擎 | 判定出现的位置 |
|---|---|---|---|
| **gene（基因）** | 标记 FASTA | blastn 基因路径 | 对筛查结果运行 `gapit typing`：`gapit.typing_result/1`（JSON）与八列 TSV/MD 表 |
| **cluster（基因簇）** | GBK/GFF 位点 | minimap2 位点路径 | 集成在 `gapit screen` 内：`gapit.cluster/1` 的 `best.phenotype` + `best.phenotype_detail`，以及仅限已分型数据库的 TSV `PHENOTYPE` 列 |

按类型区分的两条硬约束：

- **基因数据库**只评估 `weighted_genes` 与 `exact_set` 规则（不存在位点可供
  `cluster_match`/`learned_linear` 使用）。
- **基因簇数据库**每个数据库只承载**一个** scheme——`gapit.cluster/1` 每个文件只有
  一个表型槽位。多 scheme 文档会在构建时报 `TYPING_MALFORMED`；应把每个 scheme 安装为
  独立数据库（见下文 O/K 模式）再在外部组合，或者在基因数据库上建模整个 panel，
  `compose` 即可原生运行。

## rule 与 scheme 原语

每条 rule 由其 `model` 判别。完整细节见 `gapit schema typing`；这里是工作参考。

### rule 模型

| 模型 | 打分 | 关键字段 |
|---|---|---|
| `weighted_genes` | 基因存在计分求和，按正权重总和归一到 [0, 1] | `weights`、`negative`（存在时反向计分）、`identity_floor`（必填）、`coverage_floor`（可选）、`require_any`（任一不满足则整条规则归零） |
| `exact_set` | 当且仅当每个 `requires` 基因都在、`requires_any` 非空时其中至少一个在、且每个 `excludes` 基因都不在时为 1.0，否则 0.0 | `requires`、`requires_any`（任一满足门控）、`excludes`、阈值默认 90/90 |
| `cluster_match` | 对覆盖最佳的位点做加权 coverage + identity + key genes 组合分，各分量带阈值 | `coverage{weight,floor}`、`identity{weight,floor}`、`key_genes{weight,genes}` |
| `learned_linear` | 命名特征的加权求和加 bias 过 sigmoid | `features`（`gene:<id>:<cov\|ident\|present>`、`cluster:<locus>:<coverage\|identity>`）、`weights`、`bias`、`trained_on` |

基因对 `weighted_genes`/`exact_set` "算作存在"，要求其最佳存活命中（基因簇路径上则是
gene call）清过该规则的阈值。`cluster_match` 规则的 `coverage + identity` 权重之和请
保持在 cutoff 以下：应当由 key genes 决定胜负，这样完整覆盖但缺 key genes 的位点无法
窃取判定。

### scheme 字段

| 字段 | 含义 |
|---|---|
| `name` | `[a-z0-9_]+`，文档内唯一；判定输出以它为键 |
| `cutoff` | 最高分规则必须达到才判定，否则走 fallback |
| `ambiguity_margin` | 第一名与第二名的分差必须达到，否则判定为 `ambiguous`（表型为 null，列出前两名） |
| `fallback` | 低置信度输出字符串（`NT`、`NG`、`OUT`、`KUT` 等） |
| `control_gene` | 门控：该基因缺失时 scheme 输出 fallback 并附 `control gene absent` 注记 |
| `unique_group` + `mixed_phenotype` | 组内多于一个基因存在时，调用 `mixed_phenotype` 并把这一对列入 `ambiguous`（混合感染信号） |
| `compose` | `"{o_group}:{k_group}"`：无规则的 scheme，在兄弟 scheme 判定后由其结果渲染——fallback 字符串原样流入，ambiguous 成分则产生携带其变体的 null 组合 |
| `notes`（规则级） | 原样呈现在获胜判定上的自由文本——警示、引文、标记 |

一个并列例外：涉及满足条件的 `exact_set` 的 1.0 并列按**声明顺序**裁决，而非
margin（标记表按构造互斥，并列意味着两行都命中，先声明者胜）。非 exact 并列保持
ambiguous。

## 实战手册（cookbook）

六份示例文档位于
[`tests/data/typing/schemes/`](https://github.com/indexofire/gapit/tree/main/tests/data/typing/schemes)，
每份配一个标记 FASTA 或 GBK 位点文件。每份文档都通过当前 schema 校验，
`tests/test_typing_schemes.py` / `tests/test_typing_schemes_cluster.py` 会构建每个夹具
数据库并断言下文的表型判定。**所有序列均为合成序列**——它们编码的是方案逻辑而非
策展内容；按各提供商的许可条款，真实策展是未来的工作。唯一的策展例外就放在它们
旁边：内置的 `lm_doumith` 面板（真实公共领域 INSDC 序列，独立测试模块），见下文
第一节。

### 李斯特菌血清群 —— 内置 lm_doumith 面板（基因库）

`lm_doumith` 是第六个内置数据库：取自公共领域 INSDC 记录的 Doumith 标记，加上
`doumith_serogroup` 方案（`gapit.typing/2`）——与 LisSero 相同的五标记多重 PCR
语义，`prs` 为属级 `control_gene`，每条规则 95/95 阈值（Doumith 等 2004，J Clin
Microbiol 42:3819）：

| 判定 | requires | excludes |
|---|---|---|
| `IVb-v` | `ORF2819`、`ORF2110`、`lmo0737` | `lmo1118` |
| `IVb`（4b、4d 或 4e） | `ORF2819`、`ORF2110` | `lmo0737`、`lmo1118` |
| `IIb`（1/2b、3b 或 7） | `ORF2819` | `lmo0737`、`lmo1118`、`ORF2110` |
| `IIc`（1/2c 或 3c） | `lmo0737`、`lmo1118` | `ORF2819`、`ORF2110` |
| `IIa`（1/2a 或 3a） | `lmo0737` | `lmo1118`、`ORF2819`、`ORF2110` |

`IVb-v`——携带谱系 II `lmo0737` 盒的 4b 变体（Huang 2011 "unusual 4b"；FDA 4bV；
ST382/ST554 新发克隆）——声明在 `IVb` **之前**，使携带 `lmo0737` 的
ORF2819+ORF2110 谱式判给变体。大写基因名（`ORF2819`、`ORF2110`）在构建、筛检与
分型全程原样传递。两阶段使用，零网络（首次筛检时自动物化）：

```console
$ gapit screen isolates.fa --db lm_doumith --output lm.tsv --quiet
$ gapit typing lm.tsv --quiet
FILE	SCHEME	PHENOTYPE	GENES	CONFIDENCE	SCORE	RUNNER_UP	NOTES
lm_doumith_iia.fa	doumith_serogroup	IIa	lmo0737;prs	high	1.0000	IVb-v (0.0000)	...
lm_doumith_ivb.fa	doumith_serogroup	IVb	ORF2110;ORF2819;prs	high	1.0000	IVb-v (0.0000)	...
```

fallback 为 `untypeable (4a/4c, atypical profile, or non-Lm Listeria)`；缺 `prs` 时
控制基因门触发（`control gene absent`）。已知局限以 notes 随判定输出：4b/4d/4e
共享所有标记、EGD-e 判为 IIc（Doumith 2004 原文记载）、ORF 标记的水平转移可能
伪造 IIb 谱式——详见数据库指南中的
[lm_doumith 一节](./databases.md)。

### 李斯特菌血清群 —— Doumith 表（基因库）

`doumith.json` + `doumith.fa`：Doumith 等 2004（J Clin Microbiol 42:3819）五标记多重
PCR 完整表——`prs` 作为 scheme 的 `control_gene`，`lmo0737`、`lmo1118`、`orf2819`、
`orf2110` 上各血清群一条 `exact_set` 规则：

| 血清群 | requires | excludes |
|---|---|---|
| `1/2a-3a` | `lmo0737` | `lmo1118`、`orf2819`、`orf2110` |
| `1/2c-3c` | `lmo0737`、`lmo1118` | `orf2819`、`orf2110` |
| `1/2b-3b-7` | `orf2819` | `lmo0737`、`lmo1118`、`orf2110` |
| `4b-4d-4e` | `orf2819`、`orf2110` | `lmo0737`、`lmo1118` |

fallback `NT`（4a/4c 与非典型谱式）。4b 规则携带警示注记，包括 4b*（IVb-v1）水平
基因转移陷阱：4b* 菌株获得了谱系 II 的 `lmo0737` 盒，不匹配任何朴素谱式，落入
`NT`——测试套件钉住了这一判定（J Clin Microbiol 2022, PMID 36472431）。

### 志贺菌 / EIEC —— ShigaTyper 语义（基因库）

`shigella.json` + `shigella.fa`：ShigaTyper（Wu 等 2019，Appl Environ Microbiol
85:e00165-19）的骨架，`ipaH_c` 为控制基因，各 wzx 变体接成一个 `unique_group`
（`mixed_phenotype: "mixed Shigella serotypes"`——ShigaTyper 的 "multiple wzx" 检查点
变成混合感染判定）。要点：

- **S. sonnei form I**（`Ss_wzx`/`Ss_wzy` 的 `exact_set`）声明在 **form II** 规则
  （`Ss_methylase` 的 `weighted_genes`）**之前**：抗原相变异的 form I 样本两条都拿
  1.0，exact 并列例外裁决为 form I。
- **S. dysenteriae 1** 要求 `Sd1_wzx` + `Sd1_rfp`，注记说明 rfp 阴性变异在完整
  panel 中需要自己的规则。
- **S. flexneri** 血清型转换是 `gtr`/`Oac` 上的精确集合（`2a`、`2b`、`3a`），外加
  `Y/novel` 兜底规则（`requires` Sf 基底、`excludes` 已列出的转换基因）——完整 panel
  下该规则会拆分为 Y 与 ShigaTyper 的 "novel serotype"。
- **EIEC** 近似 ShigaTyper 检查点 3：`exact_set` 要求 `EclacY` 并 `excludes` 豁免基因
  `Sb9_wzx`/`Sb15_wzx`；cadA 与 ipaB 质粒覆盖度细节作为未建模项写入注记。

### 脑膜炎奈瑟菌血清群 —— meningotype panel（基因库）

`meningotype.json` + `meningotype.fa`：`ctrA` 控制基因，每个血清群一条
`weighted_genes` 规则，基于 Mothershed 等 2004（J Clin Microbiol 42:320–328）的实时
PCR panel（`sacB`/A、`synD`/B、`synE`/C、`xcbB`/X、`synF`/Y，按 meningotype 工具的
命名）——W/Y 则通过 **EX7E 等位基因探针** `synG_EX7E_P` / `synG_EX7E_G` 以
`identity_floor` 99.5 判定（见下文）。fallback `NG`。

### 副溶血弧菌 O 与 K —— Kaptive 模式（基因簇库）

`vp_ok.json` + `vp_ok.gbk`：仿照 van der Graaf-van Bloois 等 2023（Microb Genom
9:mgen001007；16 个 O 与 71 个 K 血清型）的 Kaptive 数据库构建的 OAgc/CPSgc 样位点。
三个 scheme：

- `o_group` —— `cluster_match` 规则（coverage 0.55 / identity 0.25 / key genes 0.2，
  阈值 95/95），fallback `OUT`（O 不可分型）。`O3/O13` 规则编码合并标签：O3 与 O13
  位点含有相同的基因，无法按基因内容区分。
- `k_group` —— K 位点上同构，fallback `KUT`（K 不可分型）。
- `serotype` —— 记录 panel 输出形状的 `"{o_group}:{k_group}"` compose scheme。

基因簇数据库只接受一个 scheme，因此可运行的安装是同一 GBK 构建的**两个数据库**
（`vp_o`、`vp_k`）——恰是 Kaptive 分别发布 O 与 K 数据库的方式——各 scheme 文档从
手册文件中抽取。测试套件同时断言拆分后的判定（`O3/O13`、`K6`）与整份文档构建触发
单 scheme 守卫；O:K 组合本身在基因数据库上原生运行（由
`typing_serotyping` golden 锁定）。

### 霍乱弧菌 O1/O139 与 Ogawa/Inaba（基因簇库）

`cholerae.json` + `cholerae.gbk`：位点 `wbe`（O1：`wbeV`、`wbeW`、`wbeT`）与 `wbm`
（O139：`wbmV`、`wbmW`、`wbfZ`），外加不可分型的 `wbc` 位点，两个 scheme：

- `serogroup` —— `cluster_match` 的 O1 / O139，fallback `non-O1/non-O139`。O139 规则
  的注记记录 **wbfZ 陷阱**：`wb*` 盒通过保守的 `gmhD`/`rjg` 交换接合点重组，环境血清
  群可能携带 O139 样盒从而使单基因 `wbfZ` 检测误报（Sozhamannan 等 1999，Infect
  Immun 67:6215）——该规则改为以整个 `wbm` 位点为键。
- `subserotype` —— 通过 wbeT 等位基因探针区分 Ogawa/Inaba：**ogawa** = `wbeT` 上的
  `weighted_genes`，`identity_floor` 99.9（该决定因素是 wbeT 甲基转移酶中的单突变；
  Stroeher 等 1992，PNAS 89:2566；BMC Microbiol 2013, 13:173）；**inaba** = `wbeV`
  正向加 `negative` 中的 `wbeT: -2.0`（Inaba 的定义：wbeT 清不过等位基因阈值），并带
  Hikojima 注记——稀有、不稳定、无稳定参照，不可判定。wbeT 突变样本在 serogroup
  数据库上判 `O1`、同时在 subserotype 数据库上判 `inaba`
  （`cluster_vc_inaba` golden）。

### 致泻性大肠埃希菌 —— 双方案判定（基因库）

`dec.json` + `dec.fa`：DEC 判定编码为**同一个基因数据库上的两个 scheme**——typing/2
多方案的示范。两个方案都用精确语义（cutoff 1.0、margin 0.0）、`uidA` 作
`control_gene`、fallback `non-DEC`，并共享一条按严重度排序的规则阶梯
（EHEC > STEC/EPEC > ETEC > EIEC > EAEC——混合株以下一名 runner_up 呈现，按声明顺序
裁决）：

| 致病型 | requires | requires_any | excludes |
|---|---|---|---|
| `EHEC` | `escV` | `stx1a`/`stx1b`/`stx2a`/`stx2b` | |
| `STEC` | | `stx1a`/`stx1b`/`stx2a`/`stx2b` | `escV` |
| `EPEC` | `escV`、`bfpB` | | |
| `EPEC_atypical` | `escV` | | `bfpB` 与全部 stx 亚基 |
| `ETEC` | | `lt`/`sth`/`stp` | |
| `EIEC` | `invE` | | |
| `EAEC`（`gb4789_6`） | | `aggR`/`pic`/`astA` | |
| `EAEC`（`risk_monitoring`） | `aggR` | | |

- **`gb4789_6`** 编码 GB 4789.6-2016 panel 语义：EAEC 为 aggR/pic/astA **任一**存在。
- **`risk_monitoring`** 编码现行食品安全风险监测方案（最新食品安全风险监测方案）：
  aggR 为必备——仅有 pic/astA 不充分。携带 pic+astA 而无 aggR 的分离株在
  gb4789_6 下判 `EAEC`、在 risk_monitoring 下判 `non-DEC`——两个方案之间的标志性
  分歧，由 `ecoli_dec_typing` golden 锁定（两条判定同处一份
  `gapit.typing_result/1` 文档；捆绑数据库携带同一份文档）。
- `EPEC_atypical` 声明在 EHEC/STEC **之后**，stx 阳性分离株永远不会落入该规则；
  escV+/bfpB−/stx− 分离株越过典型 EPEC 规则后落到它上面。
- 夹具镜像了权利方数据库中的重复记录（真实库 pic ×2、sth ×3；夹具各带一份重复）：
  重复记录按基因名折叠——blastn 的 `-culling_limit 1` 保留每个查询跨度的最佳
  subject，引擎按名称折叠出每基因一条 call。

## 等位基因探针技巧（及其后继）

手册中的两个方案用存在/缺失原语判别**等位基因**而非基因：

- **EX7E（meningotype）**：synG 的 EX7E 基序肽在 W 血清群读 P、Y 读 G、双重 W/Y 读 S
  （meningotype 的 `menwy` 检查）。夹具把两条等位基因序列作为独立探针记录、
  `identity_floor` 99.5：只有完全匹配的等位基因能清过阈值，因为 ~120 bp 探针上单个
  SNP 就会跌到 99.5% 以下。
- **wbeT（cholerae）**：Ogawa 与 Inaba 之差是 wbeT 中的单突变，所以 ogawa 要求 wbeT
  达到 99.9，而 inaba 是否定规则——失活突变（98% 一致性）仍能比对、仍算位点规则的
  present 基因，却过不了等位基因阈值。

两个实操约束：

- 探针记录必须在数据库中**独立存在**。gapit 以 abricate 的 `-culling_limit 1` 运行
  blastn（字节级兼容是特性），包含探针区的更长记录会遮蔽探针的比对——所以
  meningotype 夹具在探针旁没有全长的 synG 记录。
- 落在探针**之间**的等位基因（双重 W/Y 的 S 基序，或新的 wbeT 错义突变）两条阈值都
  清不过，落到 fallback——等位基因级别的歧义，以低置信度判定呈现。

当未来的 `allele_match` 原语落地（模式 + 等位基因表，像 meningotype 翻译其 EX7E 窗口
那样在比对到的查询序列上裁决）时，两个技巧都归于单条显式等位基因的规则。

## 校准与 learned_linear 展望

对照带标签的装配调阈值，正是
[`scripts/cluster_calibration.py`](https://github.com/indexofire/gapit/blob/main/scripts/cluster_calibration.py)
（开发者工具）的用途：指向一个已分型的基因簇数据库和标签 TSV，它打印每个期望表型的
分数分布、判定×期望一致性矩阵、以及分歧清单：

```console
$ python scripts/cluster_calibration.py mycps --datadir ./db labels.tsv
```

这份报告正是未来 `learned_linear` 训练流程要消费的循环：特征向量
（`gene:<id>:cov|ident|present`、`cluster:<locus>:coverage|identity`）恰是校准列，
训练出的模型作为又一条规则放进同一个 scheme 即可。

## 构建并分型手册方案

```console
$ gapit db build doumith doumith.fa --typing doumith.json
{"db":"doumith","records":5,"dbtype":"nucl","destination":".../doumith"}

$ gapit screen lm_4b.fa --db doumith --output lm_4b.tsv --quiet
$ gapit typing lm_4b.tsv --quiet
FILE	SCHEME	PHENOTYPE	CONFIDENCE	SCORE	RUNNER_UP	NOTES	GENES
lm_4b.fa	doumith	4b-4d-4e	high	1.0000	1/2a-3a (0.0000)	...	orf2110;orf2819;prs

$ gapit screen vc_inaba.fa --db vc_subserotype --quiet
FILE      BEST_LOCUS  TYPE  PHENOTYPE  COVERAGE  IDENTITY  PRESENT  PARTIAL  MISSING_IDS
vc_inaba.fa  wbe      O1    inaba      100.00    99.44     3        0        -
```

第二个代码块是基因簇例外：typed 基因簇数据库在 `gapit screen` 内直接判定（没有
第二阶段——不存在可供重读的筛查表）。

基因簇 O/K 模式按 scheme 拆分 panel：

```console
$ gapit db build vp_o vp_ok.gbk --typing vp_o_group.json
$ gapit db build vp_k vp_ok.gbk --typing vp_k_group.json
$ gapit screen sample.fa --db vp_o --quiet   # → O3/O13
$ gapit screen sample.fa --db vp_k --quiet   # → K6
```
