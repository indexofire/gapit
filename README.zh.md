# gapit

[English](README.md) | [中文](README.zh.md)

对 contig 和 reads 中的已知基因进行大规模筛查——覆盖 AMR、毒力、血清型、质粒复制子、
物种靶标或任意自定义数据库。这是 [abricate](https://github.com/tseemann/abricate) 的
面向智能体（agent-first）Python 重实现：提供字节级兼容的 TSV，同时把 JSON 和 Markdown
作为一等输出。

## 为什么选择 gapit

- **冻结的 abricate 兼容基线。** contig 筛查保留与 abricate 相同的 BLAST 流程、命中规则
  和 abricate 格式 TSV，作为现有流水线的冻结基线。从 v0.5.0 起，gapit 走自己的契约；
  abricate 字节级一致性不再是发布门槛——parity 对照工具
  （`pixi run -e parity parity`）仍然作为这套冻结表面的回归参照。
- **默认机器可读。** 版本化输出 schema（`gapit.report/1`、`gapit.reads/1`、
  `gapit.summary/1`），可通过 `gapit schema` 自描述；stderr 输出类型化 JSON 错误信封，
  并有文档化的退出码。智能体无需阅读文档即可发现完整契约。
- **不只支持 contig，也支持 reads。** `gapit screen --r1/--r2` 通过 minimap2 筛查
  FASTQ；`--r1` 也直接接受 assembly FASTA（按内容识别，强制 `map-ont`），用于快速
  存在性调查。abricate 无法筛查原始 reads。

## 安装

从 PyPI 安装（需要 PATH 上有外部 BLAST+ 和 minimap2 二进制——见
[docs/installation.zh.md](docs/installation.zh.md)）：

```bash
pip install gapit
```

或者用 [pixi](https://pixi.sh) 从 git clone 安装；pixi 会连同外部二进制一起管理环境：

```bash
git clone https://github.com/indexofire/gapit.git
cd gapit
pixi install
pixi run gapit --version
```

也可以用仓库里的 `environment.yml`，把 gapit 和外部二进制装进 conda 环境：

```bash
git clone https://github.com/indexofire/gapit.git
cd gapit
conda env create -f environment.yml
conda activate gapit-env
gapit --version
```

本包只支持 Python 3.12（仅 3.12）；开发环境固定在 3.12。

## 快速开始

```bash
# 筛查 contig（stdout 输出 abricate 兼容 TSV；--db 必填）
gapit screen contigs.fa --db ncbi

# 面向智能体和人类的输出
gapit screen contigs.fa --db ncbi --format json
gapit screen contigs.fa --db ncbi --format md

# 两阶段判定：screen 写基因表，typing 从表中做判定
# （带 typing 的基因数据库——内置 ecoli_dec 携带双方案 DEC panel）
gapit screen isolates.fa --db ecoli_dec --output dec.tsv
gapit typing dec.tsv
# 脑膜炎奈瑟菌：血清群 + ctrA + sodC 内参（不含 MLST）
gapit db fetch nme_sero
gapit screen isolate.fna -d nme_sero

# 筛查 FASTQ reads；--read-type 选择 minimap2 preset
gapit screen --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --read-type sr --db card
gapit screen --r1 ont_reads.fastq.gz --db ncbi --read-type map-ont --format json

# 把报告表汇总成基因存在/缺失矩阵
gapit summary *.tsv

# 数据库：七个随 wheel 内置、开箱即用；其余在 fetch 时从上游下载
gapit db fetch all         # 默认集合（card、vfdb）
gapit db fetch ncbi
gapit db list
gapit db outdated                # 标记超过过期阈值的数据库
gapit db search "tet(M)"         # 在所有已安装数据库中查找基因
gapit db build mydb my_genes.fa --tsv my_meta.tsv   # 从任意 FASTA 构建自定义库

# 自省
gapit schema report     # gapit.report/1 的 JSON Schema
```

contig 筛查沿用 abricate 默认值（`--minid 80`、`--mincov 80`）。reads 模式的存在判定
默认要求 90% 比对广度（`--min-breadth 90`）。

## 数据库

数据库是 datadir 下的一个目录。datadir 的解析顺序是 `--datadir`，然后
`$GAPIT_DATADIR`，然后 `$XDG_DATA_HOME/gapit/db`，最后
`~/.local/share/gapit/db`（可用 `gapit db path` 查看）。全新安装时，
`gapit db list` 会引导创建 datadir 并列出目录清单（内置库行在首次使用前即可见）；
`gapit db setup` 会创建 datadir 并物化全部七个内置数据库。

目录中共有二十五个数据库——十三个基因库提供商、八个基因簇提供商，以及七个内置快照
（`cje_cps`、`ecoh`、`ecoli_dec`、`lmo_sero`、`ncbi`、`resfinder`、`vpa_ok`；
其中 `ecoh`/`ncbi`/`resfinder` 同时也是提供商名）。每个提供商都通过
`gapit db fetch` 构建，从上游下载——包括 `nme_sero`（脑膜炎奈瑟菌血清群 + ctrA +
sodC，非 MLST）——另有七个经过内容来源审计的库随包内置（`ecoli_dec`、`lmo_sero`
以及 `ncbi`、`resfinder`、`ecoh` 快照——公共领域 / Apache-2.0 / BSD-3-Clause / MIT——
以及 MIT 的 typed-cluster 快照 `vpa_ok`），首次使用时物化进 datadir，零网络；内置
`ecoli_dec` 还携带逐基因一致性下限（`floors.json`，`pic` 为 90%），两个引擎在做存在
判定前都会应用；`cje_cps` 则同时携带 floors 和每数据库的 `screen.json` 默认值
（MAPQ 20、depth 5、10% lpxA 相对深度）。同一套 `screen.json`（`gapit.screen/1`）
机制允许任何数据库声明自己的 reads/contigs 阈值，因此不会把一套全局参数强加给扩增子
长度不同的 panel。若干上游许可（CARD 的 McMaster 条款、VFDB 的 CC BY-NC、Kaptive 的
GPL-3.0）禁止随 MIT 许可的发行版再分发，所以其他内容都按需抓取——`db fetch <name>`
会从最新上游刷新由提供商支持的内置库名（`ecoh`、`ncbi`、`resfinder`）。数据库结构、
datadir 共享、conda/Docker 指南和安全注意事项见[数据目录](docs/datadir.zh.md)。

| 名称 | 内容 | dbtype |
|---|---|---|
| `ncbi` | NCBI AMRFinderPlus 精选 AMR | nucl |
| `card` | CARD 蛋白同源耐药模型 | nucl |
| `resfinder` | CGE ResFinder 获得性耐药基因 | nucl |
| `argannot` | ARG-ANNOT 获得性耐药基因 | nucl |
| `plasmidfinder` | CGE PlasmidFinder 复制子 | nucl |
| `megares` | MEGARes 抗微生物耐药基因 | nucl |
| `ecoh` | E. coli O 与 H 抗原（srst2 EcOH） | nucl |
| `ecoli_dec` | 致泻性大肠埃希菌标记 panel（GB 4789.6 + 风险监测判定） | nucl |
| `vfdb` | VFDB 毒力因子（set A，核酸） | nucl |
| `ecoli_vf` | E. coli 毒力因子 | nucl |
| `lmo_sero` | 单核细胞增生李斯特菌血清群（Doumith 2004） | nucl |
| `bacmet2` | BacMet2 杀菌剂/耐药基因（蛋白） | prot |
| `victors` | Victors 毒力因子 | nucl |
| `upec_expec_vf` | UPEC/ExPEC 毒力基因 | nucl |
| `vpa_ok` | 副溶血弧菌 O 与 K 抗原位点（OAgc/CPSgc，vpasero；cluster，MIT，一次扫描 O:K 血清型） | nucl |
| `vch_oagc` | 霍乱弧菌 O 抗原基因簇（190 个 OL locus + ompW 内参；cluster，按需抓取） | nucl |

gapit 也能读取 abricate 自己构建的 datadir（旧式 `~~~` 表头）。反过来不成立：
gapit 原生数据库使用 `gapit/v1` 表头格式（见 SPEC.md §11），abricate 无法读取。
`bacmet2` 等蛋白数据库通过 blastx 筛查。

基因**簇**数据库是第二种类型：`gapit db build NAME loci.gbk|gff3` 构建
`kind: cluster` 数据库（通过 minimap2 做 locus call，可选 `--typing FILE` 做表型
评分），七个 Kaptive 基因簇数据库（kpsc_k、kpsc_o、kosc_k、kosc_o、ab_k、ab_o、
ecoli_kps）在安装时抓取 Kaptive 抗原 locus 参考（GPL-3.0 内容，因此不内置——使用
结果时请引用 Kaptive/Wyres et al. 2020）。带 typing 的基因簇内置库 `vpa_ok` 把
vpasero 的 OAgc + CPSgc 参考（MIT）合并成一个内置库——这是 ecoli_dec/lmo_sero 模式
的纯内置库，不需要 `db fetch`——单次 `gapit screen --db vpa_ok` 即可渲染
`gapit.cluster/2`，给出 `o_group`、`k_group` 以及组合的 `serotype` O:K 判定
（RIMD 2210633 → O3:K6），`pip install gapit` 后零网络即可工作。
`gapit db list --json` 会在提供商声明了内容许可证时暴露该许可证。

## 输出契约

- **stdout 纯净。** 数据写 stdout，诊断写 stderr，始终如此。`--quiet` 只静默 stderr。
- **确定性。** 稳定的排序、固定的工具参数、数据载荷内没有 wall-clock 时间戳。
- **错误**以 JSON 信封写到 stderr，并以非零码退出：

```text
{"schema": "gapit.error/1", "code": "...", "message": "...", "context": {...}}
```

| 退出码 | 含义 |
|---|---|
| 0 | 成功 |
| 1 | 未预期错误 |
| 2 | 用法错误 |
| 3 | 缺少依赖 |
| 4 | 数据库错误 |
| 5 | 输入错误 |

## MCP 服务器

gapit 内置 MCP（Model Context Protocol）stdio 服务器，智能体运行时无需解析 CLI 输出
即可筛查 assembly：`gapit mcp` 或 `gapit-mcp` console script 在 stdin/stdout 上使用
换行分隔的 JSON-RPC 2.0（无额外依赖——协议为手写实现）。它暴露十二个工具：
`screen`（默认 `gapit.report/1`；`aligner minimap2` 用于快速 assembly 调查）、
`screen_reads`（通过 minimap2 处理 FASTQ，`gapit.reads/1`）、`summary`
（`gapit.summary/1`）、`typing`（从筛查表做判定，并给出 W/Y EX7E 检查等后续提示）、
`schema`、`db_list`，以及数据库工具 `db_path`、`db_setup`、`db_fetch`、`db_build`、
`db_search` 和 `db_outdated`，让智能体可以在会话中自行准备和检查数据库。
工具失败时返回 `isError: true`，文本内容是 `gapit.error/1` 信封。
用 MCP 客户端注册：

```json
{"mcpServers": {"gapit": {"command": "gapit-mcp"}}}
```

## 可选的 gseek 引擎

`GAPIT_ENGINE=legacy|gseek|shadow` 可以把 contig BLAST、蛋白 `blastx-fast`、reads
minimap2 和基因簇 minimap2 路由到 [gseek](https://github.com/indexofire/gseek)。
`legacy` 仍是默认值，并在所有路径上保留；`shadow` 同时运行两者并返回 legacy 结果。
gseek 需要 Python 3.12，可通过 pixi `default` 环境（`py312` feature）安装，或
`pip install "gapit[gseek]"`。见 [gseek 引擎](docs/gseek-engine.zh.md)。

## 开发

| 任务 | 运行内容 |
|---|---|
| `pixi run lint` | ruff check |
| `pixi run fmt` | ruff format |
| `pixi run fmt-check` | ruff format --check |
| `pixi run typecheck` | basedpyright（strict） |
| `pixi run test` | pytest，离线单元/CLI 套件 |
| `pixi run -e parity parity` | 与真实 abricate 做筛查字节级对比 |
| `pixi run -e parity summary-parity` | 与真实 abricate 做 summary 字节级对比 |

用户文档在 [`docs/index.zh.md`](docs/index.zh.md)，渲染于
<https://indexofire.github.io/gapit/>。

`SPEC.md` 是 parity 契约，`AGENTS.md`/`CONTRIBUTING.md` 是贡献规则，
`SECURITY.md` 是信任边界与报告策略，`CITATION.cff` 是引用元数据，
`CHANGELOG.md` 是变更历史。`PLAN.md` 作为历史路线图保留（当前架构见
`AGENTS.md`、`SPEC.md` 和 `docs/datadir.zh.md`）。

## 许可证

gapit 使用 MIT 许可。它是 abricate（GPL-2.0）的行为重实现，不复制任何 Perl 代码；
abricate 本身仍为 GPL-2.0。包内唯一的数据库内容是 `src/gapit/data/dbs/` 下七个经过
审计、许可宽松的内置库（公共领域 / Apache-2.0 / BSD-3-Clause / MIT）；其他每个提供商
都在抓取时从上游下载，遵循其各自许可证（SPEC.md §9）。
