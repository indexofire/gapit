# gapit

[English](README.md) | [中文](README.zh.md)

对 contig 和 reads 中的已知基因进行大规模筛查——覆盖 AMR、毒力、血清型、质粒复制子、
物种靶标或任意自定义数据库。contig 筛查输出 abricate 兼容 TSV，同时把 JSON 和
Markdown 作为一等输出。支持 Python 3.12。

## 安装

```bash
pip install gapit
```

需要在 `PATH` 上提供 BLAST+（`blastn`、`blastx`、`makeblastdb`、`blastdbcmd`）和
`minimap2`；pixi / conda 安装方式会自动管理这些二进制。pixi、conda、Docker 和固定版本
安装见[安装文档](docs/installation.zh.md)。

```bash
# 从 git clone 安装（自动管理外部二进制）
git clone https://github.com/indexofire/gapit.git
cd gapit
pixi install
pixi run gapit --version
```

## 快速开始

内置数据库首次使用时离线物化；其他数据库通过 `gapit db fetch` 安装。

```bash
# contig 筛查 -> stdout 输出 abricate 兼容 TSV
gapit screen contigs.fa --db eco_dec

# 面向智能体和人类的 JSON / Markdown
gapit screen contigs.fa --db eco_dec --format json
gapit screen contigs.fa --db eco_dec --format md

# 通过 minimap2 筛查 FASTQ reads
gapit screen --r1 R1.fastq.gz --r2 R2.fastq.gz --read-type sr --db eco_dec

# 两阶段基因判定：screen 写表，typing 从表判定
gapit screen isolate.fa --db eco_dec --output result.tsv
gapit typing result.tsv

# 汇总多份报告为基因矩阵
gapit summary *.tsv

# 数据库
gapit db fetch -d ncbi
gapit db list

# 查看输出契约
gapit schema report
```

contig 默认值为 `--minid 80 --mincov 80`；reads 存在判定默认要求 90% 比对广度。
完整参数与示例见 [contig 筛查](docs/screen.zh.md)、[reads](docs/reads.zh.md)、
[typing](docs/typing.zh.md)、[summary](docs/summary.zh.md)。

## 文档

| 主题 | 页面 |
|---|---|
| 安装、快速上手 | [安装](docs/installation.zh.md)、[快速上手](docs/quickstart.zh.md) |
| contig 筛查、cluster 数据库 | [筛查](docs/screen.zh.md) |
| reads 与 assembly 调查 | [reads](docs/reads.zh.md) |
| 基因矩阵 | [summary](docs/summary.zh.md) |
| typing 方案与表型 | [typing](docs/typing.zh.md) |
| JSON/Markdown 契约、错误信封 | [输出](docs/outputs.zh.md)、`gapit schema` |
| 数据库目录、提供商、内置库 | [数据库](docs/databases.zh.md) |
| datadir、共享、Docker | [数据目录](docs/datadir.zh.md) |
| 自定义数据库 | [自定义数据库](docs/custom-db.zh.md) |
| 智能体与 MCP 服务器 | [智能体](docs/agents.zh.md)、[MCP](docs/mcp.zh.md) |
| 可选 gseek 引擎 | [gseek](docs/gseek-engine.zh.md) |

完整手册：<https://indexofire.github.io/gapit/>。

## 开发与许可证

贡献流程与门禁见 [CONTRIBUTING.md](CONTRIBUTING.md)；内部规则见
[AGENTS.md](AGENTS.md)；信任边界见 [SECURITY.md](SECURITY.md)；变更历史见
[CHANGELOG.md](CHANGELOG.md)。

gapit 使用 MIT 许可。它是 abricate（GPL-2.0）的行为重实现，不复制任何 Perl 代码。
七个由 gapit 策划的内置数据库均为宽松许可；其他提供商在 `db fetch` 时从上游下载，
遵循各自许可证。使用 gapit 时请引用 [CITATION.cff](CITATION.cff)。
