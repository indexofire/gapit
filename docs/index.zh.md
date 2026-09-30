# gapit 文档

gapit 筛查 contig 文件和 FASTQ reads 中的已知基因：AMR（抗微生物药物耐药性）、毒力、
血清型、质粒复制子、物种靶点，或任意自定义数据库。它是 [abricate](https://github.com/tseemann/abricate)
的 agent 优先 Python 重实现：stdout 输出 abricate 兼容的 TSV，同时提供带版本号 schema
的一等 JSON 与 Markdown 输出。

## 页面

| 页面 | 内容 |
|---|---|
| [安装](./installation.md) | PyPI wheel 与 pixi 两条安装路径、外部二进制依赖、安装验证、数据库引导、shell 补全 |
| [快速开始](./quickstart.md) | 在离线环境下用测试夹具数据库完成一次完整的首个会话 |
| [contig 筛查](./screen.md) | 对 FASTA/GBK/EMBL 输入运行 `gapit screen`：阈值、过滤器、输出格式 |
| [reads（FASTQ）筛查](./reads.md) | FASTQ 与 assembly FASTA 走 minimap2：`--r1`/`--r2`、预设参数、基于广度的存在判定、两阶段普查 |
| [汇总报告](./summary.md) | `gapit summary`：把多份报告表折叠成基因存在/缺失矩阵 |
| [数据库](./databases.md) | 数据库目录（NAME 与上游 PROVIDER）、`gapit db fetch/list/install`、数据目录、原生 `gapit/v1` 格式、基因簇数据库（GBK/GFF、typing.json、kaptive） |
| [自定义数据库](./custom-db.md) | `gapit db build` 实战演练：任意 FASTA 变成可筛查的数据库，完整示例 |
| [输出](./outputs.md) | TSV/JSON/Markdown 格式、schema、错误信封、退出码 |
| [MCP 服务器](./mcp.md) | `gapit mcp`：通过 stdio JSON-RPC 为 agent 运行时提供分析 + 数据库工具 |
| [面向 agent 的指南](./agents.md) | 自主 agent 如何消费 gapit：schema、自省、错误处理 |
| [常见问题](./faq.md) | 常见问题、与 abricate 的差异、故障排查 |

## 命令

| 命令 | 作用 | 文档 |
|---|---|---|
| `gapit screen` | 筛查 contig 文件或 FASTQ reads 中的已知基因 | [筛查](./screen.md) |
| `gapit summary` | 把一份或多份报告表汇总成基因存在/缺失矩阵 | [汇总报告](./summary.md) |
| `gapit db fetch` | 抓取并构建数据库到数据目录（kaptive 基因簇数据库在抓取时下载） | [数据库](./databases.md) |
| `gapit db list` | 列出已知数据库（NAME、上游 PROVIDER、安装状态） | [数据库](./databases.md) |
| `gapit db install` | 校验 SHA256 后安装本地文件 | [数据库](./databases.md) |
| `gapit db build` | 从基因 FASTA 构建数据库，或从 GBK/GFF 位点构建基因簇数据库（可选配 `--typing` 表型评分规范） | [自定义数据库](./custom-db.md) |
| `gapit setupdb` | 为数据目录下的所有数据库构建 BLAST 索引 | 本页 |
| `gapit schema` | 打印 gapit 输出文档的 JSON Schema（`report`、`reads`、`reads2`、`cluster`、`summary`、`error`、`version`、`features`、`typing`） | [输出](./outputs.md) |
| `gapit mcp` | 运行 MCP stdio 服务器（同时安装为 `gapit-mcp` console script） | [MCP 服务器](./mcp.md) |

## 项目文档

| 文件 | 内容 |
|---|---|
| [README](https://github.com/indexofire/gapit/blob/main/README.md) | 项目概览、安装、快速上手、输出契约、数据库表 |
| [SPEC.md](https://github.com/indexofire/gapit/blob/main/SPEC.md) | abricate 行为规范；比对一致性的唯一事实来源 |
| [PLAN.md](https://github.com/indexofire/gapit/blob/main/PLAN.md) | 按阶段划分的开发路线图 |
| [AGENTS.md](https://github.com/indexofire/gapit/blob/main/AGENTS.md) | 贡献者与 agent 指南：目录结构、约定、验证 |
| [CHANGELOG.md](https://github.com/indexofire/gapit/blob/main/CHANGELOG.md) | 变更历史 |
