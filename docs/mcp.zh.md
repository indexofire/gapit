# MCP 服务器

gapit 内置一个 Model Context Protocol（MCP）stdio 服务器，把 CLI 暴露为工具。agent
运行时可以筛查 assembly 和 reads、汇总报告、自省 schema、自行准备数据库，无需
shell 外调或解析终端输出。协议是手写的换行分隔 JSON-RPC 2.0，零额外依赖：stdin 每
行一个请求，stdout 每行一个响应。

## 入口

| 命令 | 是什么 |
|---|---|
| `gapit mcp` | CLI 的 `mcp` 子命令。 |
| `gapit-mcp` | 随包安装的 console script（同一份代码）。 |

两者一直服务到 stdin 关闭（EOF 干净地结束进程）。stdout 只承载协议行；stderr 留给
协议内部错误。

## 客户端配置

用 MCP 客户端注册 console script，Claude 风格写法：

```json
{"mcpServers": {"gapit": {"command": "gapit-mcp"}}}
```

如果 `gapit-mcp` 不在客户端的 PATH 上，用绝对路径，或以
`{"command": "gapit", "args": ["mcp"]}` 调用模块 CLI。

## 工具

九个工具：五个只读分析工具，四个数据库自备工具（`db_fetch`、`db_build`、
`db_search`、`db_outdated`），与 `gapit db ...` 完全对应。下表参数列表转写自线上的
`tools/list` `inputSchema` 对象。

| 工具 | 参数 | 返回 |
|---|---|---|
| `screen` | `files`（字符串数组，必填）、`db`（string，必填 —— 无默认值）、`minid`（number）、`mincov`（number）、`format`（string：`json` \| `tsv` \| `md`，默认 `json`）、`aligner`（string：`blastn` \| `minimap2`，默认 `blastn`）、`mergeFragments`（boolean，默认 `false`；仅 blastn：跨 contig 片段合并）、`min_breadth`（number 0-100，默认 `90`；仅 minimap2）、`min_identity`（number 0-100，默认 `0`；仅 minimap2）、`min_mapq`（integer ≥ 0，默认 `0`；仅 minimap2）、`minGeneCov`（number 0-100，默认 `90`；仅 cluster 数据库）、`minGeneId`（number 0-100，默认 `90`；仅 cluster 数据库）、`minClusterCov`（number 0-100，默认 `96`；仅 cluster 数据库）、`datadir`（string） | 报告文本：默认 `gapit.report/1` JSON，按 `format` 也可为 TSV 或 Markdown；`aligner minimap2` 跑 reads 引擎并输出 `gapit.reads/1`。cluster 类型的 `db` 输出 `gapit.cluster/1`：最佳位点判定、逐基因判定，数据库携带 `typing.json` 时还有带 `phenotype_detail` 分数分解的表型判定 |
| `screen_reads` | `r1`（字符串数组，必填，每条 lane 一个路径）、`r2`（字符串数组，数量与 `r1` 相同）、`read_type`（string：`sr` \| `map-ont` \| `map-hifi`，默认 `sr`）、`min_breadth`（number 0-100，默认 `90`）、`min_identity`（number 0-100，默认 `0`）、`min_mapq`（integer ≥ 0，默认 `0`）、`format`（string：`json` \| `md`，默认 `json`）、`db`（string，必填 —— 无默认值）、`datadir`（string） | 默认 `gapit.reads/1` JSON，按 `format` 可为 Markdown；`min_identity`/`min_mapq` > 0 切换为 `gapit.reads/2` |
| `summary` | `files`（报告表路径数组，必填）、`identity`（boolean）、`nopath`（boolean） | `gapit.summary/1` JSON |
| `schema` | `name`（string，必填，取值 `cluster`、`error`、`features`、`reads`、`reads2`、`report`、`summary`、`typing`、`typing_result`、`version` 之一） | 该输出文档的 JSON Schema |
| `db_list` | 无 | `gapit.dblist/1`：数据库名、安装状态、记录数 |
| `db_fetch` | `name`（string，必填——数据库名，或 `all` 表示 card+vfdb 默认集合）、`datadir`（string）、`force`（boolean，默认 `false`） | 每个数据库一行 JSON 回执（`db`、`records`、`dbtype`、`destination`）。每个数据库都从其上游提供商下载，安装需要网络且可能耗时数分钟 |
| `db_build` | `name`（string，必填）、`fasta`（string，必填，本地文件系统路径）、`tsv`（string）、`dbtype`（string：`nucl` \| `prot`）、`description`（string）、`datadir`（string）、`force`（boolean，默认 `false`） | 一行 JSON 回执（`db`、`records`、`dbtype`、`destination`） |
| `db_search` | `term`（string，必填）、`db`（string）、`field`（string：`gene` \| `accession` \| `function` \| `product` \| `any`，默认 `any`）、`exact`（boolean，默认 `false`）、`limit`（integer ≥ 0，默认 `100`；`0` = 不限）、`datadir`（string） | TSV 命中行，列为 `DB`、`GENE`、`ACCESSION`、`FUNCTION`、`PRODUCT`、`LENGTH` |
| `db_outdated` | `days`（integer ≥ 0，默认 `90`）、`datadir`（string） | TSV 行，列为 `NAME`、`FETCHED_AT`、`AGE_DAYS`、`STATUS` |

省略的 `screen` 阈值回退到 CLI 默认值（`minid` 80、`mincov` 80；
`aligner minimap2` 时：`min_breadth` 90、`min_identity` 0、`min_mapq` 0），与
`gapit screen` 一致。`min_breadth`/`min_identity`/`min_mapq` 参数在默认 `blastn`
路径上被拒绝（它们是 reads 引擎参数）；`minimap2` 忽略 `minid`/`mincov`。

注意事项：

- 文件路径相对服务器的工作目录解析；绝对路径最稳妥。
- 数据目录来自 `GAPIT_DATADIR` 环境变量，然后 `~/.local/share/gapit/db`（见
  `./databases.md`）。db 工具也接受每次调用显式的 `datadir` 参数。
- `csv` 输出仅限 CLI：MCP `screen` 提供 `json`/`tsv`/`md`，`screen_reads` 提供
  `json`/`md`（reads 模式拒绝表格式输出）。
- `db_build` 的 `fasta` 必须是服务器文件系统上的路径：agent 提供本地路径，而非文件
  内容。
- 工具失败不使用 JSON-RPC 错误，而是返回 `isError: true`，文本内容为序列化的
  `gapit.error/1` 信封。

## 示例会话

用仓库的测试夹具搭一个用完即弃的数据目录（与测试套件同一配方），然后直接与服务器
对话：

```bash
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
export GAPIT_DATADIR=/tmp/gapit-mcp-demo/datadir
mkdir -p "$GAPIT_DATADIR"
cp -r tests/data/db/tinyamr "$GAPIT_DATADIR/"
makeblastdb -in "$GAPIT_DATADIR/tinyamr/sequences" -dbtype nucl \
  -out "$GAPIT_DATADIR/tinyamr/sequences" > /dev/null
```

向 `gapit mcp` 喂三行请求：

```bash
cat <<'EOF' | gapit mcp
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18"}}
{"jsonrpc":"2.0","id":2,"method":"tools/list"}
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"screen","arguments":{"files":["tests/data/contigs/full.fa"],"db":"tinyamr"}}}
EOF
```

响应第 1 行，原样：

```text
{"jsonrpc":"2.0","id":1,"result":{"protocolVersion":"2025-06-18","capabilities":{"tools":{}},"serverInfo":{"name":"gapit","version":"0.5.0"}}}
```

响应第 2 行（真实输出，中段省略；每个工具携带完整的 `inputSchema`）：

```text
{"jsonrpc":"2.0","id":2,"result":{"tools":[{"name":"screen","description":"Screen contig files for known genes (json = gapit.report/1; aligner minimap2 = fast assembly survey emitting gapit.reads/1).","inputSchema":{"type":"object","properties":{"files":{"type":"array","items":{"type":"string"}},"db":{"type":"string"}, ... },"required":["files","db"]}}, {"name":"screen_reads","description":"Screen FASTQ reads for genes via minimap2 (json = gapit.reads/1; min_identity/min_mapq > 0 emits gapit.reads/2). Returns the rendered document.","inputSchema":{"type":"object","properties":{"r1":{"type":"array","items":{"type":"string"}},"r2":{"type":"array","items":{"type":"string"}},"read_type":{"type":"string","enum":["sr","map-ont","map-hifi"],"default":"sr"}, ... },"required":["r1","db"]}}, {"name":"summary", ...}, {"name":"schema", ...}, {"name":"db_list", ...}, {"name":"db_fetch", ...}, {"name":"db_build", ...}, {"name":"db_search", ...}, {"name":"db_outdated", ...}]}}
```

响应第 3 行（真实输出，text 内容省略）。screen 结果以 JSON 字符串的形式位于
`result.content[0].text`：

```text
{"jsonrpc":"2.0","id":3,"result":{"content":[{"type":"text","text":"{\n  \"schema\": \"gapit.report/1\", ... }"}],"isError":false}}
```

解码后，那段文本是一份完整的 `gapit.report/1` 文档（同一次运行的真实输出）：

```json
{
  "schema": "gapit.report/1",
  "tool": {"name": "gapit", "version": "0.5.0"},
  "created_at": "2026-09-19T01:15:44Z",
  "params": {"db": "tinyamr", "minid": 80.0, "mincov": 80.0, "threads": 1},
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

### 筛查 reads

`screen_reads` 对应 `gapit screen --r1/--r2`：每条 lane 一个 minimap2，存在判定看
比对广度。数据目录约定与 `screen` 相同，但不需要 BLAST 索引：minimap2 直接读取
sequences 文件。对仓库 `tinyreads` 夹具的临时副本的真实会话：

```bash
mkdir -p /tmp/gapit-mcp-demo/datadir
cp -r tests/data/reads_db/tinyreads /tmp/gapit-mcp-demo/datadir/
cat <<'EOF' | gapit mcp
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"screen_reads","arguments":{"r1":["tests/data/reads/tetx_full.fq"],"db":"tinyreads","datadir":"/tmp/gapit-mcp-demo/datadir"}}}
EOF
```

原样响应（text 内容省略）：

```text
{"jsonrpc":"2.0","id":3,"result":{"content":[{"type":"text","text":"{\n  \"schema\": \"gapit.reads/1\", ... }"}],"isError":false}}
```

解码后，文本是一份完整的 `gapit.reads/1` 文档（同一次运行的真实输出）：

```json
{
  "schema": "gapit.reads/1",
  "tool": {"name": "gapit", "version": "0.5.0"},
  "created_at": "2026-09-22T00:07:35Z",
  "params": {
    "db": "tinyreads",
    "read_type": "sr",
    "min_breadth": 90.0,
    "threads": 1
  },
  "files": [
    {
      "reads": [
        "tests/data/reads/tetx_full.fq"
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

`r2` 接受 mate 路径（数量与 `r1` 相同），`read_type` 选择 minimap2 预设，非零的
`min_identity`/`min_mapq` 把文档切换为 `gapit.reads/2`（语义见 `./reads.md`）。
assembly FASTA 请改用 `aligner: "minimap2"` 调 `screen`：普查会自动解析 `map-ont`
预设。

### 失败形态

缺失的输入文件返回工具错误，而非协议错误（真实行）：

```text
{"jsonrpc":"2.0","id":5,"result":{"content":[{"type":"text","text":"{\"schema\":\"gapit.error/1\",\"code\":\"INPUT_NOT_FOUND\",\"message\":\"input file not found or unreadable: /tmp/gapit-mcp-demo/nope.fa\",\"context\":{\"file\":\"/tmp/gapit-mcp-demo/nope.fa\"}}"}],"isError":true}}
```

把 `content[0].text` 解析为 JSON，按 `code` 分支（错误码与退出码映射：
`./outputs.md`）。

### 自备数据库：先构建，后筛查

db 工具让 agent 在会话中途自备数据库。对新数据目录的真实会话（`my_genes.fa` 是一
条 240 bp 合成基因，`query.fa` 是它的 200 bp 子串）：

```bash
mkdir -p /tmp/gapit-mcp-demo/datadir
export GAPIT_DATADIR=/tmp/gapit-mcp-demo/datadir
```

两行请求及其原样响应：

```text
{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"db_build","arguments":{"name":"myamr","fasta":"/tmp/gapit-mcp-demo/my_genes.fa"}}}
{"jsonrpc":"2.0","id":2,"result":{"content":[{"type":"text","text":"{\"db\":\"myamr\",\"records\":1,\"dbtype\":\"nucl\",\"destination\":\"/tmp/gapit-mcp-demo/datadir/myamr\"}"}],"isError":false}}
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"screen","arguments":{"files":["/tmp/gapit-mcp-demo/query.fa"],"db":"myamr"}}}
{"jsonrpc":"2.0","id":3,"result":{"content":[{"type":"text","text":"{ ... gapit.report/1 ... }"}],"isError":false}}
```

解码后，screen 的文本是一份完整的 `gapit.report/1` 文档，唯一命中正是 agent 刚构建
进数据库的那个基因（真实输出）：

```json
{
  "schema": "gapit.report/1",
  "tool": {"name": "gapit", "version": "0.5.0"},
  "created_at": "2026-09-20T13:38:39Z",
  "params": {"db": "myamr", "minid": 80.0, "mincov": 80.0, "threads": 1},
  "files": [
    {
      "file": "/tmp/gapit-mcp-demo/query.fa",
      "hits": [
        {
          "sequence": "contig1",
          "start": 1,
          "end": 200,
          "strand": "+",
          "gene": "demov2",
          "coverage": "1-200/240",
          "coverage_map": "=============..",
          "gaps": "0/0",
          "coverage_pct": 83.33,
          "identity_pct": 100.0,
          "database": "myamr",
          "accession": "",
          "product": "demo beta-lactamase variant 2",
          "resistance": ""
        }
      ]
    }
  ]
}
```

同样的会话模式也适用于 `db_fetch`（`name` 传 `all` 时安装 card+vfdb 默认集合）以及用于
巡检的 `db_search`/`db_outdated`。自定义数据库构建规则（表头类型、`--tsv` 元数
据）：`./custom-db.md`。

## 协议注意事项

- **分帧。** stdin 每行一个 JSON-RPC 2.0 消息，stdout 每个请求一行响应，UTF-8、LF
  行尾。响应使用紧凑分隔符。
- **initialize。** 原样回显请求的字符串 `protocolVersion`；缺失或为空时以内置默认
  `2025-06-18` 应答。结果携带 `capabilities.tools` 和 `serverInfo`（`name: gapit`，
  包版本）。
- **方法。** 只处理三个：`initialize`、`tools/list`、`tools/call`。其他一律返回
  JSON-RPC 错误 `-32601`：

  ```text
  {"jsonrpc":"2.0","id":4,"error":{"code":-32601,"message":"method not found: resources/list"}}
  ```

- **未知工具。** `tools/call` 指向未注册的工具，或参数非法时，返回 `-32602`：

  ```text
  {"jsonrpc":"2.0","id":6,"error":{"code":-32602,"message":"unknown tool: nope"}}
  ```

- **不支持批处理数组。** JSON 数组行不是合法帧，会被丢弃；请每行一个请求。
- **通知静默。** 无 `id` 的帧（例如 `notifications/initialized`）从不产生响应。
- **垃圾容忍。** 非 JSON 行被忽略，不回解析错误，循环继续服务。EOF 干净终止。
