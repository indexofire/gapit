# 面向 agent 的 gapit

gapit 生来就是被程序驱动的。每个输出都以机器可读为设计目标，整个契约可以从二进制
本身发现，失败是类型化的 JSON 而不是散文。本页就是契约：集成时你不需要其他任何文
档。

## 无需文档的发现

三条命令就能回答"这东西能干什么"，不用读手册。

版本，一行 JSON：

```console
$ gapit --version --json
{"schema":"gapit.version/1","name":"gapit","version":"0.5.0"}
```

输出 schema。十个文档可自省：`report`、`typing_result`、`reads`、`reads2`、`cluster`、
`summary`、`error`、`version`，外加数据库侧文档 `features` 和 `typing`。每个都打印完整的
JSON Schema：

```console
$ gapit schema report          # gapit.report/1 (contig screening)
$ gapit schema typing_result   # gapit.typing_result/1 (gapit typing 判定)
$ gapit schema reads           # gapit.reads/1  (FASTQ screening)
$ gapit schema reads2          # gapit.reads/2  (filtered FASTQ screening)
$ gapit schema cluster         # gapit.cluster/1 (cluster-database screening)
$ gapit schema summary         # gapit.summary/1
$ gapit schema error           # gapit.error/1
$ gapit schema version         # gapit.version/1
$ gapit schema features        # gapit.features/1 (cluster db feature table)
$ gapit schema typing          # gapit.typing/2 (db scoring spec; /1 still reads)
```

例如 `gapit schema version`（真实输出）：

```json
{
  "description": "gapit.version/1 \u2014 compact one-line self-description.",
  "properties": {
    "schema": {
      "const": "gapit.version/1",
      "default": "gapit.version/1",
      "title": "Schema",
      "type": "string"
    },
    "name": {
      "default": "gapit",
      "title": "Name",
      "type": "string"
    },
    "version": {
      "title": "Version",
      "type": "string"
    }
  },
  "required": ["version"],
  "title": "VersionDocument",
  "type": "object"
}
```

数据库：数据库目录与安装状态（`gapit db list --json`；真实输出，十二个中的前两
个，有删节）：

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
    }
  ]
}
```

## 确定性输出

相同输入，相同字节。成立的原因：

- 排序固定且有文档（`./screen.md`、`./summary.md`）。
- 数据载荷内没有墙钟时间戳。唯一的时间字段是文档元数据（`created_at`）。
- LF 行尾、UTF-8、每个 schema 稳定的键顺序。
- **stdout 纯净**：数据走 stdout，诊断走 stderr，永远如此。把 stdout 重定向到文
  件，你得到的就是文档本身，别无他物。

## 退出码即控制流

退出时的整数在你读 stderr 的任何字节之前就告诉你失败类别。

| 退出码 | 含义 |
|---|---|
| 0 | 成功 |
| 1 | 意外错误 |
| 2 | 用法错误 |
| 3 | 缺少依赖 |
| 4 | 数据库错误 |
| 5 | 输入错误 |

信封细节与错误码：`./outputs.md`。

## JSON 稳定性策略

- 每个文档用版本字符串自报家门：`gapit.report/1`、`gapit.reads/1`、
  `gapit.cluster/1`、`gapit.summary/1`、`gapit.error/1`、`gapit.version/1`。先查
  `schema`，再按它分发。
- schema 遵循 semver。minor 变更绝不重命名或改型既有字段；新字段可能出现，所以请
  忽略未知键而不是拒绝它们。
- 键是 snake_case，单位显式（`identity_pct`、`coverage_pct`）。
- 权威的机器契约就是 schema 本身：`gapit schema <doc>` 永远是最新的，本页可能滞
  后。

## 错误信封解析

任何失败都会向 stderr 打印恰好一行 JSON：

```json
{"schema":"gapit.error/1","code":"...","message":"...","context":{...}}
```

真实示例，对一个不存在的数据库筛查（退出码 4）：

```text
{"schema":"gapit.error/1","code":"DATABASE_NOT_FOUND","message":"Database nosuchdb is not in /tmp/gapit-mcp-demo/datadir. Available: tinyamr","context":{"db":"nosuchdb","datadir":"/tmp/gapit-mcp-demo/datadir"}}
```

缺外部二进制（退出码 3）：

```text
{"schema":"gapit.error/1","code":"MISSING_DEPENDENCY","message":"required binary not found on PATH: blastn","context":{"binary":"blastn"}}
```

解析建议：

- 进度诊断（例如 `Processing: sampleA.fa`）也走 stderr。信封是 stderr 中唯一能解析
  为带 `"schema":"gapit.error/1"` 的 JSON 的行。
- `code` 是稳定的字符串枚举（`USAGE_ERROR`、`MISSING_DEPENDENCY`、
  `DATABASE_NOT_FOUND`、`INPUT_NOT_FOUND` 等）。按 `code` 分支，不要按 `message`。
- `context` 携带机器可读细节（文件路径、数据库名），可以安全地记录或呈现。

## 推荐工作流

发现、筛查、解析、汇总。下面的会话对仓库测试夹具数据目录运行（搭建配方见
`./mcp.md`）。

1. 发现可用数据库（`gapit db list --json`，或 `db_list` MCP 工具），选一个 `db`
   名。

2. 以 JSON 筛查每个样本，直接从文档解析命中：

   ```console
   $ gapit screen sampleA.fa --db tinyamr --format json > sampleA.json
   $ gapit screen sampleB.fa --db tinyamr --format json > sampleB.json
   ```

   ```python
   import json

   for path in ("sampleA.json", "sampleB.json"):
       doc = json.load(open(path))
       assert doc["schema"] == "gapit.report/1"
       for f in doc["files"]:
           for hit in f["hits"]:
               print(f["file"], hit["gene"], hit["identity_pct"], hit["coverage_pct"])
   ```

   真实输出：

   ```text
   sampleA.fa tetA 100.0 100.0
   ```

3. 需要跨样本矩阵时，保留默认的 TSV 报告，喂给 `gapit summary`：

   ```console
   $ gapit screen sampleA.fa --db tinyamr > sampleA.tsv
   $ gapit screen sampleB.fa --db tinyamr > sampleB.tsv
   $ gapit summary sampleA.tsv sampleB.tsv
   #FILE	NUM_FOUND	tetA
   sampleA.tsv	1	100.00
   sampleB.tsv	0	.
   ```

筛查选项：`./screen.md`。汇总语义：`./summary.md`。

## 工具运行时集成

如果你的运行时讲 MCP，把 gapit 注册为服务器，用 `screen`、`screen_reads`、
`summary`、`schema` 和 `db_list` 作为工具而不是子进程：`./mcp.md`。服务器也暴露
`db` 命令（`db_fetch`、`db_build`、`db_search`、`db_outdated`），agent 可以在同一个
会话里准备好所需数据库并完成筛查。
