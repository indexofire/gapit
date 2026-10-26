# 数据目录与数据库布局

gapit 刻意把**工具环境**和**数据库状态**分开。本页说明数据库放在哪里、一个数据库目录里有什么、
为什么默认不放进 conda 环境，以及多环境共享和容器场景的规范做法。

如果只想知道路径：`gapit db path` 打印当前生效的 datadir；`gapit db list` 会先创建该目录，
再列出数据库目录清单。

## 快速选择表

| 场景 | 推荐 datadir |
|---|---|
| 单用户、多个 conda/env | 默认用户级目录（共享缓存）；升级 gapit 后按需刷新 |
| 项目需要可复现的数据库版本 | 每个项目显式 `--datadir` 或 `GAPIT_DATADIR` |
| 长期工作站 | `$XDG_DATA_HOME/gapit/db`，否则 `~/.local/share/gapit/db` |
| CI | 缓存或挂载 workspace 数据目录；固定数据库版本 |
| Docker / Podman | named volume 或 bind mount + `GAPIT_DATADIR`；不要依赖容器内 `$HOME` |
| 多用户 HPC | 每用户独立 datadir；共享只读参考库需要显式权限模型 |
| 严格安全边界 | 不要跨用户或跨不可信任务共享可写 datadir |

## 什么是 gapit 数据库？

一个数据库就是 datadir 下的一个**子目录**，一个数据库名一个目录：

```text
<datadir>/<name>/
├── sequences                  # 必需：BLAST/minimap2 使用的 FASTA
├── gapit-manifest.json        # 来源与认证信息（gapit.manifest/1）
├── records.jsonl              # 基因 truth store（gapit 构建的 gene 库）
├── features.json              # 位点/基因表（cluster 库）
├── typing.json                # 可选：gapit.typing/1 或 /2 规则
├── floors.json                # 可选：按基因的 reads 一致性下限
└── sequences.n* / sequences.p* / sequences.nhr ...   # BLAST 索引
```

`sequences` 是筛查输入。它的 header 携带基因元数据：gapit 原生库使用带标签的 `gapit/v1`
形式（`gapit|db=...|gene=...|acc=...|func=...`，空格后是 product）；旧式 abricate 库使用
`~~~` 分隔的 header。gapit 两种都能读。

BLAST 索引（核酸 `.nin/.nsq/.nhr`，蛋白 `.pin/.psq/.phr`）是生成物，不是源数据。索引缺失或
过期时，`gapit setupdb` 或正常构建流程会重建。

`gapit-manifest.json` 最后写入，用于认证产物（`sha256`、`n_records`、`dbtype`、`kind`、
`source`、`typing_schema`、许可证、工具版本）。发现和列表会使用它，但为了兼容 abricate，
没有 manifest 的旧式目录仍然可以读取。

### Gene 数据库

gene 数据库每个已知基因/标记一条 FASTA 记录。contig 筛查使用 `blastn`（蛋白库使用
`blastx`）；reads 筛查使用 `minimap2`。gene 库可以携带可选 `typing.json`，让
`gapit typing` 从筛查表判定表型。

### Cluster 数据库

cluster 数据库（`kind: cluster`）每个**位点**一条 FASTA 记录，另有 `features.json` 描述
位点到基因及坐标的映射。`gapit screen` 使用 minimap2 cluster 引擎，根据 `cs` walk 计算基因/
位点评分。多 scheme 的 `typing.json` 输出 `gapit.cluster/2`；否则输出 `gapit.cluster/1`。
内置的 `vp_ok` 是带 typing 的 cluster 示例。

### 数据库来源

| 来源 | 含义 |
|---|---|
| `db fetch` | Provider 目录：抓取时从上游下载 |
| Bundled | 随 wheel 发布的内容，首次使用时物化 |
| `db build` | 自定义 FASTA、GBK 或 GFF/GBK cluster 数据 |
| `db install` | 经 SHA256 校验的本地字节，原子安装 |
| abricate 数据目录 | 旧式 `~~~` header，无需 manifest |

内置库物化后也是用户数据：它们被复制并建索引到 datadir，而不是直接从 wheel 执行。

## datadir 的位置

优先级从高到低：

1. `--datadir PATH` / `-D PATH`
2. `$GAPIT_DATADIR`
3. `$XDG_DATA_HOME/gapit/db`（仅当 `XDG_DATA_HOME` 是绝对路径）
4. `~/.local/share/gapit/db`

`gapit db path` 打印生效路径，且**不会创建**它。加 `--json` 输出 `gapit.datadir/1` 文档，
包含 `datadir`、`exists` 和 `source`（`cli`、`env`、`xdg`、`default`）——适合 Docker
entrypoint 和 agent 工作流。

不同命令对缺失目录的行为：

| 命令 | 缺失 datadir 时 |
|---|---|
| `gapit db list` | 创建空根目录，然后列出目录清单 |
| `gapit setupdb` | 创建根目录、物化内置库、为全部库建索引 |
| `gapit db fetch` / `db build` | 作为写路径的一部分创建根目录 |
| `screen` / `typing` 指向内置库名 | 创建根目录并物化该内置库 |
| `screen` / `typing` 指向已 fetch 的名字 | 要求目录存在（`DATADIR_NOT_FOUND`，退出码 4） |
| `db search` / `db outdated` | 要求目录存在（`DATADIR_NOT_FOUND`，退出码 4） |

`db list` 和 `setupdb` 是显式的首次使用引导路径；`db path` 只做路径查询。

## 为什么不默认放在 `$CONDA_PREFIX`？

conda 环境是**工具链环境**，不是用户数据仓库：

- 环境可能是只读或共享的；
- 环境经常被删除重建；
- 运行时往 `$CONDA_PREFIX` 写数据通常不被鼓励；
- 每个环境都会复制一份大体积参考库；
- 环境里安装的 wheel 版本，不等于运行时 datadir 里的数据库版本（内置库只物化一次）。

如果参考数据是**随 conda 包发布的只读、固定版本内容**，放 `$CONDA_PREFIX/share/...` 可以是
合理的。但 gapit 的数据库是运行时下载、构建、刷新和自定义的，属于用户状态，应该放在用户数据
目录。

如果确实要环境内数据，显式选择：

```bash
export GAPIT_DATADIR="$CONDA_PREFIX/db/gapit"
gapit db list
```

代价正好和默认相反：删除环境会删除数据库；每个环境一份；嵌套 shell 需要恢复或取消变量。

## 跨环境共享 datadir

共享是特性而非意外：只下载一次、只建一次索引，多个环境看到同一份参考内容。但它**不是**隔离
或安全边界。

### 安全的部分

- 对不可变数据库文件的并发**读**（`screen`、`typing`、`db list`）；
- 同一数据库同一时间只有一个写者；
- 不同的自定义数据库名可以共存。

### 需要小心

**版本漂移。** 环境 A 的 wheel 可能带比 datadir 更新的内置快照。`db list` 显示安装状态，但
`db outdated` 不会比较 wheel 内置快照与已安装副本。升级某个环境的 gapit 后，要么显式刷新
（`gapit db fetch <name> --force`），要么使用项目级 datadir。

**命名冲突。** 两个环境构建同名自定义库时，没有 `--force` 不会静默覆盖，但 `--force` 会影响
所有共享该 datadir 的环境。使用唯一名字或显式项目 datadir。

**写并发。** 并发 `db fetch`/`db build` 现在会通过数据库级隐藏 `flock` 串行化，并且先构建到
隐藏 staging 目录，再 rename 切换。读者通常看到旧库或新库；但 canonical 路径替换仍有一个很短
的 rename 窗口。如果共享 datadir 上的读者不能容忍这个窗口，不要并发写同一数据库。
`gapit setupdb` 会就地重建 BLAST 索引，因此不要对同一数据库并发运行它和
`db fetch`/`db build`。

**权限。** 新 datadir 和事务 staging 默认 `0700`。多用户共享需要显式的 group/ACL 策略；默认
是每用户私有。

**信任。** 部分 provider 仍然只有 HTTP，且下载没有 checksum pin。共享 datadir 会放大一次被污染
fetch 的影响。高保障场景应固定数据库版本/校验和，并优先使用只读共享参考目录。

## Docker

容器内的 `$HOME` 和 `$CONDA_PREFIX` 写入通常是**临时的**：它们位于容器可写层，容器删除后消失。
只读根文件系统会让默认路径直接失败。

使用显式挂载和环境变量：

```bash
docker volume create gapit-data

docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -v gapit-data:/data/gapit \
  -e GAPIT_DATADIR=/data/gapit \
  gapit:0.6.1 gapit db list
```

CI 中挂载或缓存 workspace 目录，并固定版本：

```bash
docker run --rm \
  -v "$PWD/.gapit-db:/data/gapit" \
  -e GAPIT_DATADIR=/data/gapit \
  gapit:0.6.1 gapit screen sample.fna --db ecoli_dec
```

可以为**只读**筛查把参考库烤进镜像：

```dockerfile
RUN gapit setupdb --datadir /opt/gapit/db
ENV GAPIT_DATADIR=/opt/gapit/db
```

运行时 `db fetch`/`db build` 的写入会落在容器层（随后消失），或者在只读根文件系统下失败。
如果容器需要写自定义库，挂载可写 volume；如果同时需要镜像内参考库，可以从镜像 seed volume，
或改用项目级 datadir。

## 安全与可复现性检查表

- 把 datadir 放在用户 home 或显式项目/volume 路径。
- entrypoint 用 `gapit db path --json`，不要猜路径。
- 已发表分析固定数据库版本，并记录 manifest。
- gapit 升级后有意识地刷新数据库。
- 没有 ACL 时不要跨用户共享可写 datadir。
- 把 provider 下载视为信任边界；尽量校验上游内容。
- 备份 `records.jsonl`/自定义输入，而不是生成的 BLAST 索引。

## 相关命令

```text
gapit db path [--json]      # 解析 datadir，不创建
gapit db list [--json]      # 目录清单 + 安装状态；会引导创建根目录
gapit setupdb               # 引导 + 物化 + 建索引
gapit db fetch NAME         # provider 下载/构建
gapit db build NAME INPUT   # 自定义 gene 或 cluster 数据库
gapit db outdated           # 已安装库年龄（严格；要求 datadir）
```
