# 从 gapit 迁移

`gspec` 就是改名后的 `gapit`。改名原因是：在文献与代码检索中，"gapit" 与 **GAPIT**
（全基因组关联分析 R 包）完全撞名——搜这个工具，出来的全是 GWAS。

改名分两步：

| 版本 | 改了什么 | 没改什么 |
|---|---|---|
| **0.8** | 项目名、PyPI 发行名、CLI 横幅、文档、内置库 vendor 字符串 | `gapit` 命令、`import gapit`、所有磁盘上的名字 |
| **0.9** | 磁盘上的全部契约（schema id、manifest 文件名、FASTA 头前缀） | 旧拼写仍可读取 |

## 变化对照

| 面 | 0.8.1 及以前 | 从 0.9.0 起 |
|---|---|---|
| PyPI 发行名 | `gapit` | `gspec`（`gapit` 兼容包依赖它） |
| 命令 | `gapit` | `gspec`；`gapit` 仍可通过兼容包使用 |
| 版本横幅 | `gapit 0.8.1` | `gspec 0.9.0` |
| `--version --json` | `gapit.version/1`、`name: gapit` | `gspec.version/1`、`name: gspec` |
| JSON schema id | `gapit.<doc>/<n>` | `gspec.<doc>/<n>` |
| manifest 文件名 | `gapit-manifest.json` | `gspec-manifest.json` |
| FASTA 头前缀 | `gapit\|db=…` | `gspec\|db=…` |
| manifest `header_format` | `gapit/v1` | `gspec/v1` |
| 数据目录环境变量 | `GAPIT_DATADIR` | `GSPEC_DATADIR`（旧名仍读取） |
| 引擎环境变量 | `GAPIT_ENGINE` | `GSPEC_ENGINE`（旧名仍读取） |
| 默认数据目录 | `~/.local/share/gapit/db` | `~/.local/share/gspec/db`；只存在旧目录时继续用旧目录 |
| Python import | `import gapit` | `import gspec`；`import gapit` 由兼容包提供 |

## 仍然不变的部分

0.9 之前的拼写**全部可读**，写出时统一用新拼写：

- **schema id**：文档里的 `"schema": "gapit.typing/2"`、`gapit.report/1`、`gapit.cluster/2`
  等仍能通过校验，并被统一成 `gspec.*` 输出。归档 JSON 不需要重写。
- **manifest 文件名**：没有 `gspec-manifest.json` 时会读取 `gapit-manifest.json`；
  0.8 时代安装的数据库原样可用。
- **FASTA 头前缀**：`gapit|db=…|gene=…|acc=…|func=…` 仍可解码，新建库写出 `gspec|…`。
- **`header_format`**：`gapit/v1` 读作 `gspec/v1`。
- **环境**：`GSPEC_*` 未设置时读取 `GAPIT_*`；新目录不存在时沿用旧的 home 数据目录。

不需要强制迁移。想让现有数据目录统一成新名字，可重建（`gspec db setup --all` 原地重建索引；
`gspec db build` 会重写 manifest 与头前缀），或直接移动目录：

```bash
mv ~/.local/share/gapit ~/.local/share/gspec
```

## 升级

```bash
pip install --upgrade gspec        # 只装新名
pip install --upgrade gapit        # 保留 import gapit 与 gapit 命令
gspec db path                      # 确认当前数据目录
gspec db list                      # 已有数据库照常解析
```

两个命令解析同一个数据目录、同一批数据库，所以仍在调用 `gapit` 的流水线可以继续运行，
脚本可以按自己的节奏迁移到 `gspec`。

## 时间线

- **0.8**：品牌改名，磁盘契约不变。
- **0.9**：磁盘契约改名，旧拼写仍可读。
- **1.0**：计划保留旧的读取别名；任何移除都需要主版本号 + 弃用周期。

## 引用

请引用 `gspec`（见 `CITATION.cff`）。使用 0.x 完成的工作引用当时使用的版本即可，
Zenodo 记录指向的是同一软件。
