# 从 gspec 迁移（0.8 之前）

`gspec` 就是改名后的 `gapit` 0.7.x，功能完全相同。改名的原因是：在文献与代码检索中，
"gspec" 与 **GAPIT**（全基因组关联分析 R 包）完全撞名——搜这个工具，出来的全是 GWAS。

数据库与输出契约没有任何变化。本页逐项列出：哪些变了、哪些作为别名保留、哪些冻结不动。

## 变化的部分

| 面 | 0.7.3 及以前 | 从 0.8.1 起 |
|---|---|---|
| PyPI 发行名 | `gapit` | `gspec`（`gapit` 变为依赖 `gspec` 的过渡壳） |
| 命令 | `gapit` | `gspec`——`gapit` 命令仍然随包安装并可用 |
| 版本横幅 | `gspec 0.7.3` | `gspec 0.8.1` |
| `--version --json` 的 `name` | `gapit` | `gspec`（`gapit.version/1` 的 schema id 不变） |
| 数据目录环境变量 | `GAPIT_DATADIR` | `GSPEC_DATADIR`；旧名仍会读取 |
| 引擎环境变量 | `GAPIT_ENGINE` | `GSPEC_ENGINE`；旧名仍会读取 |
| 默认数据目录 | `~/.local/share/gapit/db` | `~/.local/share/gspec/db`；若机器上只存在旧目录则继续沿用旧目录 |
| 文档、README、CITATION | `gapit` | `gspec` |

## 刻意不改的部分

以下都是冻结契约——改名会破坏已有数据库、归档结果和流水线，因此保留 0.8 之前的拼写：

- **JSON schema id**：`gapit.report/1`、`gapit.typing/1|2|3`、
  `gapit.typing_result/1`、`gapit.cluster/1|2`、`gapit.error/1` 等全部
  `gspec.<document>/<version>` 名称。
- **数据库清单文件名**：`gapit-manifest.json`。
- **原生 FASTA 头部前缀**：`gapit|db=…|gene=…|acc=…|func=…`，以及清单里的
  `header_format: gapit/v1`。
- **Python import 包名**：`import gspec`（发行名是 `gspec`，两者本来就相互独立）。
- **已有数据库目录**：磁盘上不改动，继续可读。

## 升级

```bash
pip install --upgrade gspec        # 或：pip install --upgrade gapit
gspec db path                      # 确认当前使用的数据目录
gspec db list                      # 已有数据库照常解析
```

两个命令解析同一个数据目录，所以仍在调用 `gapit` 的流水线可以继续运行，脚本可以慢慢迁移到 `gspec`。

想显式换到新目录：

```bash
mv ~/.local/share/gapit ~/.local/share/gspec
```

## 时间线

- **0.8.x**：两个名字并存；旧环境变量继续生效；schema 冻结。
- **1.0**：计划在 `gspec.*` 之外**并行**发布 `gspec.*` schema id（是别名而非替换），
  并保留 `gapit` 命令与 import 作为兼容层。任何移除都需要主版本号 + 弃用周期。

## 引用

请引用 `gspec`（见 `CITATION.cff`）。使用 0.7.x 完成的工作，旧名指代的是同一软件；
0.7.3 的 Zenodo 记录仍可引用，0.8.1 的记录会声明其为前一版本。
