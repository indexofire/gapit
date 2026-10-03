# 安装 gapit

两条安装路径：PyPI wheel（外部二进制由你提供），或者用 [pixi](https://pixi.sh) 从
git 克隆安装（二进制由它代管）。

## 从 PyPI 安装

```bash
pip install gapit
```

wheel 里打包了 Python 包和 `gapit`/`gapit-mcp` console script，外加六个经过审计、
许可宽松的数据库内置包，首次使用时物化进数据目录（见数据库指南的
[内置数据库一节](./databases.zh.md)）；其余每个提供商都在抓取时从
上游下载，因为多家上游许可禁止再分发。它也**不**打包
BLAST+ 和 minimap2：请先装好它们（见下文 "外部二进制" 一节，例如
`conda create -n gapit-env -c bioconda blast minimap2`）。

## 从源码安装（pixi）

### 前置条件

- [git](https://git-scm.com) 和 [pixi](https://pixi.sh)。macOS/Linux 上：
  `curl -fsSL https://pixi.sh/install.sh | bash`
- 如果在 pixi 之外用 pip/pyproject 安装本包，需要 Python 3.11+。用 pixi 则无需操心：
  环境自带 Python（开发环境锁定 3.14）。

### 安装

```bash
git clone https://github.com/indexofire/gapit.git
cd gapit
pixi install
```

通过 pixi 运行 CLI，或者把环境的 bin 目录加进 PATH：

```bash
pixi run gapit --version
# 或者，等价地：
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
gapit --version
```

## 验证

```console
$ gapit --version
gapit 0.5.3
$ gapit --version --json
{"schema":"gapit.version/1","name":"gapit","version":"0.5.3"}
```

## 外部二进制

gapit 会调用 BLAST+（`blastn`、`blastx`、`makeblastdb`、`blastdbcmd`）和 `minimap2`。
两者都来自 pixi 环境（conda-forge 和 bioconda），无需手动安装。如果你在 pixi 之外运行
gapit，请自行确保这些二进制在 PATH 上。输入归一化（plain/gz/bz2 的 FASTA、FASTQ、
GenBank、EMBL）是原生的，不需要 `any2fasta`。

## 数据库引导

筛查需要数据目录（`$GAPIT_DATADIR`，然后是 `~/.local/share/gapit/db`；可用
`--datadir` 逐次覆盖）里至少有一个数据库。六个经过审计、许可宽松的数据库随 wheel
发布且无需网络：对某个内置库名的第一次筛查会自动物化它，`gapit setupdb` 则一次性
物化全部六个。其余每个提供商都在抓取时从上游下载（需要网络），因为多家上游许可
——CARD 的 McMaster 条款、VFDB 的 CC BY-NC、Kaptive 的 GPL-3.0——禁止随 MIT 许可的
包再分发。完整的数据库表与内置库说明见[数据库](./databases.zh.md)。

```bash
gapit db fetch all         # installs the default set (card, vfdb) into the default datadir
```

同样的安装动作放进一个临时数据目录，看看 fetch 会打印什么。进度行走 stderr，每个
数据库一行 JSON 回执走 stdout：

```console
$ gapit db fetch all --datadir /tmp/gapit-docs/dd
gapit: downloaded 1 source file(s)
gapit: read 6059 records from card
gapit: generated /tmp/gapit-docs/dd/card/sequences
gapit: self-check passed for card
gapit: BLAST index built (nucl)
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/gapit-docs/dd/card"}
gapit: downloaded 1 source file(s)
gapit: read 4769 records from vfdb
gapit: generated /tmp/gapit-docs/dd/vfdb/sequences
gapit: self-check passed for vfdb
gapit: BLAST index built (nucl)
{"db":"vfdb","records":4769,"dbtype":"nucl","destination":"/tmp/gapit-docs/dd/vfdb"}
```

随时用 `gapit db list`（或 `gapit db list --json`）确认数据库目录和安装状态。

## Shell 补全

```bash
gapit --install-completion   # bash, zsh, or fish; installs for the current shell
gapit --show-completion      # print the completion script to copy or customize
```

## 升级

```bash
git pull
pixi install
```

数据库内容不会随代码更新：重新运行 `gapit db fetch <name> --force`，从上游重建
已安装的数据库。
