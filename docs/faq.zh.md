# 常见问题

对实践中常出现的问题给出简短回答。每个论断都有源码或真实运行支撑。

## gapit 是 abricate 的即插即用替代品吗？

对 contig 筛查，是的：同一 BLAST 流水线、同样的命中规则、stdout 输出 abricate 格式
TSV，并在一个语料上对真实 abricate 做过逐字节一致性校验。它能原样读取 abricate 格
式的数据库（遗留 `~~~` 表头）。反过来不成立：abricate 读不了 gapit 原生数据库
（`gapit/v1` 标记表头）。细节：`./databases.md`。

## 显示 80.00% 覆盖度的命中为什么被过滤了？

`--mincov 80` 的比较在未取整浮点上运行；TSV 列是显示取整。真实案例：1000 nt 基因
在 25000 nt 参考上有 1 nt 缺口，得到 100*(20000-1)/25000 = 79.996%，打印成 `80.00`
但被丢弃。这与 abricate 完全一致。背景见 `./screen.md`。

## 为什么我得到零命中？

检查阈值。默认 `--minid 80`（在 blastn 内部经 `-perc_identity` 强制）和
`--mincov 80`（对未取整覆盖度做后置过滤，见上一个问题）。分歧大或部分的基因请调
低。注意零匹配的输入不是错误：表头行照常打印，退出码为 0。

```console
$ gapit screen none.fa --db tinyamr
Processing: none.fa
Found 0 genes in none.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
$ echo $?
0
```

## 如何更新数据库？

用 `--force` 覆盖旧库重新抓取：

```console
$ gapit db fetch ncbi --force
```

不带 `--force`，抓进已存在的数据库会失败而不是覆盖。

## 我的数据目录在哪里？

解析顺序：命令行 `--datadir`，然后 `$GAPIT_DATADIR`，然后
`~/.local/share/gapit/db`。解析出的路径不存在是错误（`DATADIR_NOT_FOUND`，退出码
4）。见 `./databases.md`。

## 能用现有的 abricate 数据库吗？

能，原样使用：把 gapit 指向 abricate 建好的数据目录，直接筛查。BLAST 索引缺失或
过期时，全部重建：

```console
$ gapit setupdb
```

只有 gapit 原生输出是单向的（abricate 读不了 `gapit/v1` 数据库）。

## 如何筛查基因簇或 Kaptive 位点？

基因簇数据库是第二种数据库，用 `gapit db build NAME loci.gbk` 从 GenBank/GFF 构建
（可选 `--typing FILE` 获得表型判定），或从七个 kaptive 提供商（`kpsc_k`、`kpsc_o`、
`kosc_k`、`kosc_o`、`ab_k`、`ab_o`、`ecoli_kps`，即 Kaptive v3 官方安装关键字）抓取。
`gapit screen assembly.fa --db kpsc_k` 会分发到
minimap2 基因簇引擎，每文件报告一个最佳位点判定，附逐基因判定，typed 数据库上还有
表型。细节：`./screen.md` 和 `./databases.md`。

## `--threads` 和 `--jobs` 有什么区别？

`--threads` 是一次筛查运行内部的 BLAST 工作线程（传给 `-num_threads`，默认 1）。
`--jobs` 是并行筛查的输入文件数（默认 1）。N 个单 contig 文件，`--jobs N` 能扩展；
一个巨型文件，靠 `--threads`。

## reads 模式为什么拒绝 `--format tsv`？

reads 结果没有逐命中的 TSV 语义（它们是逐基因的广度/深度行），所以 TSV/CSV 是用法
错误。请用 `json`（默认）或 `md`：

```console
$ gapit screen --r1 reads.fq --read-type sr --db tinyamr --format tsv
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"--format tsv|csv is not available in reads mode (use json or md)","context":{}}
$ echo $?
2
```

reads 模式：`./reads.md`。

## 这个 "MISSING_DEPENDENCY" 错误是什么？

某个外部二进制（blastn、makeblastdb、blastdbcmd、minimap2）不在 PATH 上。gapit 以
退出码 3 结束，信封指明是哪个二进制：

```text
{"schema":"gapit.error/1","code":"MISSING_DEPENDENCY","message":"required binary not found on PATH: blastn","context":{"binary":"blastn"}}
```

安装 BLAST+ 及其伙伴（pixi 会替你装好），并确保它们在 PATH 上。

## 数据库适用哪些许可证？

gapit 本身以 MIT 许可。包里不内置任何数据库内容：每个提供商都在抓取时从上游下载，
内容保留其原始许可证（NCBI 公有领域、CARD 的 McMaster 非商业条款、VFDB 的
CC BY-NC、CGE、Kaptive 的 GPL-3.0 等），gapit 不重新授权。声明了许可证的提供商会在
`gapit db list --json` 中暴露它。细节见 SPEC.md §9。

## 如何向默认集合添加新数据库？

写一个提供商模块并按名字抓取：每个提供商都在抓取时从上游下载、变换记录、本地构建
—— 绝不内置任何内容（CARD、VFDB 等上游许可禁止随 MIT 许可的发行版再分发）。表头格
式、变换和构建流水线见 SPEC.md §11。
