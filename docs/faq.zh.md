# 常见问题

对实践中常出现的问题给出简短回答。每个论断都有源码或真实运行支撑。

## gapit 是 abricate 的即插即用替代品吗？

contig 筛查刻意保留 abricate 的界面：同一 BLAST 流水线、同样的命中规则、stdout 输出
abricate 格式 TSV，并在一个语料上对真实 abricate 做过逐字节校验。从 v0.5.0 起，这条
基线被冻结为历史基线面，而不再是发布门槛——gapit 走自己的契约（原生
JSON/Markdown、reads 与基因簇引擎、分型），parity 校验（`pixi run -e parity parity`）
只作为冻结面的回归参照。它能原样读取 abricate 格
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

## 如何判定致病型或血清群？

基因数据库上的判定是两阶段的：先用 typed 数据库筛查，再把结果表读回来。内置的
`ecoli_dec` 面板携带双方案 DEC 判定：`gapit screen -d ecoli_dec sample.fna -o
result.tsv`，然后 `gapit typing result.tsv`——或者直接用管道把一步送进另一步。基因
簇数据库是例外：其判定集成在筛查本身（`PHENOTYPE` 列）。细节：`./typing.md`。

## `--threads` 和 `--jobs` 有什么区别？

`--threads` 是一次筛查运行内部的 BLAST 工作线程（传给 `-num_threads`，默认 1）。
`--jobs` 是并行筛查的输入文件数（默认 1）。N 个单 contig 文件，`--jobs N` 能扩展；
一个巨型文件，靠 `--threads`。

## reads 模式输出什么格式？

默认 TSV —— tsv 是 gapit 所有界面上的人类默认值，json/md 是给智能体的显式选项。
reads 结果是逐基因的广度/深度行，流式表格把样本名放在第一列 `#SAMPLE`；
`--format json`（或 `md`）则给出版本化的 `gapit.reads/1` 文档：

```console
$ gapit screen --r1 tests/data/reads/tetx_full.fq --read-type sr --db tinyreads
#SAMPLE	GENE	BREADTH%	DEPTH	READS	PRESENT	DATABASE	ACCESSION	PRODUCT
tests/data/reads/tetx_full.fq	tetX	100.00	2.30	12	yes	tinyreads	SYN-001	extended resistance determinant tetX
$ gapit screen --r1 tests/data/reads/tetx_full.fq --read-type sr --db tinyreads --format json
{"schema": "gapit.reads/1", ...}
```

reads 模式：`./reads.md`。

## 能用通配符一次筛查多个 FASTQ 文件吗？

能：当所有位置参数文件都是 FASTQ 时，`gapit screen` 进入 reads 模式，并按文件名
（`_R1`/`_R2`、`_1`/`_2` 约定）自动配对样本，所以 `gapit screen -d ecoli_dec *.fq.gz
-j 4` 一条命令就能并发筛查全部样本。FASTA 与 FASTQ 位置参数混用是用法错误；配对
约定与按样本输出见 `./reads.md`。

## 什么是逐基因一致性下限（floors）？

数据库侧的 `floors.json`（`gapit.floors/1`）为 reads 模式的存在判定声明每个基因的
最小比对一致性：低于下限的比对在广度/深度聚合之前被丢弃，于是"接近但不完全相同"
的同源基因不再过量判定。催生它的真实案例：内置 `ecoli_dec` 面板中的 SPATE 同源基因
`pic` 在约 86.6% 一致性、97.6% 广度时越过了 90% 存在阈值，而 blastn contig 路径
正确地拒绝了它；随库发布的下限（`pic` 90）让 reads 模式与之一致。用 `gapit db
build --floors FILE` 安装，或直接把文件放进数据库目录；没有该文件的数据库筛查输出
与之前逐字节一致。细节：`./reads.md` 和 `./databases.md`。

## 这个 "MISSING_DEPENDENCY" 错误是什么？

某个外部二进制（blastn、makeblastdb、blastdbcmd、minimap2）不在 PATH 上。gapit 以
退出码 3 结束，信封指明是哪个二进制：

```text
{"schema":"gapit.error/1","code":"MISSING_DEPENDENCY","message":"required binary not found on PATH: blastn","context":{"binary":"blastn"}}
```

安装 BLAST+ 及其伙伴（pixi 会替你装好），并确保它们在 PATH 上。

## 数据库适用哪些许可证？

gapit 本身以 MIT 许可。六个经过审计、许可宽松的内置库随 wheel 发布（公有领域的
`ecoli_dec`、`lm_doumith` 面板，加上 `ncbi`、`resfinder`、`ecoh`、`upec_expec_vf`
快照：Apache-2.0 / BSD-3-Clause / MIT），首次使用时物化进数据目录；其余每个提供商
都在抓取时从上游下载，内容保留其原始许可证（CARD 的 McMaster 非商业条款、VFDB 的
CC BY-NC、CGE、Kaptive 的 GPL-3.0 等），gapit 不重新授权。声明了许可证的提供商会在
`gapit db list --json` 中暴露它。细节见 SPEC.md §9 与数据库指南的
[内置数据库一节](./databases.zh.md)。

## 如何向默认集合添加新数据库？

写一个提供商模块并按名字抓取：每个提供商都在抓取时从上游下载、变换记录、本地构建。
除六个经审计的内置库之外，wheel 不携带任何内容（CARD、VFDB 等上游许可禁止随 MIT
许可的发行版再分发）。表头格式、变换和构建流水线见 SPEC.md §11。
