# contig 筛查

`gapit screen` 用 BLAST 把 contig 文件比对到基因数据库并报告哪些基因存在，流水线与
命中规则与 abricate 相同。默认输出 abricate 兼容的 TSV；JSON 和 Markdown 是一等
替代格式（见 [./outputs.md](./outputs.md)）。

数据库的搭建见 [./databases.md](./databases.md)；完整的上手演练在
[./quickstart.md](./quickstart.md)。

## 选项

转写自 `gapit screen --help`（gapit 0.5.0）。标记为 *reads 模式* 的参数只在传入
`--r1`/`--r2` 时生效，记录在 [./reads.md](./reads.md)。

每个选项都接受 **短** 列给出的单横线短形式（例如 `--db` 的 `-d`）；长形式仍是规范写法，
布尔否定形式（如 `--no-merge-fragments`）保持只有长形式。`--db` 还单独带一个词风格的
短别名：`-db ncbi`（需空格分隔；连写 `-dbncbi` 仍会绑定到 `-d` 加值）。位置参数若是一个
已知数据库名但在磁盘上不存在，将以指明 `--db` 的用法错误被拒绝。

| 参数 | 短 | 类型 | 默认值 | 说明 |
|---|---|---|---|---|
| `FILE...` | — | 路径 | 必填* | 要筛查的 contig 输入文件。全 FASTQ 的通配符进入 reads 模式并按文件名自动配对样本（见 [./reads.md](./reads.md)）。 |
| `--db` | `-d`, `-db` | str | 必填 | 用于筛查的数据库（数据目录的子目录）。没有默认值：一次筛查绝不会悄悄物化内置数据库 —— 请显式选择（`gapit db list`）。 |
| `--datadir` | `-D` | path | `$GAPIT_DATADIR`，然后 `~/.local/share/gapit/db` | 数据库目录。 |
| `--minid` | `-i` | float | `80.0` | 最小一致性百分比，`0 < x <= 100`。在 BLAST 内部通过 `-perc_identity` 强制执行。 |
| `--mincov` | `-c` | float | `80.0` | 最小覆盖度百分比，`0 <= x <= 100`。对未取整的浮点值做后置过滤。 |
| `--threads` | `-t` | int | `1` | BLAST 工作线程数（传给 `-num_threads`）。 |
| `--jobs` | `-j` | int | `1` | 并发筛查 N 个输入文件。输出顺序始终是输入顺序。 |
| `--merge-fragments` | `-m` | 开关 | 关闭 | 合并被 contig 边界拆开的基因片段（gapit 扩展；见下文）。 |
| `--fofn` | `-F` | path | 无 | 文件名列表文件；取代位置参数 FILE。 |
| `--quiet` | `-q` | 开关 | 关闭 | 静默 stderr 诊断。 |
| `--noheader` | `-n` | 开关 | 关闭 | 不输出 `#FILE ...` 表头行。 |
| `--nopath` | `-p` | 开关 | 关闭 | FILE 列只保留文件名。 |
| `--debug` | `-v` | 开关 | 关闭 | 详细的 stderr 诊断；回显每条外部命令行。 |
| `--format` | `-f` | tsv\|csv\|json\|md | `tsv` | 输出格式（tsv 在所有路径上都是默认值：reads 模式每完成一个文件或样本就流出对应的表格行/分块）。json/md 是给 agent 的显式选项；json 在结束时一次性写出（单一文档）。 |
| `--output` | `-o` | path | stdout | 把报告写入 PATH 而不是 stdout（截断已存在的文件）。流式格式逐文件刷新；此时 stdout 不输出任何数据。 |
| `--aligner` | `-a` | blastn\|minimap2 | 按输入决定 | 比对引擎（默认：contig 文件用 blastn，`--r1`/`--r2` reads 用 minimap2）。`--aligner minimap2` 把位置参数给出的 FASTA 装配体送进 minimap2 引擎（要求 FASTA 内容；见 [./reads.md](./reads.md)）。 |
| `--r1` | `-1` | str | 无 | *reads 模式。* reads 或 assembly FASTA 文件，逗号分隔，每条 lane 一个。 |
| `--r2` | `-2` | str | 无 | *reads 模式。* 逗号分隔的 mate FASTQ 文件；数量必须与 `--r1` 一致。 |
| `--read-type` | `-x` | sr\|map-ont\|map-hifi | FASTQ 用 `sr`，FASTA 用 `map-ont` | *reads 模式。* minimap2 预设；省略时按检测到的输入解析。 |
| `--min-breadth` | `-b` | float | `90.0` | *reads 模式。* 判定存在的最小广度百分比。 |
| `--min-gene-cov` | `-g` | float | `90.0` | *仅基因簇数据库。* 基因判定为 `present` 的最小覆盖度百分比。 |
| `--min-gene-id` | `-G` | float | `90.0` | *仅基因簇数据库。* 基因判定为 `present` 的最小一致性百分比。 |
| `--min-cluster-cov` | `-C` | float | `96.0` | *仅基因簇数据库。* 最佳位点判定的最小位点覆盖度百分比。 |

\* 位置参数 FILE 或 `--fofn`，或经 `--r1` 进入 reads 模式。位置文件与 `--r1`/`--r2`
互斥。

## 输入文件

归一化是原生的（不依赖外部 `any2fasta`）：plain FASTA、gzip 和 bzip2 压缩的 FASTA、
FASTQ、GBK、EMBL 都可以。转换出的 FASTA 缓冲在内存里（基因组规模的装配体也只有几
MB），经 stdin 喂给 blastn。如果归一化失败（不是序列文件），gapit 在 stderr 打印
`gapit.error/1` 信封并以退出码 5 结束。

归一化之前有一条例外：当**每个**位置参数文件都是 FASTQ（`.fastq`/`.fq` ± `.gz`，
或扩展名含糊时按内容嗅探为 FASTQ）且未指定 `--aligner` 时，整条命令改走 reads 引
擎 —— `gapit screen -d ecoli_dec *.gz` 把 glob 按自动配对的样本筛查，而不是当作
contig。FASTA 与 FASTQ 位置参数混用是用法错误，错误信息点名 reads 文件。配对规则
与按样本输出见 [./reads.md](./reads.md)。

`--fofn` 文件每行列一个路径，完全取代位置参数。

## 实战示例

`tinyamr` 夹具数据库在仓库测试数据里（三个短 AMR 基因）。把它复制到临时数据目录、
建索引，再把 `$GAPIT_DATADIR` 指过去。实验时千万别写进
`~/.local/share/gapit/db`。

```console
$ mkdir -p /tmp/gapit-demo/datadir/tinyamr
$ cp tests/data/db/tinyamr/sequences /tmp/gapit-demo/datadir/tinyamr/sequences
$ gapit setupdb --datadir /tmp/gapit-demo/datadir
Indexed tinyamr (3 sequences, nucl)
$ export GAPIT_DATADIR=/tmp/gapit-demo/datadir
```

### 默认 TSV

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr
Processing: tests/data/contigs/full.fa
Found 1 genes in tests/data/contigs/full.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
```

`Processing:` 和 `Found N genes` 行走 stderr；表头和命中行走 stdout。丢弃 stderr 后
输出就是纯数据：

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr 2>/dev/null
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
```

### JSON

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr --format json
{
  "schema": "gapit.report/1",
  "tool": {
    "name": "gapit",
    "version": "0.5.0"
  },
  "created_at": "2026-09-19T01:12:04Z",
  "params": {
    "db": "tinyamr",
    "minid": 80.0,
    "mincov": 80.0,
    "threads": 1
  },
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

字段名是 snake_case，单位显式（`identity_pct`、`coverage_pct`）。schema 可自省：
`gapit schema report`。见 [./outputs.md](./outputs.md)。

### `--mincov` 如何过滤

`tests/data/contigs/partial.fa` 只携带 88 bp 的 `blaTEM-1` 参考的前 44 个碱基，覆盖度
是 50%。默认阈值会丢弃它：

```console
$ gapit screen tests/data/contigs/partial.fa --db tinyamr --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
```

降低阈值就能保留：

```console
$ gapit screen tests/data/contigs/partial.fa --db tinyamr --mincov 50 --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/partial.fa	contig1	1	44	+	blaTEM-1	1-44/88	========.......	0/0	50.00	100.00	tinyamr	J01749.1:1-861	class A beta-lactamase TEM-1	BETA-LACTAM
```

### 压缩输入、CSV、路径与表头控制

```console
$ gzip -c tests/data/contigs/gap.fa > /tmp/gapit-demo/gap.fa.gz
$ gapit screen /tmp/gapit-demo/gap.fa.gz --db tinyamr --quiet --nopath
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
gap.fa.gz	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
$ gapit screen tests/data/contigs/full.fa --db tinyamr --format csv --quiet --nopath
#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE
full.fa,contig1,1,79,+,tetA,1-79/79,===============,0/0,100.00,100.00,tinyamr,NC_000913.3:100-900,tetracycline efflux pump TetA,TETRACYCLINE
```

### 多文件：`--fofn` 与 `--jobs`

```console
$ printf '%s\n' tests/data/contigs/full.fa tests/data/contigs/gap.fa > /tmp/gapit-demo/files.txt
$ gapit screen --fofn /tmp/gapit-demo/files.txt --db tinyamr --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
tests/data/contigs/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
$ gapit screen tests/data/contigs/full.fa tests/data/contigs/gap.fa --db tinyamr --jobs 2 --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
tests/data/contigs/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
```

用 `--jobs 2` 时两个文件并发筛查，但 stdout 仍按输入顺序排列，输出与串行运行完全
一致。

### typed 基因数据库：判定是第二条命令

当 `--db` 指向带 `typing.json` 构建的基因数据库（捆绑的 `ecoli_dec`，或任何
`gapit db build --typing` 构建）时，筛查是纯基因检测：每种格式的输出都是逐字节相同
的冻结 15 列 abricate 表——typed 与 untyped 基因数据库在筛查面上无法区分。判定来自
两阶段流水线的第二阶段：用 `--output`（或重定向）写出结果表，再对它运行
[`gapit typing`](./typing.zh.md#两阶段判定工作流)——或者直接经管道送入，规范一行式：

```console
$ gapit screen dec_s3_stx2a_escV_aggR_uidA.fasta dec_s2_pic_astA_uidA.fasta --db ecoli_dec --output dec.tsv --nopath --quiet
$ gapit typing dec.tsv --quiet
FILE	SCHEME	PHENOTYPE	CONFIDENCE	SCORE	RUNNER_UP	NOTES	GENES
dec_s3_stx2a_escV_aggR_uidA.fasta	gb4789_6	EHEC	high	1.0000	EAEC (1.0000)	GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE); severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up	aggR;escV;stx2a;uidA
dec_s3_stx2a_escV_aggR_uidA.fasta	risk_monitoring	EHEC	high	1.0000	EAEC (1.0000)	GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE); severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up	aggR;escV;stx2a;uidA
dec_s2_pic_astA_uidA.fasta	gb4789_6	EAEC	high	1.0000	EHEC (0.0000)	GB 4789.6: any of aggR/pic/astA	astA;pic;uidA
dec_s2_pic_astA_uidA.fasta	risk_monitoring	non-DEC	low	0.0000	STEC (0.0000)	GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE); severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up	astA;pic;uidA
```

管道形式完全省去中间文件——screen 的 stdout 就是 typing 的 stdin（两侧都保持
stdout 纯净：数据走管道，诊断信息走 stderr）：

```console
$ gapit screen dec_s3_stx2a_escV_aggR_uidA.fasta --db ecoli_dec --nopath --quiet | gapit typing --quiet
```

只要 stdin 不是终端，`gapit typing` 就从中读取；`gapit typing -` 是显式 stdin
标记（挂着终端也有效——读取直到 EOF）。完整的输入契约见
[typing 页面](./typing.zh.md#两阶段判定工作流)。

同一条管道也能喂给矩阵视图。一次筛查多个 assembly 会把所有行汇入一张表（每行
携带同一个 `DATABASE`；`FILE` 列横跨整个批次），因此
[`gapit summary`](./summary.zh.md) 可以直接从管道读取并折叠成基因×文件存在矩阵：

```console
$ gapit screen -d ecoli_dec *.fna --quiet | gapit summary
```

基因簇数据库是例外：其判定集成在筛查本身（typed 基因簇 TSV 的 `PHENOTYPE` 列，
见[下文](#基因簇数据库-kind-cluster)）。

## 流式输出与 `--output`

长时间的批量筛查可以边跑边反馈，报告也可以直接落盘：

- **tsv/csv/md 按文件流式输出。** 表头（或 Markdown frontmatter）在筛查开始时打印
  一次，每个文件的行/小节在该文件完成的那一刻打印 —— 顺序始终是输入顺序。
  `--jobs N > 1` 时采用队头阻塞式发射：第 *i* 个文件的输出要等文件 1..*i* 全部完成
  才流出（线程池按输入顺序产出），因此 stdout 字节与串行运行完全一致。位置参数
  FASTQ 通配符路径（`screen -d db *.fastq.gz`，见 [./reads.md](./reads.md)）按样本遵
  循同一契约：`--jobs` 并行各样本，md 先流出静态 frontmatter、每完成一个样本即流
  出其 `## <样本>` 小节，批量中途失败时 `--output` 保留已流出的前缀。
- **md 流式输出；json 在结束时一次性写出。** Markdown frontmatter 是**静态**元数据
  （schema、tool、`created_at`、db、阈值 —— 不含运行总数），因此可以先行输出，每个
  文件的小节随文件完成即时流出，与 tsv 行完全一致。基因簇数据库的每个文件小节自带
  该文件的摘要行和基因表。运行总数（`files`、`hits` 等）由 JSON 文档承载 —— json 是
  单一对象，在结束时一次性写出。
- **`--output PATH`**（所有引擎，含 reads 模式）把报告写入 PATH 而不是 stdout：文
  件在第一个输出字节时打开（截断已存在的文件 —— v1 为覆盖语义，不追加），每个流
  式分块都刷新落盘，此时 stdout **不输出任何数据**（stderr 诊断不变）。
- **批内出错。** 若文件 1..*k-1* 已流出后第 *k* 个文件失败，已发射的输出原样保留
  —— stdout 上已打印；`--output` 文件里保留表头/frontmatter 和文件 1..*k-1* ——
  然后在 stderr 打印类型化的 `gapit.error/1` 信封并按文档退出码退出。在任何输出产
  生之前就失败的运行（用法、依赖、数据库错误）根本不会创建 `--output` 文件。

### `--output`

```console
$ gapit screen tests/data/contigs/full.fa tests/data/contigs/gap.fa --db tinyamr --output report.tsv
Processing: tests/data/contigs/full.fa
Found 1 genes in tests/data/contigs/full.fa
Processing: tests/data/contigs/gap.fa
Found 1 genes in tests/data/contigs/gap.fa
$ cat report.tsv
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
tests/data/contigs/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
```

`Processing:`/`Found` 行走 stderr；stdout 为空，`report.tsv` 的内容与同一运行打印
的字节完全相同。

### `--debug`

把原生归一化步骤和每条外部命令行回显到 stderr：

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr --debug 2>&1 >/dev/null
Processing: tests/data/contigs/full.fa
gapit: normalize: tests/data/contigs/full.fa (fasta)
gapit: run: blastn -task blastn -dust no -perc_identity 80.0 -db /tmp/gapit-demo/datadir/tinyamr/sequences -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen stitle' -num_threads 1 -evalue 1E-20 -culling_limit 1 -max_target_seqs 10000
Found 1 genes in tests/data/contigs/full.fa
```

`--debug` 不改变 stdout；它是纯 stderr 参数。

## 片段合并（`--merge-fragments`）

草图装配体常把一个基因切断在 contig 边界上：每个 contig 只携带基因的部分片段，每
个片段单独都过不了 `--mincov 80`，基因于是漏报。`--merge-fragments`（gapit 扩展，
默认关闭；仅限 blastn contig 模式：reads 模式或 `--aligner minimap2` 下是用法错
误）正好挽救这类基因：

- 常规的单命中过滤之后，gapit 把某个基因的低于 `--mincov` 的片段归组：这些片段已在
  `(contig, start, end)` 上去重，并通过了 `--minid`。
- 如果片段的受试区间**并集**覆盖了基因的 `>= mincov`（`100 * union_length / slen`，
  按未取整浮点比较，边界规则与单命中过滤相同），**且**至少涉及两个片段，就报告一
  条合并命中。
- 已经有一条单独及格命中的基因永远不会被合并，它的零散部分会被忽略。并集仍低于
  `mincov` 的组依旧什么都不报。

合并行的语义：

- `%COVERAGE` 是并集覆盖度；`%IDENTITY` 是按比对长度加权的片段一致性均值；`GAPS`
  和 `COVERAGE_MAP` 的 `/` 标记来自求和后的缺口数，map 按区间并集分箱（与普通行相
  同的 15 字符算术）。
- `SEQUENCE` 是携带这些片段的 contig 名的逗号串接（已排序；TSV 没问题，合并行请避
  免 `--format csv`，字段内逗号会有歧义）。
- `START`/`END`/`STRAND` 取自**锚点**片段：比对长度最大的那个（平局取 `(contig, start)`
  顺序的第一个）。`COVERAGE` 是并集的包围跨度。
- JSON（`gapit.report/1`）额外携带 `"merged": true` 和一个 `fragments` 数组，每个片
  段一项 `{contig, start, end, strand, identity_pct, coverage_pct}`：可选的累加字
  段，非合并行一律不出现。Markdown 会在表下为每个合并基因加一行细节。

不带这个参数时整条路径都不会运行：默认输出与 abricate 逐字节一致（下文的注意事项
仍然全部适用）。

## 基因簇数据库（`kind: cluster`）

当 `--db` 指向一个基因簇数据库时（经 [`gapit db build`](./databases.md) 从 GBK/GFF
构建，或从 kaptive 提供商抓取），`gapit screen` 会分发到**基因簇引擎**而不是
blastn 基因流水线。基因路径没有任何变化：cluster 分发是累加式的，gene 类筛查与
cluster 功能引入前逐字节一致。

工作方式：

- 每个输入文件一次 **minimap2 `asm20`** 调用（query 是文件的 contig，target 是数据
  库的位点 FASTA，`--cs` 短格式比对轨迹）。
- 每个**位点**的并集覆盖度与一致性来自跨全部主比对的 cs walk：来自多条 contig 的
  记录做**并集**，所以被 contig 边界打碎的位点仍能筛查为完整存在。
- 每个注释的**基因**得到一个判定：`present`（覆盖度 ≥ `--min-gene-cov` 且一致性 ≥
  `--min-gene-id`）、`partial`（覆盖度 ≥ 50% 但未达 present 阈值）或 `absent`。
- 位点按覆盖度、一致性、覆盖碱基数、id 排序；排名第一的位点在覆盖度达到
  `--min-cluster-cov`（默认 96，kaptive 风格的置信下限）时成为该文件的最佳判定。
  低于下限则不做判定。

守卫（用法错误，退出码 2）：`--minid`/`--mincov`/`--merge-fragments`/`--jobs`/`--aligner`
是基因引擎参数，在基因簇数据库上会被拒绝；三个 `--min-gene-*`/
`--min-cluster-cov` 参数是 cluster 专属，在其他地方一律拒绝，包括 reads 模式
（v1 中基因簇筛查只接受 assembly FASTA）。

### 表型判定（`typing.json`）

基因簇数据库可以携带一份 `gapit.typing/1` 评分规范（用 `gapit db build --typing
FILE` 安装）。有它时，每个文件的最佳判定都会注释上一个**表型**：

- 各规则独立打分：`weighted_genes`（基因存在性，带一致性下限和可选负标记）、
  `cluster_match`（位点覆盖度/一致性/关键基因的加权分量，各分量带下限）或
  `learned_linear`（在命名特征上训练出的 sigmoid 模型）；决策层再应用文档的
  `cutoff` 和 `ambiguity_margin`：胜者明确时以 `high` 置信度给出判定，两条规则落在
  容差内则给出 `ambiguous` 判定（表型为 null，两个候选都列出），最佳分低于 cutoff
  时回退到文档的 fallback 字符串。
- JSON（`gapit.cluster/1`）携带 `best.phenotype` 以及一个累加的
  `best.phenotype_detail` 块：`{score, confidence, components[], runner_up,
  ambiguous[]}`，即逐规则、逐分量的可解释分解。
- TSV/CSV 表头在 TYPE 之后多出 **PHENOTYPE** 列（歧义或未判定时为 `-`）；无 typing
  的基因簇数据库表头保持逐字节不变。
- 规则引用了数据库没有的基因/位点时，会在开头失败并给出 `TYPING_UNKNOWN_GENE`
  类型化错误（`db build --typing` 时同样强制检查）。

目前无 typing 的基因簇数据库（各 kaptive 抓取）只报告位点判定，即 kaptive 风格的
输出，`phenotype` 为 null。

## 注意事项

- **覆盖度过滤基于未取整浮点。** `%COVERAGE = 100 * (length - gaps) / slen` 在取整
  前与 `--mincov` 比较，显示用 `%.2f`。79.996% 的命中会打印成 `80.00` 但过不了默认
  `--mincov 80`，从而被丢弃。这是刻意的 abricate 一致性，不是取整 bug。
- **一致性从不过滤。** `%IDENTITY` 就是 BLAST 的 `pident`，按 `%.2f` 打印。
  `--minid` 在 blastn 内部经 `-perc_identity` 强制执行。
- **去重，不合并。** 拥有相同 `(contig, start, end)` query 跨度的 BLAST 行会折叠，
  第一行胜出。*不同*跨度的重叠命中全部报告；gapit 不合并重叠区间。这与 abricate 完
  全一致，是特性而非缺陷。（上文可选择的 `--merge-fragments` 是唯一经认可的例外；
  默认路径永不合并。）
- **蛋白质数据库。** `dbtype prot` 的数据库（例如 `bacmet2`）经 `blastx` 筛查，而
  blastx 不接受 `-perc_identity`。gapit 会在 stderr 打印
  `--minid is not applied to protein databases (abricate parity)` 并继续运行。
- **排序。** 单个文件内，行按 SEQUENCE 再按 START 稳定排序。一个文件的行在该文件完
  成后输出（或刷新到 `--output`）；即使 `--jobs > 1`，文件仍按输入顺序输出（队头
  阻塞：第 i 个文件等文件 1..i）。
- **退出码。** 2 用法错误，3 缺依赖，4 数据库错误，5 输入错误，1 意外错误。失败时
  在 stderr 打印 `gapit.error/1` JSON 信封；见 [./outputs.md](./outputs.md)。
