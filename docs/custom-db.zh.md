# 自定义数据库

`gapit db build` 把任意参考基因 FASTA 一条命令变成完整构建的 gapit 原生数据库：
`records.jsonl`、`sequences` 投影、BLAST 索引、manifest。本页是一次带真实会话输出
的实战巡览。提供商、`db fetch` 和数据目录布局见[数据库](./databases.md)。

## 什么时候需要它

- 你要筛查没有提供商覆盖的基因：私有标记面板、实验室自筹的变异。
- 你想要一个只含上报基因的小数据库，让输出保持简短。
- 你手里有 abricate 时代的 `~~~` FASTA，想把它变成原生 gapit 数据库。

## 命令参考

```text
gapit db build NAME FASTA [OPTIONS]
```

| 选项 | 值 | 含义 |
|---|---|---|
| `NAME`（位置参数） | 文本 | 目标数据库名（创建在数据目录下）。 |
| `FASTA`（位置参数） | path | 输入 FASTA：plain、abricate `~~~` 或 `gapit\|` 表头，逐记录检测（接受 `.gz`/`.bz2`）。 |
| `--tsv` | path | 元数据 TSV：带 `gene`/`accession`/`function` 列的表头行。 |
| `--dbtype` | `nucl` 或 `prot` | 强制分子类型（默认：从序列自动检测）。 |
| `--datadir` | path | 数据库目录（默认：`$GAPIT_DATADIR`，然后 `~/.local/share/gapit/db`）。 |
| `--description` | 文本 | FASTA 表头没有描述文本的记录的默认产物。 |
| `--force` | 开关 | 数据库已存在时覆盖。 |
| `--quiet` | 开关 | 静默 stderr 诊断。 |

构建进度行走 stderr，一行 JSON 回执走 stdout。回执携带 `db`、`records`、`dbtype`
和 `destination`，字段与 `db fetch` 相同。下文的筛查块只展示 stdout；`gapit screen`
把 `Processing:` 行打印到 stderr（见[筛查](./screen.md)）。

下面的每个示例都在同一个临时目录里的一个会话中执行，所以回执里引用的是它的绝对路
径。要复现，请从仓库根目录运行：

```bash
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
mkdir -p /tmp/opencode/customdb-docs && cd /tmp/opencode/customdb-docs
```

输入基因是合成的，用 Python 的 `random` 模块生成（固定种子让本页可复现；永远不要
用不能公开的真实序列）：

```python
import random

rng = random.Random(11)
g1 = "".join(rng.choice("ACGT") for _ in range(240))
rng = random.Random(12)
g2 = "".join(rng.choice("ACGT") for _ in range(240))
open("genes.fa", "w").write(
    f">labcur1 synthetic tetracycline efflux pump\n{g1}\n"
    f">labcur2 synthetic macrolide esterase\n{g2}\n"
)
```

```console
$ grep '>' genes.fa
>labcur1 synthetic tetracycline efflux pump
>labcur2 synthetic macrolide esterase
```

## 实战示例

### 1. plain FASTA 到第一条命中

把名为 `labgenes` 的数据库建进本地数据目录，然后筛查一条携带 `labcur1` 前 200 个
碱基加侧翼碱基的 query contig：

```console
$ gapit db build labgenes genes.fa --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/labgenes/sequences
gapit: self-check passed for labgenes
gapit: BLAST index built (nucl)
{"db":"labgenes","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labgenes"}
$ gapit screen query1.fa --db labgenes --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labgenes		synthetic tetracycline efflux pump	
```

命中行解码了 plain 表头：基因 `labcur1`，产物来自表头描述，`ACCESSION` 和
`RESISTANCE` 为空因为 plain 表头两者皆无。`%COVERAGE` 是 83.33，因为 query 只携带
基因 240 个碱基中的 200 个；coverage map 尾部的 `..` 正是那 40 个未覆盖碱基的示
意。

### 2. 用 --tsv 挂上元数据

元数据 TSV 为每个基因补充 accession 和功能类别。表头行必填，`gene` 列必填，
`accession` 和 `function` 按文件可选：

```console
$ cat meta1.tsv
gene	accession	function
labcur1	LAB-0001	tetracycline
labcur2	LAB-0002	macrolide
$ gapit db build labmeta genes.fa --datadir ./db --tsv meta1.tsv
gapit: generated /tmp/opencode/customdb-docs/db/labmeta/sequences
gapit: self-check passed for labmeta
gapit: BLAST index built (nucl)
{"db":"labmeta","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labmeta"}
$ gapit screen query1.fa --db labmeta --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labmeta	LAB-0001	synthetic tetracycline efflux pump	tetracycline
```

基因属于多个类别时 `function` 列用 `;` 分隔。两个类别会作为一个字符串落入
`RESISTANCE` 列：

```console
$ cat meta2.tsv
gene	accession	function
labcur1	LAB-0001	virulence;marker
labcur2	LAB-0002	macrolide
$ gapit db build labmulti genes.fa --datadir ./db --tsv meta2.tsv
gapit: generated /tmp/opencode/customdb-docs/db/labmulti/sequences
gapit: self-check passed for labmulti
gapit: BLAST index built (nucl)
{"db":"labmulti","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labmulti"}
$ gapit screen query1.fa --db labmulti --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labmulti	LAB-0001	synthetic tetracycline efflux pump	virulence;marker
```

### 3. abricate 风格的 ~~~ 表头

从 abricate 约定导出的文件原样构建。字段按 `~~~` 切分，最后一个字段之后的描述成为
产物：

```console
$ grep '>' legacy.fa
>oldlab~~~tetA_lab~~~SYN-100~~~TETRACYCLINE synthetic tetracycline pump
>oldlab~~~ermX~~~SYN-101~~~MACROLIDE synthetic methylase
$ gapit db build oldlab legacy.fa --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/oldlab/sequences
gapit: self-check passed for oldlab
gapit: BLAST index built (nucl)
{"db":"oldlab","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/oldlab"}
$ gapit screen query_tet.fa --db oldlab --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query_tet.fa	contigB	1	240	+	tetA_lab	1-240/240	===============	0/0	100.00	100.00	oldlab	SYN-100	synthetic tetracycline pump	TETRACYCLINE
```

第一个 `~~~` 字段是原始数据库名，会被丢弃：`db` 始终是你在构建的 NAME。基因、
accession、耐药类别全部流入筛查输出。

### 4. 经 gapit| 表头往返

原生记录是可移植的。从示例 1 写出的 `records.jsonl` 再生一份 FASTA（每条记录一个
`gapit|` 表头，序列和产物原样），然后把它重建为新的数据库：

```console
$ grep '>' ported.fa
>gapit|db=labgenes|gene=labcur1|acc=|func= synthetic tetracycline efflux pump
>gapit|db=labgenes|gene=labcur2|acc=|func= synthetic macrolide esterase
$ gapit db build labgenes_rt ported.fa --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/labgenes_rt/sequences
gapit: self-check passed for labgenes_rt
gapit: BLAST index built (nucl)
{"db":"labgenes_rt","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labgenes_rt"}
$ gapit screen query1.fa --db labgenes_rt --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labgenes_rt		synthetic tetracycline efflux pump	
```

筛查输出与示例 1 一致（表头中的 `db` 值被重定向到新的 NAME）。空的 `acc=` 和
`func=` 段没问题，携带 accession 和 function 值的 `gapit|` 表头同样原样往返。

### 5. 蛋白质数据库（blastx）

氨基酸输入会被自动检测。注意回执写的是 `"dbtype":"prot"`，BLAST 索引行写的是
`(prot)`：

```console
$ grep '>' toxins.faa
>toxA synthetic pore-forming toxin
$ gapit db build toxins toxins.faa --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/toxins/sequences
gapit: self-check passed for toxins
gapit: BLAST index built (prot)
{"db":"toxins","records":1,"dbtype":"prot","destination":"/tmp/opencode/customdb-docs/db/toxins"}
$ gapit db build toxins2 toxins.faa --datadir ./db --dbtype prot
gapit: generated /tmp/opencode/customdb-docs/db/toxins2/sequences
gapit: self-check passed for toxins2
gapit: BLAST index built (prot)
{"db":"toxins2","records":1,"dbtype":"prot","destination":"/tmp/opencode/customdb-docs/db/toxins2"}
```

筛查蛋白质数据库运行 `blastx`，所以 query 文件必须是核苷酸。这里 query contig 携
带该毒素的编码序列，按每个氨基酸一个密码子反向翻译，带侧翼碱基：

```console
$ gapit screen query_toxin.fa --db toxins --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query_toxin.fa	contigT	31	390	+	toxA	1-120/120	===============	0/0	100.00	100.00	toxins		synthetic pore-forming toxin
```

读这一行：`START`/`END` 是 contig 上的核苷酸坐标（30 bp 侧翼，随后是 360 bp 的
CDS），而 `COVERAGE` 和 `%COVERAGE` 按 120 个残基的蛋白统计氨基酸。`--dbtype
nucl|prot` 为的是字母表本身会误导启发式的边缘情况；常规输入下自动检测与显式参数的
结论一致。

### 6. 压缩输入

`.gz`（和 `.bz2`）输入与 plain 文件的构建方式完全一致：

```console
$ gapit db build labgz genes.fa.gz --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/labgz/sequences
gapit: self-check passed for labgz
gapit: BLAST index built (nucl)
{"db":"labgz","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labgz"}
```

### 7. 迭代：用 --force 重建

往 FASTA 里加第三个基因，重建示例 1 的 `labgenes`。已有的数据库绝不被静默覆盖：

```console
$ grep '>' genes_v2.fa
>labcur1 synthetic tetracycline efflux pump
>labcur2 synthetic macrolide esterase
>labcur3 synthetic fosfomycin thiolase
$ gapit db build labgenes genes_v2.fa --datadir ./db
{"schema":"gapit.error/1","code":"DB_ALREADY_EXISTS","message":"won't overwrite existing database labgenes (use --force)","context":{"db":"labgenes"}}
$ echo $?
4
```

信封就是标准的 `gapit.error/1` 形状（退出码 4，数据库错误）。

`--force` 删除后原地重建。示例 2 的元数据 TSV 依然生效：

```console
$ gapit db build labgenes genes_v2.fa --datadir ./db --tsv meta1.tsv --force
gapit: generated /tmp/opencode/customdb-docs/db/labgenes/sequences
gapit: self-check passed for labgenes
gapit: BLAST index built (nucl)
{"db":"labgenes","records":3,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labgenes"}
$ gapit screen query1.fa --db labgenes --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labgenes	LAB-0001	synthetic tetracycline efflux pump	tetracycline
```

回执现在报告 `records: 3`。编辑 `records.jsonl` 或 FASTA，用 `--force` 重建，每个
下游产物都会再生。

### 8. 元数据的边缘规则

同一基因的重复行：**第一行胜出**，stderr 警告（`--quiet` 静默），构建成功。筛查输
出证明哪一行存活：

```console
$ cat dup.tsv
gene	accession	function
labcur1	KEEP-1	tetracycline
labcur1	LOST-2	macrolide
$ gapit db build labdup genes.fa --datadir ./db --tsv dup.tsv
WARNING: duplicate gene 'labcur1' in metadata TSV: keeping the first row
gapit: generated /tmp/opencode/customdb-docs/db/labdup/sequences
gapit: self-check passed for labdup
gapit: BLAST index built (nucl)
{"db":"labdup","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labdup"}
$ gapit screen query1.fa --db labdup --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labdup	KEEP-1	synthetic tetracycline efflux pump	tetracycline
```

TSV 里出现 FASTA 不携带的基因会警告并跳过。FASTA 是存在与否的唯一事实：

```console
$ cat ghost.tsv
gene	accession	function
labcur1	LAB-0001	tetracycline
ghostgene	G-999	virulence
$ gapit db build labghost genes.fa --datadir ./db --tsv ghost.tsv
WARNING: gene 'ghostgene' in metadata TSV not found in FASTA: skipped
gapit: generated /tmp/opencode/customdb-docs/db/labghost/sequences
gapit: self-check passed for labghost
gapit: BLAST index built (nucl)
{"db":"labghost","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labghost"}
$ grep -c '"gene"' db/labghost/records.jsonl
2
```

### 9. 不带 --datadir 时构建落在哪

数据目录按调用解析：`--datadir`，然后 `$GAPIT_DATADIR`，然后
`~/.local/share/gapit/db`。本例把 `GAPIT_DATADIR` 指向一个临时目录（不设置时，同一
命令会写进默认的 `~/.local/share/gapit/db`）：

```console
$ export GAPIT_DATADIR=/tmp/opencode/customdb-default/db
$ gapit db build labenv genes.fa
gapit: generated /tmp/opencode/customdb-default/db/labenv/sequences
gapit: self-check passed for labenv
gapit: BLAST index built (nucl)
{"db":"labenv","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-default/db/labenv"}
$ gapit db outdated
NAME	FETCHED_AT	AGE_DAYS	STATUS
labenv	2026-09-27T14:46:27Z	0.00	ok
```

与筛查不同，构建会创建缺失的数据目录而不是失败，所以新机器第一次构建即可自举。这
个临时数据目录在截图后已删除。

### 10. 收尾流水线：两个样本，一个矩阵

用示例 2 的 `labmeta` 筛查两个 contig 文件，保存报告表，再用
[gapit summary](./summary.md) 折叠成汇总矩阵：

```console
$ gapit screen sample1.fa --db labmeta --datadir ./db > sample1.tsv
$ gapit screen sample2.fa --db labmeta --datadir ./db > sample2.tsv
$ cat sample1.tsv
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
sample1.fa	SAM001	41	280	+	labcur1	1-240/240	===============	0/0	100.00	100.00	labmeta	LAB-0001	synthetic tetracycline efflux pump	tetracycline
$ cat sample2.tsv
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
sample2.fa	SAM002	16	255	+	labcur2	1-240/240	===============	0/0	100.00	100.00	labmeta	LAB-0002	synthetic macrolide esterase	macrolide
$ gapit summary sample1.tsv sample2.tsv
#FILE	NUM_FOUND	labcur1	labcur2
sample1.tsv	1	100.00	.
sample2.tsv	1	.	100.00
```

## 规则与边缘行为

### 表头自动检测

表头类型逐记录检测，所以一份 FASTA 可以混用全部三种格式：

| FASTA 表头 | Gene | Accession | Function | Product |
|---|---|---|---|---|
| `>mixplain plain header gene` | `mixplain` | 无 | 无 | `plain header gene` |
| `>oldlab~~~fromtilde~~~SYN-100~~~TETRACYCLINE tilde gene` | `fromtilde` | `SYN-100` | `TETRACYCLINE` | `tilde gene` |
| `>gapit\|db=elsewhere\|gene=fromtag\|acc=SYN-200\|func=ampicillin;gentamicin tagged gene` | `fromtag` | `SYN-200` | `ampicillin`、`gentamicin` | `tagged gene` |
| `>secondplain`（裸） | `secondplain` | 无 | 无 | `--description`，否则 `secondplain` |

plain 表头有描述文本时，id 之后的文本就是产物。没有时填 `--description TEXT`，两
者皆无时基因名成为产物：

```console
$ grep '>' bare.fa
>baregene
$ gapit db build bare_desc bare.fa --datadir ./db --description "lab-curated reference gene"
gapit: generated /tmp/opencode/customdb-docs/db/bare_desc/sequences
gapit: self-check passed for bare_desc
gapit: BLAST index built (nucl)
{"db":"bare_desc","records":1,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/bare_desc"}
$ grep -o '"product":"[^"]*"' db/bare_desc/records.jsonl
"product":"lab-curated reference gene"
$ gapit db build bare_plain bare.fa --datadir ./db --force
gapit: generated /tmp/opencode/customdb-docs/db/bare_plain/sequences
gapit: self-check passed for bare_plain
gapit: BLAST index built (nucl)
{"db":"bare_plain","records":1,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/bare_plain"}
$ grep -o '"product":"[^"]*"' db/bare_plain/records.jsonl
"product":"baregene"
```

以 `gapit|` 开头但不是合法标记语法的表头使构建失败，退出码 4：

```console
$ gapit db build labbroken broken_tag.fa --datadir ./db
{"schema":"gapit.error/1","code":"HEADER_MALFORMED","message":"malformed gapit/v1 sequence header: segment_without_key","context":{"seqid":"gapit|nonsense","reason":"segment_without_key"}}
```

上表那份混合文件按各自的规则逐记录构建，生成的 `sequences` 投影展示了重定向
（`db` 始终是构建的 NAME，原始的 `elsewhere` 标记被丢弃）：

```console
$ gapit db build labmix mixed.fa --datadir ./db 2>/dev/null
{"db":"labmix","records":4,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labmix"}
$ grep '>' db/labmix/sequences
>gapit|db=labmix|gene=mixplain|acc=|func= plain header gene
>gapit|db=labmix|gene=fromtilde|acc=SYN-100|func=TETRACYCLINE tilde gene
>gapit|db=labmix|gene=fromtag|acc=SYN-200|func=ampicillin%3Bgentamicin tagged gene
>gapit|db=labmix|gene=secondplain|acc=|func= secondplain
```

记录保持输入顺序，所以 `records.jsonl` 的第 N 行就是第 N 条 FASTA 记录。允许重复
基因名；携带两次 `twingene` 的 FASTA 构建出两记录的数据库：

```console
$ gapit db build labtwin dupname.fa --datadir ./db --force
gapit: generated /tmp/opencode/customdb-docs/db/labtwin/sequences
gapit: self-check passed for labtwin
gapit: BLAST index built (nucl)
{"db":"labtwin","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labtwin"}
$ grep -c '"gene":"twingene"' db/labtwin/records.jsonl
2
```

### dbtype 解析

给了 `--dbtype` 就以它为准；否则 abricate 的分子类型启发式从输入序列本身判断（默
认 `nucl`，除非字母表明是蛋白质）。蛋白质数据库得到 `.pin` BLAST 索引，筛查经
`blastx` 以核苷酸 query 进行（示例 5）。

### 数据库目录里会落下什么

```console
$ ls db/labmeta
gapit-manifest.json
records.jsonl
sequences
sequences.ndb
sequences.nhr
sequences.nin
sequences.njs
sequences.not
sequences.nsq
sequences.ntf
sequences.nto
```

`records.jsonl` 是可编辑的事实来源，`sequences` 是它的 `gapit/v1` 表头投影；
`sequences.n*` 文件是 BLAST 索引。两个文件契约的字段级细节见[数据库](./databases.md)。
序列原样存储，不做提供商风格的归一化：这是你策展的事实，manifest 以
`source_urls: ["local"]` 为构建背书。

## 故障排查

| 症状 | 错误（退出码） | 修复 |
|---|---|---|
| 构建拒绝已存在的数据库 | `DB_ALREADY_EXISTS`（4） | 用 `--force` 重建（示例 7）。 |
| 元数据 TSV 被拒 | `METADATA_MALFORMED`（5） | 首行必须是含 `gene` 列的表头。 |
| 空输入导致构建失败 | `BUILD_INVALID`（4） | FASTA 解析出零条记录；检查文件是否有 `>` 表头。 |
| FASTA 或 TSV 路径缺失 | `INPUT_NOT_FOUND`（5） | 信封的 `context.file` 指出不可读的路径。 |
| `gapit\|` 表头导致构建失败 | `HEADER_MALFORMED`（4） | 修正标记表头；plain 和 `~~~` 表头不能用 `gapit\|` 语法。 |
| 用自定义库筛查一无所获 | 无错误，退出码 0 | 空结果表是正常输出。放宽 `--minid`/`--mincov`，或检查 query 与链方向是否匹配；见[筛查](./screen.md)。 |

错误信封记录在[输出](./outputs.md)；完整的筛查参数集见[筛查](./screen.md)。
