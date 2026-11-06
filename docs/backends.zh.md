# 多后端分型架构

本指南解释 gapit 如何在同一个分型体系中组合比对与 k-mer 后端，以及如何设计
利用正确后端解决各自判别问题的自定义数据库。构建机制见
[自定义数据库](./custom-db.zh.md)；本页是架构配套文档。

## 总览

```text
                        gapit screen
                            │
              ┌─────────────┼─────────────┐
              │             │             │
         FASTA 输入    FASTQ 输入    FASTQ 输入
              │             │             │
         ┌────┴────┐   ┌────┴────┐   ┌────┴─────┐
         │ 基因库  │   │ 基因库  │   │ cluster库│
         │ BLASTn  │   │minimap2 │   │+kmer.json│
         │         │   │  (sr)   │   │ k-mer    │
         └────┬────┘   └────┬────┘   │containment│
              │             │        └────┬─────┘
              │             │             │
              └─────────────┼─────────────┘
                            │
                     gapit typing
                            │
                 ┌──────────┼──────────┐
                 │          │          │
            exact_set  weighted_  cluster_match
                       genes      (覆盖+同一性
                                   +key_genes
                                   +key_kmers)
```

两种库型、两种输入、三条算法路径——每种组合都可用：

| 库型 | FASTA（组装） | FASTQ（reads） |
|---|---|---|
| **gene** | BLASTn 比对 | minimap2 `-x sr` 比对 |
| **cluster** | minimap2 `-x asm20` 比对 + k-mer 命中率 | k-mer containment（无比对） |

**判别力**因路径而异：

| 判别问题 | 比对 | k-mer | 优胜者 |
|---|---|---|---|
| 基因存在/缺失 | ✅ | ✅ | 比对（基因调用、覆盖度） |
| 位点家族归属 | ✅ 覆盖/同一性 | ✅ containment | 均可 |
| **SNP 级对**（6A vs 6B） | ❌（~0.001 同一性差） | ✅（1.0 vs 0.0） | **key_kmers** |
| 多拷贝基因解析 | BLASTn `--count-same-copy` | N/A | 比对 |
| 断裂位点组装 | `merge_fragments` | 天然支持（k-mer 可加性） | k-mer |

## 库型

### 基因库（gene）

参考序列平面面板（扩增子、基因、蛋白）。筛查报告每基因的存在/
缺失及覆盖/同一性。分型规则是布尔型：

- `exact_set`：所有 `requires` 在场 + 所有 `excludes` 缺失 → 1.0
- `weighted_genes`：每基因带权重；加和为得分

**适用**：PCR 式标记面板（Penner 多重 PCR、Doumith 血清群）、毒力
因子面板、小基因簇。

**内参基因**：基因库支持 reads 模式相对深度门控——一个单拷贝对照
基因（`lpxA`、`tuf`）。见 [Reads](./reads.zh.md)。

### 位点簇库（cluster）

完整基因组位点（10–50 kb）带 CDS 注释。筛查通过 minimap2 比对报告
每位点覆盖/同一性/基因调用；分型规则在这些之上打分：

- `cluster_match`：每位点覆盖 + 同一性 + `key_genes` + `key_kmers`
  分量加权和；`key_genes.groups` 可列出互斥的备选基因集（如近乎相同
  的 O1 Inaba/Ogawa LPS 位点），取命中的最佳组计分

**适用**：多糖生物合成位点（荚膜、LPS、O 抗原）、毒力岛、任何
基因组成重要的大型可变区。

**配置**：`cluster.json` 设定覆盖模式（`primary_secondary` 用于多
同源位点）和质量门。

## 算法后端详解

### 1. BLASTn（gene，FASTA）

主力引擎：每基因最佳 HSP、覆盖/同一性阈值、可选 `--merge-fragments`
处理断裂组装。可靠、广泛理解、abricate 兼容基线。

### 2. minimap2 `-x sr`（gene，FASTQ）

短读比对用于基因面板。每基因深度/广度；`screen.json` reads 段配置
MAPQ、深度和 `min_relative_depth`（以对照基因为深度锚）过滤假阳性。

### 3. minimap2 `-x asm20 --cs`（cluster，FASTA）

组装到位点比对：跨 contig 联合覆盖、每基因判定（present/partial/
absent）、插入/缺失图谱。`primary_secondary` 覆盖模式包含次级比对
（对多同源位点如 Penner CPS 复合群至关重要）。

### 4. k-mer containment（cluster，FASTQ）

完全不做比对。查询的规范 k-mer（2-bit 打包、滚动计数）与每位点
的全 k-mer 集求交：

- **覆盖度** = 位点 k-mer 在查询中存在的比例
- **同一性** = `containment^(1/k)`（mash 式 ANI 估计）
- **key_kmers** = 位点*差异*集在查询中的命中率

此路径彻底消除竞争映射问题（k-mer 存在性是集合成员判定，不是最佳
比对）且直接使用原始 reads。

### 5. k-mer 差异集（cluster，FASTA 和 FASTQ 通用）

`key_genes` 的 k-mer 粒度类比物：一个位点对其家族兄弟独有
的 k-mer。

```json
{"key_kmers": {"weight": 0.25, "k": 31, "set": "06A.diff"}}
```

对 SNP 级对（6A 与 6B 差异为一个核苷酸），差异集包含跨越该 SNP
的 ~31 个 k-mer：

- 6A 查询：`key_kmers(06A) = 1.0`，`key_kmers(06B) = 0.0`
- 6B 查询：`key_kmers(06B) = 1.0`，`key_kmers(06A) = 0.0`

二值分离——任何歧义边际都吞不掉这个差距。

## 构建自定义数据库：决策表

| 你的判别问题 | 库型 | 分型模型 | 关键分量 |
|---|---|---|---|
| PCR 式标记（存在/缺失） | gene | `exact_set` | — |
| 梯度基因贡献 | gene | `weighted_genes` | — |
| 大位点，基因组成变化 | cluster | `cluster_match` | `key_genes` |
| 大位点，SNP 级对 | cluster + kmer.json | `cluster_match` | `key_kmers` |
| 大位点，基因 + SNP | cluster + kmer.json | `cluster_match` | `key_genes` + `key_kmers` |
| 仅 reads 工作流（免组装） | cluster + kmer.json | `cluster_match` | `key_kmers`（containment 路径） |

### 分步：带 k-mer 底座的 cluster 库

```bash
# 1. 准备位点 FASTA（或带 CDS 注释的 GBK/GFF3）
# 2. 设计分型规则（gapit.typing/3）
# 3. 声明 k-mer 集（gapit.kmerdb/1）
# 4. 构建

gapit db build my_sero loci.gff3 \
  --typing typing.json \
  --kmerdb kmer.json \
  --description "My organism serotyping"

# 5. 筛查（FASTA）
gapit screen -d my_sero -C 60 -G 90 assembly.fna

# 6. 筛查（FASTQ —— 免组装）
gapit screen -d my_sero --r1 reads_R1.fastq --r2 reads_R2.fastq
```

`kmer.json` 声明一个 k 大小和命名集合：

```json
{
  "schema": "gapit.kmerdb/1",
  "k": 31,
  "sets": {
    "typeX.diff": {"of": "typeX", "against": ["typeY", "typeZ"]}
  }
}
```

构建时 gapit 计算每位点的全 k-mer 集和每个声明的差异集，以排序
uint64 数组存储于 `kmers/`。分型规则引用集合名；引擎在构建时校验
（TYPING_UNKNOWN_KMERSET）。

### 多路径规则的权重设计

必须同时适用于 FASTA 和 FASTQ 的规则应考虑 reads 模式不产生基因
调用（`key_genes` 得 0 分）：

```json
{
  "coverage":  {"weight": 0.45, "floor": 60},
  "identity":  {"weight": 0.20, "floor": 90},
  "key_genes": {"weight": 0.10, "genes": ["geneA", "geneB"]},
  "key_kmers": {"weight": 0.25, "k": 31, "set": "typeX.diff"}
}
```

- FASTA 路径最大值：0.45 + 0.20 + 0.10 + 0.25 = **1.00**
- FASTQ 路径最大值：0.45 + 0.20 + 0 + 0.25 = **0.90**
- 设 `cutoff` ≤ 0.75 使两条路径均可达到可信判定

## 实战示例

| 库名 | 物种 | 库型 | 规则 | k-mer | SNP级 | 输入 |
|---|---|---|---|---|---|---|
| `cje_cps` | *C. jejuni* | gene | 38 `exact_set` | — | — | FASTA+FASTQ |
| `cje_cpsv2` | *C. jejuni* | gene | 31 `exact_set` | — | — | FASTA+FASTQ |
| `cje_cpsv2c` | *C. jejuni* | cluster | 31 `cluster_match` | — | — | FASTA |
| `lpn_sero` | *L. pneumophila* | cluster | 13 `cluster_match` | — | — | FASTA |
| `lpn_vf` | *L. pneumophila* | gene(蛋白) | —（仅筛查） | — | — | FASTA |
| `str_lancefield` | *链球菌* | gene | 5 `exact_set` + tuf | — | — | FASTA+FASTQ |
| `spn_sero` | *S. pneumoniae* | cluster | 100 `cluster_match` | ✅ k=31 | ✅ | **FASTA+FASTQ** |

### spn_sero：多后端完整展示

100 条 SeroBA cps 位点参考、100 个单型规则带差异 k-mer 集。同一份
`gapit.typing/3` 文档在组装（比对路径：覆盖+同一性+key_genes+
key_kmers）和 reads（k-mer containment 路径：估计覆盖/同一性+
key_kmers）上评估完全一致。SNP 级对（6A/6B、9A/9V、15B/15C、
18B/18C…）两条路径均可拆分。

## 来源规则（强制）

每个自定义数据库必须在 manifest 中写明来源元数据，在 `info.json`
中写简短描述。完整规则见
[自定义数据库 → 来源与许可](./custom-db.zh.md#来源与许可元数据强制)。

## 延伸阅读

- [自定义数据库](./custom-db.zh.md) — 构建机制、来源规则
- [分型](./typing.zh.md) — 规则模型详解
- [Reads](./reads.zh.md) — 基因库 FASTQ 筛查
- [数据库](./databases.zh.md) — 内置/提供商数据库、命名
