# gapit v0.5.6 代码质量深度分析报告

> **分析对象**：`/home/mark/repos/github/gapit`，commit `1924743`（tag `v0.5.6`，2026-10-20）  
> **分析日期**：2026-10-22  
> **分析环境**：pixi 0.65.0；Python 3.14.7（默认）/ 3.11.16；ruff 0.16.8；basedpyright 1.40.1（pyright 1.1.414）；pytest 9.1.1；rich 15.0.0  
> **分析方式**：96 个模块逐层阅读；运行仓库自身门禁；AST 度量（复杂度/注解/导入图）；定向复现；解析器模糊测试；CLI 冒烟与契约探测；git/文档/发布流程审查。  
> **说明**：本报告是工程与契约视角的静态 + 动态分析，不是生物信息学正确性认证；parity 环境（真实 abricate）本次未安装，未执行字节级差分。

> **修复进展（2026-10-22，后续提交）**：本报告列出的 P1-1/P1-2/P1-3，以及 P2-1～P2-6 已完成代码与回归测试修复；P3-9 的文档计数、错误码大小写、`stderr.txt`、`fmt -- --check` 也一并处理。修复后全量测试为 **1209 passed, 9 skipped**（1218 collected），`lint`/`typecheck`/`fmt-check` 全绿。P2-7（`cluster/2` 契约版本）与 P3 的架构重构/未版本化 JSON/MCP 对齐仍保留为后续工作；P2-6 的 checksum pinning 需要上游哈希，当前落地为 HTTP 警告 + 1 GiB 下载上限。下文保留修复前证据，作为问题记录。

---

## 1. 执行摘要

**总体判断：这是一个工程质量明显高于同类 0.5.x 开源项目的代码库。** 严格的类型检查、1183 个测试、golden/parity 资产、双语文档、许可证审计与确定性输出设计，都说明作者对“公共契约”有清醒认识。

但当前版本存在 **3 个必须优先修复的 P1 问题**：

1. **强制重建失败会留下“旧 manifest 认证新数据”的静默不一致**，破坏数据库 provenance 与后续筛查结果的可信度。
2. **CSV 输出不是合法 CSV**：使用裸 `join`，没有引号/转义；成对 reads 的 sample 字段甚至必然多出一列。
3. **提交锁定的依赖下测试门禁是红的**：默认环境全量 4 failed；py311 环境仅跑 `test_progress.py` 也同样 4 failed。CI 最近一次成功很可能是缓存环境掩盖。

此外还有若干 P2：click 参数解析错误绕过 JSON 错误信封、畸形 GBK/GFF 触发未包装异常、控制字符导致 TSV/Markdown 列错位、`make_blast_db` 先删索引后检查依赖、reads 模式 PAF 全量物化与文档“流式”声明不符、三个 provider 走明文 HTTP 且无完整性校验、`gapit.cluster/2` 在 patch 版本中被破坏性重塑。

### 1.1 评分卡（主观加权，10 分制）

| 维度 | 分数 | 简评 |
|---|---:|---|
| 架构与模块化 | 8.0 | 分层清晰、无运行时循环；但命令层臃肿、长参数列表、批处理逻辑重复 |
| 类型安全 | 9.5 | strict 0 error；443 个函数仅 1 个缺返回注解；`Any` 被限制在 JSON/MCP 边界 |
| 错误处理与健壮性 | 7.0 | 语义错误信封优秀；但 click 解析错误、解析器崩溃、写路径 OSError 未统一 |
| 测试体系 | 7.5 | 数量与结构优秀；但当前红、无覆盖率门禁、缺少对抗性/属性测试 |
| 文档与契约 | 7.5 | SPEC/AGENTS/docs 质量高；但版本、数量、错误码大小写、包内资产存在漂移 |
| 安全与供应链 | 8.0 | 无 shell、argv 调用、路径防护、SHA256 安装；但 HTTP 源无校验、无下载上限 |
| 性能与可扩展性 | 7.0 | 小数据很快；PAF/BLAST 全量物化是规模瓶颈；启动 0.42s |
| 工程流程/发布 | 8.0 | CI 矩阵、parity、docs strict、CHANGELOG 都好；release 无测试门禁、版本多处维护 |
| **综合** | **7.9（B+）** | **适合研究/内部试用；不建议当前状态用于无人值守关键流水线** |

### 1.2 最值得肯定的六件事

1. **类型纪律罕见地好**：`basedpyright --strict` 0 error；没有 `# type: ignore`、没有 `cast` 消音；`Any` 只出现在 MCP JSON-RPC 与 CARD 异构 JSON 的受控边界，并附有理由注释。
2. **测试不是摆设**：1183 个测试项、约 3356 个 assert、61 个 golden 文件、独立 parity harness；测试 docstring 普遍采用 Given/When/Then，能看出 TDD 痕迹。
3. **“parse, don't validate”落实到位**：BLAST/PAF/FASTA/GBK/GFF/typing/manifest 都在边界变成 pydantic 模型；89 个模型不是装饰，而是真实数据契约。
4. **确定性设计**：稳定排序、LF、UTF-8、schema 先行；JSON 的 `created_at` 只放在元数据块；`--jobs` 下仍按输入顺序输出，并有字节一致性测试。
5. **关键写入有原子性意识**：`sequences` 用临时文件 + `os.replace`，`db install` 先校验 SHA256 再原子替换；构建后有 header/序列自检。
6. **许可证与包体纪律**：19 个 provider 全部 fetch-on-demand；7 个 bundle 有 30 MiB 预算测试；CARD/VFDB/Kaptive 内容不随 MIT wheel 分发。这在同类工具里很少见。

---

## 2. 项目画像

| 指标 | 数值 |
|---|---:|
| 源码 LOC（`src/gapit`，含空行/注释） | 13,227 |
| 测试 LOC（`tests`） | 23,713 |
| 测试/源码比 | 1.79 : 1 |
| Python 模块（`src/gapit`，含包初始化） | 96 |
| 函数/方法 | 443 |
| pydantic 模型 | 89 |
| 测试项（`pytest --collect-only`） | 1,183 |
| 测试文件 | 99 |
| golden / parity / test-data 文件 | 61 / 8 / 111 |
| git 提交 | 78，单作者，2026-09-15 → 2026-10-20 |
| 版本 | 0.5.6（Beta） |
| 包内数据库内容 | 20.6 MiB（7 个 bundle，预算 30 MiB） |
| CLI 启动耗时（`gapit --version`） | ≈0.42 s |

**架构分层（实际依赖方向）**：

```
CLI (cmd_*/cli.py)
  → use-cases (screening*, typing_ops, db_*_ops, bundled)
    → engines / domain (blast, hits, cluster, typing_*)
      → boundary models + external tools (db, fasta, seqconvert, paf, dbcodec)
        → leaf utilities (errors, dispatch, proctools, engines)
```

整体是“命令层薄、use-case 厚、领域纯函数”的结构；`errors.py`、`dispatch.py`、`engines.py` 作为叶子模块避免了大部分循环导入。没有发现真正的运行时循环依赖（AST 图中看到的两个 SCC 分别由 `TYPE_CHECKING` 和函数内延迟导入造成）。

---

## 3. 质量门禁与验证结果

### 3.1 仓库自带门禁

| 检查 | 命令 | 结果 |
|---|---|---|
| Lint | `pixi run lint` | ✅ `All checks passed!` |
| 格式 | `pixi run fmt --check`（直接传给 ruff） | ✅ 253 files already formatted |
| 格式（CI/文档写法） | `pixi run fmt -- --check` | ❌ 本机 exit 2：`Failed to format --check`（pixi 把 `--` 也传了进去） |
| 类型 | `pixi run typecheck` | ✅ `0 errors, 0 warnings, 0 notes` |
| 测试 | `pixi run test` | ❌ **4 failed, 1170 passed, 9 skipped in 219.97s** |
| 测试收集 | `pytest --collect-only -q` | 1183 tests collected |

### 3.2 4 个失败测试

```
tests/test_progress.py::test_tty_stderr_activates_the_bar
tests/test_progress.py::test_screen_bar_keeps_stdout_byte_identical_and_suppresses_notes
tests/test_progress.py::test_reads_bar_keeps_streaming_header_first_and_suppresses_notes
tests/test_progress.py::test_fetch_bar_advances_per_db_and_keeps_receipts
```

在 `py311` 环境单独跑 `tests/test_progress.py`，同样 **4 failed, 5 passed**。

**直接原因**（不是 flaky）：

- `_LiveBatch.describe()` 直接 `progress.update(task_id, description=text)`，会把初始的 `"Screening"` / `"Fetching"` 覆盖成当前文件/样本/数据库名。
- 测试在 `describe("b.fa")` 之后仍断言 `"Screening" in export_text(...)`；锁定环境 rich 15.0.0 的 `export_text(clear=False)` 只保留最终帧，于是只剩 `b.fa`。
- 这不是偶发时序：`test_bar_renders_current_description_and_counts` 断言 `"full.fa"` 是通过的，说明同一套渲染在“断言当前 description”时正常，失败的是“断言初始 description 仍在”的过时假设。

**影响**：项目文档反复强调“每次改动必须 lint/typecheck/test 全绿”，但当前提交在锁定依赖下不满足；CI 最近一次 `gates` 成功与本机结果矛盾，最可能是 CI 复用了旧缓存环境（日志显示 py311 cache restore 失败，default 可能命中）。需要在**无缓存干净环境**重跑确认。

### 3.3 补充扫描（非仓库门禁）

| 扫描 | 结果 | 解读 |
|---|---|---|
| `ruff --select S` | 15 项：9×S101 assert、3×S603 subprocess、2×S310 urlopen、1×S607 部分路径 | S603/S607 是“无 shell”设计的误报；S101 中有 4 处位于 GBK 解析器，见 P2 |
| `ruff --select TRY,BLE` | 114×TRY003、1×BLE001 | TRY003 是项目自定义异常类的正常用法；BLE001 仅 `dispatch` 顶层，合理 |
| `ruff --select C901`（阈值 10） | 9 个函数超阈值 | 与 AST 度量一致 |
| 解析器 fuzz | 7 类解析器 × 300 个变异输入，317 个未包装异常 | 140 BLAST、157 PAF、20 GBK；GFF 首轮 0，但手工 start=0 触发 ValidationError |
| PAF 内存 | 100k 行（含 500bp `cs`）= **150.9 MB** | 约 1.5 KB/行，见 P2 性能问题 |
| CLI 启动 | `gapit --version` ≈0.42 s，`import gapit.cli` ≈377 ms | 所有 provider/命令被急切导入 |

---

## 4. 主要优点（展开）

### 4.1 类型与边界建模

- `basedpyright` strict 覆盖 `src` 与 `tests`，0 error。
- 443 个函数中只有 `typing_ops._resolve_document` 缺返回注解；没有缺失的参数注解。
- `Any` 只出现在 `mcp*.py` 的 JSON-RPC 容器、`providers/card.py` 的异构 JSON、`dispatch.py` 的 click 装饰器 TypeVar，并都有注释说明。
- `BlastRow`、`PafRecord`、`Record`、`Manifest`、`TypingDocument`、`ClusterReport`、各类 `*Document` 在边界完成解析/校验，核心逻辑不碰裸字符串。
- pydantic 的 `frozen=True` 贯穿模型，减少意外共享可变状态。

### 4.2 测试结构

- 测试按“纯函数 / 格式化 / CLI / provider / 集成 / golden / parity”分层；`tests/data` 只有 1.4 MB，离线可跑。
- 61 个 golden 文件覆盖 TSV/CSV/JSON/MD；8 个 parity 脚本保留 abricate 基线回归。
- 测试普遍先写 docstring 描述 Given/When/Then，且大量断言是“字节完全相等”，不是松散包含。
- `--jobs` 顺序、stdout 纯净、流式 chunk 拼接字节一致等并发/流式契约都有专门测试。
- 9 个 skip 均为“本地大样本/真实基因组/差分二进制不存在”，不是掩盖失败。

### 4.3 确定性、原子性与自检

- 排序键明确（`summary` 按 row key，`report` 按 `(sequence, start)`，cluster 按 coverage/identity/bp/id）。
- `generate_sequences`/`_write_locus_fasta` 使用 `NamedTemporaryFile` + `os.replace`，失败不留下半截 `sequences`。
- `build_database` 在写 manifest 前做 header/sequence 自检；`db install` 做流式 SHA256 后原子替换。
- 输出层统一 LF、UTF-8，JSON `schema` 先行；`--nopath` 只影响显示列。
- `tests/test_streaming_chunks.py` 对 chunk 拼接与缓冲返回值做字节一致性断言，这是很多 CLI 项目缺失的。

### 4.4 安全基线

- 所有外部工具调用都是 argv list，无 `shell=True`、无 `os.system`、无 `eval/exec/pickle/yaml.load`。
- minimap2 输入路径显式 `absolute()`，避免 `-d` 之类文件名被解析为选项；`db build` 明确拒绝含 `/`、`\` 的名字，以及恰好为 `.`/`..` 的名字。
- `db install` 是“校验 SHA256 后原子安装”，`download_file` 用 `.part` + `os.replace`，失败清理临时文件。
- ZIP 解压使用标准库 `extractall`（现代 Python 会清理绝对路径/`..`），tar 只 `extractfile` 指定成员；未发现直接路径穿越 RCE。
- 许可证政策（PD/Apache/BSD/MIT 才可入 wheel；GPL/非商业内容 fetch-on-demand）在代码、测试、文档三处一致，属于加分项。

### 4.5 文档与工程资产

- `SPEC.md`（32 KB）把 abricate 行为、公式、去重规则、边界 quirk 写成了冻结契约；`AGENTS.md`（20 KB）给出 contributor/agent 规则；`PLAN.md`、`CHANGELOG.md`（51 KB）齐全。
- `docs/` 13 个页面 × 中英双语（26 个文件）+ `mkdocs --strict` 的 Docs workflow。
- CI 矩阵 3.11/3.13/3.14 + 独立 parity job；`ruff`/`basedpyright`/`pytest` 三条门禁都接入了 CI。
- 版本发布使用 PyPI trusted publishing（OIDC），不是长期 API token。

---

## 5. 问题清单（按严重度）

> **P1**：破坏公共契约/数据完整性，或使质量门禁失效；建议 1.0 前必须修。  
> **P2**：明确缺陷或特定场景风险；建议近期修。  
> **P3**：可维护性、文档、流程改进；可排期。

### P1-1 强制重建失败会留下“旧 manifest 认证新数据”的静默不一致

**证据**

- `src/gapit/dbbuild.py:206-239`：manifest 在 `sequences` → self-check → SHA → makeblastdb 之后才写，设计意图是“last certifies”。
- `src/gapit/db_build_ops.py:154-203`：`--force` 时**不删除/不失效旧 manifest**，直接覆盖 `records.jsonl`、`sequences`、`typing.json`/`floors.json`，再调用 `build_database`。
- `src/gapit/providers/common.py:175-209`、`src/gapit/clusterbuild.py:239-248` 同样没有旧证书失效逻辑。
- `src/gapit/db.py:109-124`：`make_blast_db` 在调用 `makeblastdb` **之前**就删除 `sequences.[np]??` 索引文件。

**复现（本机）**

先正常构建一个 DB，再用修改后的序列 `--force` 重建，但让 PATH 中没有 `minimap2`：

```text
$ PATH=<only makeblastdb/blastn/blastdbcmd> gapit db build forced db2.fa --datadir dd --force
exit=3
{"schema":"gapit.error/1","code":"MISSING_DEPENDENCY","message":"required binary not found on PATH: minimap2",...}

manifest present: yes
manifest sha: 528f9c9a...   # 旧序列的 sha
sequences sha: 60766174...  # 新序列的 sha
index files: 8              # 新索引已建好
```

如果 PATH 中连 `makeblastdb` 都没有，则索引被删除、旧 manifest 仍在：

```text
exit=3  required binary not found on PATH: makeblastdb
manifest after: present
sequences sha: 60766174...  # 新
manifest sha:  528f9c9a...  # 旧
index files: 0
```

**影响**

- `find_database` 只要求 `sequences` 可读，不校验 manifest 或 `sha256`；`db list` 对 registry 名称以 manifest 存在判定 installed，对本地目录只看 `sequences` 与 manifest 字段。两者都不会发现“manifest 旧、sequences 新”的不一致。
- 若失败发生在索引建好之后（如本复现的 minimap2 缺失），后续 `screen` 会**成功使用新序列/新索引**，但 `db list`/`db outdated`/manifest 仍报告旧 `sha256`、旧 `fetched_at`、旧 `n_records`——provenance 静默错误。
- 若失败发生在索引建好之前，目录“看起来已安装”但实际不可用；`screen` 会以 `DATABASE_NOT_INDEXED` 提示 `gapit setupdb`，但需要用户手动修复。
- 与项目自己“manifest written last certifies the artifacts”的核心不变量直接冲突。

**建议**

1. **构建到 staging 目录**（`<db>.build.<pid>/`），成功后先 rename 旧目录到 backup、再把 staging rename 到目标（或使用版本化目录 + symlink 切换）；失败清理 staging，旧目录完全不动。
2. 或者：开始 `--force` 前删除/重命名旧 manifest（tombstone），失败时恢复；并在 discovery 时校验 manifest 的 `sha256` 与 `sequences` 实际 hash（可缓存 stat）。
3. `make_blast_db` 先做依赖 preflight（`shutil.which`）再删除旧索引；索引也写到临时前缀再替换。
4. 增加“强制重建中途失败”的 fault-injection 测试：分别让 makeblastdb/minimap2/typing 校验失败，断言旧 DB 字节不变或 manifest 被正确失效。

---

### P1-2 CSV 输出不是合法 CSV（公共格式契约破损）

**证据**

所有 CSV 渲染器都是“把 TSV 分隔符替换成逗号”，没有 RFC 4180 的引号/转义：

- `src/gapit/formats/tsv.py:29-56`：`sep.join((... hit.gene, hit.accession, hit.function))`
- `src/gapit/formats/reads_tsv.py:52-55`：`sample = ", ".join(report.reads)`，随后 `sep.join`
- `src/gapit/formats/summary.py:55-59`：矩阵行直接 `sep.join`
- `src/gapit/formats/cluster.py:191-205`、`cluster2.py:203-208`：同样裸 join
- `tests/test_tsv.py:83-87` 的 docstring 明确假设 “no field ever contains a ',' (product cleanup strips them)”——**这个假设只对 product 成立**，对 gene/accession/function/file/sample 不成立。

**端到端复现（screen）**

用 `db build` 建一个 gene 含逗号的库：

```text
$ gapit screen contig.fa --db comma --datadir dd --format csv
#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE
/tmp/.../contig.fa,ctg1,1,160,+,gene,with,commas,1-160/160,===============,0/0,100.00,100.00,comma,A,1,genewithcommas,class,one
```

表头 15 列，数据行 17 个字段；任何 CSV 解析器都会把 `gene,with,commas` 拆成三列。

**端到端复现（summary）**

```text
$ gapit summary r.tsv --format csv
#FILE,NUM_FOUND,a,b
in.fa,1,+
```

表头 4 列，数据行 3 列；gene 名 `a,b` 把表头拆成了两列。

**配对 reads 的必然错误**

`reads_tsv_chunk` 把 sample 写成 `", ".join(report.reads)`。成对/多 lane 样本的第一个字段是 `s1, s2`，在 CSV 中必然被解析为两个字段（表头 9 列，数据行 10 个字段）。这不是“只有脏数据才触发”，而是普通 paired-end 就会触发。

**影响**

- 下游 `csv.reader`/pandas/Excel 会静默错列，且不会报错；对 agent 来说，列名与值错位比直接失败更危险。
- `db search --json` 等面是安全的（JSON 转义），但 CSV 是公开 `--format` 契约，不能以“TSV 兼容”为借口。
- 与 abricate 的 parity 只约束 TSV；CSV 是 gapit 自有扩展，可以独立修复，不需要破坏 parity。

**建议**

1. 用标准库 `csv.writer`（`QUOTE_MINIMAL`，`lineterminator="\n"`）统一渲染所有 CSV；输入侧（`summary`/`typing_input`）改用 `csv.reader` 以支持带引号字段。
2. 如果为了“可 grep 的简单格式”坚持裸 join，则必须在 `db build`/`screen` 边界拒绝 `,`（并明确文档），但这会损失基因名/样本名的表达能力，不如正确转义。
3. 增加对抗性测试：gene/accession/function/product/file/sample 分别含逗号、双引号、换行；paired-end sample 的 CSV 列数必须等于表头。
4. 同时修复 `tests/test_tsv.py:83-87` 的错误假设，否则测试会继续为 bug 背书。

---

### P1-3 锁定依赖下测试门禁为红，CI 绿色不可复现

**证据**

- 默认环境（Python 3.14.7，rich 15.0.0，pixi.lock 锁定）运行完整 `pytest`：**4 failed, 1170 passed, 9 skipped, 219.97s**。
- py311 环境单独运行 `tests/test_progress.py`：**4 failed, 5 passed**。
- 失败断言集中在：
  - `tests/test_progress.py:120` `assert "Screening" in text`
  - `tests/test_progress.py:177` `assert "Screening" in bar_text and "2/2" in bar_text`
  - `tests/test_progress.py:241` 同上
  - `tests/test_progress.py:271` `assert "Fetching" in bar_text and "2/2" in bar_text`
- 实际捕获的 bar 文本是 `b.fa ━━━ ... 2/2 100%`，没有 `Screening`/`Fetching`。
- `gh run view` 显示最近一次 release commit 的 `gates` 全部成功；与本机同 commit + 同 lock 的结果矛盾。

**根因**

`_LiveBatch.describe()` 覆盖 task description；“Screening”只存在于 `add_task` 的初始帧。rich 15 的 `export_text(clear=False)` 在录制 console 上只留下最终帧，因此断言必然失败。测试把“初始标题”当成“会一直保留的前缀”，与实现语义不一致。

**影响**

- “每次改动必须绿”的工程承诺当前不成立；任何基于 `pixi run test` 的发布/合并判断都不可靠。
- CI 若靠缓存环境“变绿”，会掩盖 lock 升级/新机器安装的问题；这是比 4 个测试本身更严重的流程风险。

**建议**

1. 修复测试或实现：要么 `describe` 写成 `"Screening " + text`（保留标题前缀），要么测试断言当前 description + 计数，不依赖初始帧。
2. CI 增加“无缓存/干净 env”周期任务（至少 release 前）；`setup-pixi` 缓存命中时也要校验 `pixi.lock` hash。
3. 给 `pixi run fmt -- --check` 增加独立 task（见 P3-8），避免 CI 命令与本地 pixi 版本耦合。

---

### P2-1 click 参数解析错误绕过 JSON 错误信封

**证据**

```text
$ gapit screen --bogus
Usage: gapit screen [OPTIONS] [FILES]...
Try 'gapit screen --help' for help.

Error: No such option '--bogus'. Did you mean '--jobs'?
exit=2

$ gapit screen tests/data/contigs/full.fa
Error: Missing option '--db' / '-d' / '-db'.
exit=2
```

而语义错误是符合契约的：

```text
$ gapit screen ... --minid 0
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"--minid must be in (0, 100]: got 0.0","context":{}}
```

**原因**

`dispatch()` 只包裹命令 body；click 在进入 body 之前就完成了 option/argument 解析并直接 `ctx.exit(2)`，不会经过 `dispatch`。

**影响**

AGENTS/README/docs 声称“Any failure prints one line of JSON to stderr”，但 agent 最容易遇到的失败——拼错 flag、漏传 `--db`——恰恰是纯文本。Agent 无法用统一的 `gapit.error/1` 解析，违背“agent-first”设计中心。

**建议**

- 自定义 `click.Group`/`Command`，在 `main()`/`invoke()` 外层捕获 `click.UsageError` 与 `click.ClickException`，渲染同一 envelope，保留 exit 2。
- 或在 entrypoint 用 `standalone_mode=False` 包一层，统一转换。
- 增加测试：未知 option、缺 required option、非法 choice、非法 int/float 都要断言 `gapit.error/1` 且 exit 2。

---

### P2-2 畸形 GBK/GFF 会抛出未包装异常，违背“typed boundary”

**证据（手工最小输入）**

```python
# LOCUS 行没有名字
"LOCUS \n//\n"                 -> IndexError: list index out of range
# 5 个空格组成的行
"LOCUS       TEST\n     \n//\n" -> IndexError: string index out of range
# 没有 LOCUS 就出现 ORIGIN
"...ORIGIN..."                 -> AssertionError
# GFF CDS start=0
start=0                        -> pydantic ValidationError (GeneFeature.start >= 1)
```

**代码位置**

- `src/gapit/gbfeatures.py:285-305`：`line.split()[1]`、`line[5]`、四处 `assert record is not None` 都直接作用于用户输入。
- `src/gapit/gffparse.py:68-70`：先构造 `GeneFeature`，再在 `_build_loci` 里做 `1 <= start <= end <= len` 检查；`start=0` 在构造时就被 pydantic 拒绝，成为 `ValidationError` 而不是 `INVALID_GFF`。
- fuzz 结果：2,100 个变异输入中 317 个未包装异常（140 BLAST、157 PAF、20 GBK）；其中 BLAST/PAF 主要是数字转换错误（`int()`/`float()`），同样是边界未包装。

**影响**

这些异常最终被 `dispatch` 兜底成 `UNEXPECTED`（exit 1），而不是 `INVALID_GFF`/`INVALID_GENBANK`（exit 5）。对 agent 来说，错误码从“输入坏”变成“工具崩了”，无法做正确的重试/修复决策；同时 `assert` 在 `python -O` 下会被移除，行为不一致。

**建议**

- 解析器里所有“用户可控输入”的 `assert` 改成显式 `raise _fail(...)`；所有索引访问先做长度/字段检查。
- GFF 在构造 `GeneFeature` 前检查 `start >= 1`、`start <= end`，或捕获 `ValidationError` 转为 `_fail("INVALID_GFF", ...)`。
- BLAST/PAF 的数字转换包 `try/except (ValueError, TypeError)`，转成 `BLAST_PARSE_FAILED`/`PAF_PARSE_FAILED`，并在 context 里带原始行。
- 把本次 fuzz 的最小样例加入 regression corpus；长期引入 `hypothesis`/property-based 测试（仅 dev 依赖）。

---

### P2-3 控制字符 / 分隔符注入破坏 TSV 的 15 列契约与 Markdown 结构

**证据（端到端）**

用 `db build` 建一个 gene 为 `gene%09tab` 的库，再 screen：

```text
$ gapit screen contig.fa --db tabdb --datadir dd --noheader
...	+	gene	tab	1-160/160	...
columns=16
```

15 列契约被破坏。原因链：

- `dbcodec.decode_seqid` 会把 `%09`/`%0A` 解码成真实的 tab/newline（`src/gapit/dbcodec.py:25-39,104-108`）。
- `db_build_ops._to_record` 把解码后的 `header.gene/accession/function` 直接写入 `Record`（`src/gapit/db_build_ops.py:60-77`）。
- `dbbuild._check_record` 只拒绝 product 中的换行和空序列，不检查 gene/accession/function（`src/gapit/dbbuild.py:79-94`）。
- `formats/tsv.py` 直接 join 这些字段，无转义。
- Markdown 的 `md_cell` 只转义 `\` 与 `|`，不处理换行（`src/gapit/formats/md.py:32-35`）；frontmatter 的 `db: {params.db}` 也未做 YAML quoting。

**影响**

- TSV/CSV 列错位，agent 按列名取值会拿到错误字段。
- Markdown 表格被换行截断；frontmatter 可被 db 名/文件名中的 `: `、`#`、换行破坏。
- 对自定义 DB（`db build`）与第三方/上游内容都可能触发；不是纯理论。
- 与 P1-2 同根：输出层没有“字段安全”抽象。

**建议**

- 在 `db build`/provider 入库边界统一校验：gene/accession/function 不得含 `\t\r\n` 及 C0 控制字符；product 不得含换行；db 名不得含控制字符。
- 若必须保留任意字符，TSV 用 `\t`/`\n` 转义约定，CSV 用标准 writer，Markdown 用 `md_cell` 处理换行（转 `<br>` 或空格）并给 frontmatter 加 YAML quoting。
- 增加“字段含 tab/newline/comma/quote”的表驱动测试，覆盖 report/reads/summary/cluster/typing 五个渲染器。

---

### P2-4 `make_blast_db` 先删除旧索引，再检查/运行 makeblastdb

**证据**

`src/gapit/db.py:109-124`：

```python
for index_file in sequences_path.parent.glob(f"{sequences_path.name}.[np]??"):
    index_file.unlink()
argv = [...]
result = run_tool(argv)  # 这里才可能抛 MISSING_DEPENDENCY
```

**影响**

- 机器上缺 `makeblastdb`（或二进制不可执行）时，`setupdb`/`--force` 会先把一个原本可用的 DB 索引删掉，然后失败。
- 与 P1-1 叠加时，目录可能既没有索引、又有旧 manifest，形成“已安装但不可用”的状态。
- 索引构建失败也不会清理半成品 `.nin/.pin`；`gapit db list` 走 `db_list_entries`，并不调用 `blastdbcmd -info`，因此半成品/损坏索引在列表里仍会显示为已安装。只有 `setupdb` 路径的 `list_databases` 会再调 `blast_db_info`。

**建议**

- 先 `shutil.which("makeblastdb")` 或 `ensure_blast()` preflight，再删除旧索引。
- 索引构建到临时前缀（`sequences.tmp.*`），成功后替换；失败清理。
- `db list` 增加“索引可被 `blastdbcmd -info` 读取”的校验，而不仅看 `.nin/.pin` 是否存在。

---

### P2-5 reads 模式把完整 PAF 物化成 `list[PafRecord]`，与“流式/不物化”声明矛盾

**证据**

- `src/gapit/minimap2_run.py:21-67` 的 docstring 写 “streaming PAF rows off stdout as lines arrive (no full-PAF materialization)”，但函数返回 `list[PafRecord]`，并在 `run_minimap2` 中 `rows.extend(...)` 汇总所有 lane。
- `src/gapit/reads.py:215-230` 又做一次 `filter_alignments` 列表、一次 `apply_gene_floors` 列表；峰值会同时持有原始/过滤后的行。
- 本机 `tracemalloc`：100k 行（500bp `cs`）= **150.9 MB**，约 1.5 KB/行；1,000 万行就是 ~15 GB 量级。

**影响**

- 对“mass screening”使命，raw reads 的 PAF 行数很容易到百万/千万级；当前实现会在聚合前把全部对齐留在内存，且 `--jobs` 会成倍放大。
- `cs` 字符串是主要内存放大项；`PafRecord` 还持有 15 个 Python 字段。
- 文档/注释声称的流式保证不成立，会误导容量规划。

**建议**

- 把 `run_minimap2` 改为生成器/回调，按 lane 逐行送入聚合器；`aggregate_coverage` 可增量维护 per-target interval union（用排序结构或分块合并），`filter_alignments`/`apply_gene_floors` 改为逐行过滤。
- cluster 引擎同理：按 locus 分组增量累积 `CsRun`，避免先 `rows` 全量。
- 如果短期不改，至少修正 docstring，并在 docs 标明内存是 O(alignments)，给出经验上限。
- 增加一个 1M 行级别的 memory regression（可标记 slow）。

---

### P2-6 三个 provider 使用明文 HTTP，且下载内容没有完整性校验

**证据**

- `src/gapit/providers/bacmet2.py:56`：`http://bacmet.biomedicine.gu.se/...`
- `src/gapit/providers/victors.py:104-105`：`http://phidias.us/...`
- `src/gapit/providers/vfdb.py:92`：`http://www.mgc.ac.cn/...`
- `src/gapit/providers/common.py:84-108`：`download_file` 只做 URL 打开 + 流式写入，没有 SHA256、签名、内容类型/大小上限。

**影响**

- 任何能 MITM 的网络位置都可以替换下载内容，进而改变基因序列/注释/分型结果；这对科研结论是完整性问题，不是“只是下载慢”。
- 没有下载上限，恶意/故障上游可填满磁盘；ZIP 解压也没有解压后大小上限（ZIP bomb）。
- 标准库 ZIP 解压会清理绝对路径/`..`，未发现直接路径穿越 RCE；主要风险是完整性与资源耗尽。

**建议**

- 能上 HTTPS 的源优先 HTTPS；确实只能 HTTP 的源，在 provider 元数据里标注 `insecure` 并在 stderr 警告。
- 对每个 provider 发布快照 SHA256（或签名），`fetch` 时校验；至少对版本化 URL（`megares_v3.00.zip`）做固定 hash。
- `download_file` 加 `Content-Length`/累计字节上限；ZIP 解压前检查成员总解压大小。
- 文档中明确“fetch 内容视为上游信任边界”，不要把 `db fetch` 当成可复现构建。

---

### P2-7 `gapit.cluster/2` 在 0.5.5→0.5.6 的 patch 版本中被破坏性重塑

**证据**

- `CHANGELOG.md:12-22` 明确写：0.5.6 “`gapit.cluster/2` TSV/CSV reshaped … replaces the packed PHENOTYPES cell and the global-best columns … Amendment of a contract shipped in 0.5.5 days earlier”。
- `src/gapit/formats/cluster2.py:74-78` 的 schema 仍是 `gapit.cluster/2`，没有 bump 到 `/3`。
- AGENTS §10.2 要求“Public contract changes (JSON schema, exit codes, CLI flags) require a schema version bump and a note in PLAN.md”。
- `pyproject` 版本 0.5.5→0.5.6 是 patch bump；语义化版本下 patch 不应含破坏性契约变更。

**影响**

- 任何在 0.5.5 消费 `/2` 的脚本/agent 在 0.5.6 会静默拿到不同的 TSV 列和 JSON 字段（虽然新增了 `winner_*`，但删除了 `BEST_LOCUS`/全局证据列）。
- 项目以“agent-first / versioned contracts”为卖点，这种先例比单个字段更伤信任。

**建议**

- 若已有下游（哪怕“no downstream adoption”是当前判断），应发布 `/3` 并保留 `/2` 一个周期；或明确把 `/2` 标为“0.5.x 实验性，不承诺跨 patch 稳定”。
- 更新 AGENTS 的“schema version bump”规则，明确 patch/minor 的边界。
- 在 CI 中加入 schema snapshot 测试：`gapit schema cluster2` 的字段集合变化必须显式更新 golden。

---

### P3-1 模块规模、函数复杂度与长参数列表

**度量**

- 10 个模块超过 AGENTS 的 250 行目标：`cmd_screen.py` 340、`gbfeatures.py` 316、`mcp_tools.py` 300、`db_ops.py` 298、`typing_models.py` 295、`typing_engine.py` 281、`screening_reads.py` 276、`formats/cluster_md.py` 264、`formats/md.py` 255、`cluster.py` 254。
- 9 个函数 C901 > 10；最高：`run_screen_reads_positional` 28、`screen_command` 26、`perform_build` 25、`screen_command.run` 25、`parse_genbank_features` 24、`run_screen` 23。
- 15 个函数超过 50 行；最长 `screen_command` 163 行、`run_screen_reads_positional` 133 行、`cmd_screen.run` 127 行、`run_screen` 124 行、`perform_build` 123 行。
- `run_screen` 19 个参数、`run_screen_assemblies` 18 个、`run_screen_reads_positional` 17 个；其中大量参数是 `float/int/bool`，类型检查器无法防止位置传参错位。

**建议**

- 为 `screen` 引入 `ScreenRequest`/`ScreenOptions` dataclass（或 pydantic 模型），CLI 只构造请求；use-case 接收请求对象。
- 把 `screen_command` 的 4 分支路由抽成 `resolve_mode()` + `run_*` 函数；把 click 装饰器保留在薄壳。
- 把 `screening.py` 的“批处理 + 进度 + 并行 + chunk sink”抽成 `BatchRunner`，供 gene/reads/cluster 三条路径复用（见 P3-2）。
- 按 AGENTS 的 250 行目标，把 `cmd_screen`、`gbfeatures`、`mcp_tools`、`db_ops` 再拆一层。

### P3-2 批处理/流式/进度逻辑在三处重复

`src/gapit/screening.py:172-195`、`src/gapit/screening_reads_positional.py:117-150`、`src/gapit/screening_cluster.py:169-175` 各自实现其中一部分：gene/reads 两条路径复制了完整的 `iter_reports()` + `ThreadPoolExecutor.map` + `tracked()`，三条路径都复制了 `chunks`/`sink` + `emit` 的累积转发，`screen_progress` 的 describe/advance 分支也重复出现。

这是 50–80 token 级别的真实重复，且“懒 header”“json 单文档”等细节已经在三处出现分歧。建议抽一个 `run_batch(items, screen_one, render, fmt, emit, quiet)`，把格式差异做成策略对象；用参数化测试锁住三条路径的字节一致性。

### P3-3 分层倒置：`records` → `formats.json`，靠延迟导入规避循环

- `src/gapit/records.py:18` 从 `gapit.formats.json` 导入 `ToolDocument`。
- `src/gapit/formats/json.py:16-18` 又导入 `hits`/`report`，形成 `db -> records -> formats.json -> report -> hits -> db` 的潜在环。
- `src/gapit/db.py:190,209` 因此用函数内 `from gapit.records import read_manifest` 规避，注释也承认这是“cycle hazard”。

**建议**：把 `ToolDocument`、`utc_timestamp`、`__version__` 相关的最小元数据移到叶子模块（如 `gapit/meta.py`）；`records` 只依赖叶子，`formats` 依赖 `records`，消除延迟导入。这样 `db.py` 可以在模块顶层导入 `records`，分层更自然。

### P3-4 JSON 表面不统一：部分输出没有 schema/未注册

- `gapit db search --json` 输出 JSONL，每行是 `SearchEntry`，**没有顶层 `schema`**（`src/gapit/db_query_ops.py:83-92`、`src/gapit/cmd_db_search.py:74`）。
- `db fetch`/`db build`/`db install` 的 receipt 也是无 schema 的裸 JSON（`db_ops.py:45-51`、`db_build_ops.py:46-57`、`cmd_db_install.py:29-34`）。
- `gapit.dblist/1`、`gapit.dboutdated/1` 有 schema 字段，但没有注册进 `SCHEMA_MODELS`（`formats/schemas.py:23-36`），所以 `gapit schema dblist` 不存在。
- 错误码是公共 API，但 `gapit schema error` 只有 envelope 形状，没有 `code` enum；agent 无法从二进制自描述全部错误码。

**建议**：给所有 stdout JSON 一个顶层 `schema`（receipt 可用 `gapit.receipt/1` + `kind`），把 `dblist`/`dboutdated`/`receipt` 注册进 schema registry；考虑从错误类生成 `ErrorCode` enum 并输出到 `gapit schema error` 的 `code.enum`。

### P3-5 MCP 表面与 CLI 不完全对齐

- `_initialize` 直接把客户端请求的 `protocolVersion` 回显（`src/gapit/mcp.py:70-77`）；若客户端发一个不支持的版本，服务器应返回自己支持的版本，而不是回声。
- 输入 schema 命名不一致：`mergeFragments`、`minGeneCov`、`minGeneId`、`minClusterCov` 是 camelCase，其他是 `min_breadth`/`min_identity`/`read_type`（`mcp_schemas.py:42-56`）。
- `db_list` 没有 `datadir` 参数，`_tool_db_list` 直接 `config.resolve_datadir(None)`（`mcp_tools.py:238-239`），无法列出自定义 datadir；若传入 `datadir` 会被静默忽略。
- `db_build` 没有暴露 CLI 的 `typing`/`floors`/`kind`，agent 无法复现 typed DB 构建（`mcp_schemas.py:103-118`）。
- `screen`/`screen_reads` 的 `format` enum 不含 `csv`，与 CLI 的四种格式不一致。

**建议**：统一 snake_case；把 `datadir` 加到 `db_list`；为 `db_build` 增加可选 `typing`/`floors`；在 `_initialize` 中校验协议版本并在不支持时返回 supported version；用协议测试覆盖。

### P3-6 错误分类/退出码在 DB 元数据上不一致

- `records.read_manifest`/`read_records` 对 malformed 内容抛 `InputError`（exit 5），例如 `src/gapit/records.py:142-157`、`66-77`。
- 但 `typing.json` malformed 用 `DatabaseError`（exit 4，`typing_models.py:245-269`），`floors.json` malformed 也用 `DatabaseError`（`gene_floors.py:53-71`），`dbcodec` header 用 `DatabaseError`。
- `db_query_ops._utc_timestamp` 对 manifest 的 `fetched_at` 又抛 `InputError`（`db_query_ops.py:167-184`）。

**影响**：用户/agent 无法用 exit code 稳定区分“我给的输入坏了”与“安装的 DB 坏了”；同一个 manifest 在不同命令下可能 4 或 5。

**建议**：把“已安装数据库的内部元数据损坏”统一为 `DatabaseError`（exit 4），用户提供的 `--tsv`/report/typing 输入保持 `InputError`（exit 5）；补一张 error-code → exit-class 的表并测试。

### P3-7 性能：启动、重复读取与不必要的依赖

- `gapit --version` ≈0.42 s，`import gapit.cli` ≈377 ms；根因是 `cli.py` 急切导入所有 `cmd_*`、`providers`（19 个模块）和 `mcp`。对 agent 高频调用，建议 lazy registry / 延迟导入 `cmd_*`。
- `reads.screen_reads` 每个 sample 都 `dict(iter_fasta_headers(database.sequences_path))`（`reads.py:228`），大 DB × 多样本会重复读整个 FASTA header；可缓存在 `Database`/run 级。
- `blast.run_blastn` 用 `subprocess.run(capture_output=True)` 把完整 BLAST 输出读入内存，再构造 `BlastRow` 列表；比 reads 的 PAF 问题轻，但同样是 O(output) 内存。
- `build_database`/`build_cluster_database` 总是调用 `tool_version_line(["minimap2","--version"])`（`dbbuild.py:234-237`、`clusterbuild.py:171-172`），因此纯 gene/contig 用户也被硬性要求安装 minimap2；cluster build 又总是建 BLAST 索引（`clusterbuild.py:159`），而 cluster screen 只用 minimap2。建议把 optional tool version 做成 `""` 容错，并让 cluster build 的 BLAST 索引可配置/延迟。

### P3-8 测试与 CI 的改进点

- **无覆盖率门禁**：仓库没有 `pytest-cov`/`coverage` 配置，CI 不收集覆盖率；无法知道 1183 个测试覆盖了哪些分支。建议至少加 `--cov=gapit --cov-report=term-missing` 并设一个可逐步提高的 baseline。
- **无属性/模糊测试**：解析器边界（BLAST/PAF/GBK/GFF/CSV）是 bug 高发区；本次 fuzz 已找到未包装异常。建议把最小样例加入 regression，长期引入 `hypothesis`（dev-only）。
- **测试对实现细节敏感**：progress 测试依赖 rich 的 `export_text` 帧行为；时间戳测试用正则；这类测试应改为断言公开行为（当前 description/计数/退出码），或用可注入的 renderer。
- **无“干净环境”CI**：缓存命中会掩盖 lock/依赖问题；建议 nightly/发布前跑一次 `pixi install --locked` 的 clean job。
- **release workflow 不跑测试**：`.github/workflows/release.yml` 在 tag 上直接 `python -m build` + 发布，没有 `needs`/test/`twine check`；建议把 gates 作为 reusable workflow 串到 release，或至少 build 后安装 wheel 做 smoke test。
- **`fmt -- --check` 不可移植**：`.github/workflows/ci.yml:26` 使用 `pixi run ... fmt -- --check`，本机 pixi 0.65 会 exit 2；`pixi run fmt --check` 正常。建议新增 `fmt-check = "ruff format --check"` task，CI 用它。

### P3-9 文档、版本与仓库卫生漂移（最新 release 尤其明显）

**最新 bundle 未同步到文档**（v0.5.6 实际有 **7** 个 bundle、**23** 个 catalog 名称）：

- `README.md:66,88,197`、`docs/installation.md:13`、`docs/faq.md:163`、`docs/databases.md:130`、`AGENTS.md:186` 仍说 “six bundles”。
- `docs/databases.md:248` 称 `vp_ok` 是 “sixth bundled database”；`docs/typing.md:166` 称 `lm_doumith` 是 “sixth”。在 7 个 bundle 的当前树上，这些 “sixth” 说法已不成立。
- `README.md:84`、`docs/agents.md:66`、`docs/quickstart.md:50`、`docs/databases.md:260` 仍说 “twenty-two databases”；`gapit db list --json` 实测 23。
- `docs/databases.md` 的 bundled 表没有 `vc_serogroup` 行；`CHANGELOG.md` 没有 `vc_serogroup`/ctxA 条目（feature commit `84cc51a` 只在 commit message 里说明）。
- `vc_serogroup/bundled.json` 没有 `license` 字段（`vp_ok` 有）；其余多数 bundle 的 `bundled.json` 也没有 license。`db list --json` 对 bundled 行也不回填 license——`_bundled_entries` 没有传 `database.metadata.license`（`db_ops.py:196-220`）。对“许可证合规”是卖点的项目，这属于机器可读 provenance 缺口。

**其他版本/文档漂移**：

- `README.md:183` 说 “1024 offline tests”；实际 1183。
- `recipe/meta.yaml:37` 说 “776 offline tests”，并称 “sdist does not ship tests/”；实际 sdist 含 279 个 tests 文件。
- `recipe/meta.yaml:5-6` 明确保留 0.3.1 的 sdist SHA；bioconda 提交前会校验失败（作者已注明 defer，但这是 release readiness 风险）。
- `docs/agents.md:19` 说 “Eleven documents are introspectable”，实际 `SCHEMA_MODELS` 有 12 个（漏了 `cluster2`）。
- `docs/outputs.md:494`/`outputs.zh.md:476` 说 error `code` 是 “snake_case”；实际是 `SCREAMING_SNAKE_CASE`（`USAGE_ERROR`、`INPUT_NOT_FOUND`）。
- `AGENTS.md:66` 仍写 “typer entrypoint”，与同文 §2/§3 的 click 迁移矛盾；`PLAN.md:10,23` 仍写 typer/any2fasta。
- `stderr.txt` 被 git 跟踪（commit `1f1ef46`），内容是一次本地 smoke 的 pixi 输出（含 `/tmp/opencode/...`），并被打进 sdist。建议删除并加 ignore。
- `src/gapit/data/typing/vp_ok.json` 与 `src/gapit/data/dbs/vp_ok/typing.json` 字节完全相同（sha256 `60bcc9c6...`），但 runtime 只使用后者；前者仅被 docs/scripts 引用，属于可漂移的重复资产，建议单一来源或删除。
- 缺少 `SECURITY.md`、`CONTRIBUTING.md`、`CODE_OF_CONDUCT.md`、issue/PR 模板、`CITATION.cff`；对 bioinformatics 工具，引用文件尤其有价值。单作者 78 commits 也意味着 bus factor 1，建议在文档中明确维护/响应预期。

### P3-10 其他边界问题

- `fragments.merge_fragment_hits` 只按 `gene` 分组（`fragments.py:125-129`），同一基因名但不同 database/accession/s_len 的 partial 会被合并；`_merged_hit` 用 `ordered[0].s_len` 作分母，若 fragments 来自不同长度记录会算错。建议分组键至少加 `database + accession + s_len`。
- `_merged_hit` 的 `identity_pct = sum(identity * aligned_len) / sum(aligned_len)`（`fragments.py:96-97`）在全部 `aligned_len == 0` 时会 ZeroDivisionError；建议显式退化。
- `cmd_screen_reads_args._READS_SUFFIXES` 包含 `.fastq.gzip`，但 `readpairs._COMPRESSION_EXTENSIONS` 只认 `.gz/.bz2`，所以 `.gzip` 文件会进入 reads 模式却无法按 sample 配对（不会崩，但会误判 single-end）。建议统一后缀表。
- `gapit.report/1` 的 `params` 不记录 `merge_fragments`（`formats/json.py:66-73`），消费方只能从 `merged: true` 行推断；若要求可复现，建议加 additive 字段。
- `parse_paf_row` 用 `model_construct` + 手写 `int()`（`paf.py:75-91`），数字损坏会 `ValueError`；`parse_blast_row` 同样（`blast.py:61-83`）。虽属“外部工具输出”，但与项目“typed boundary”目标不一致，建议包装。

---

## 6. 架构与可维护性评价

### 6.1 做得好的部分

- **职责边界的命名很诚实**：`proctools` 只做 argv，`fasta` 只做流式读，`seqconvert` 只做归一化，`dbcodec` 只做 header 编解码，`minimap` 只做 15 字符 bar 算术。读者能快速定位。
- **纯函数核心**：`hits.process_rows`、`minimap`、`cluster_math`、`typing_engine`/`typing_decide`、`fragments` 都是纯函数，测试无需 I/O；这是复杂生物信息逻辑最需要的结构。
- **契约先行**：`SPEC.md` 把 abricate 的 quirk（去重键、minus 交换、`%COVERAGE` 未舍入过滤、bar 的 1-based 伪影）写成可执行测试；新 surface（reads/cluster/typing）各自有文档化的 schema/版本。
- **错误类型分层**：`GapitError` 基类 + 4 个 exit-code 子类 + `render_error`，在语义错误路径上非常一致。

### 6.2 主要结构性债务

1. **命令层过重**：`cmd_screen.py` 同时做 click 声明、mode 路由、flag 互斥、输出交付；`screen_command` 的 4 分支嵌套在一个 127 行的 `run()` 里。它是当前最可能因新 flag 而回归的模块。
2. **长参数列表**：`run_screen`/`run_screen_reads_positional`/`run_screen_assemblies` 的形参大量是同型 `float/int/bool`，调用点靠位置对齐；新增一个阈值时很容易插错位置。`ScreenRequest` 是低风险高收益的重构。
3. **批处理逻辑重复**：见 P3-2；三条路径的“流式/顺序/进度/JSON 单文档”语义应共享一个 runner，否则未来新增第四种引擎会复制第三份。
4. **分层倒置**：见 P3-3；`records` 作为持久化模型不应依赖 presentation 的 `ToolDocument`。
5. **`gbfeatures.py` 是“解析器 + 模型 + writer”三合一**（316 行）：模型（`GeneFeature`/`LocusFeatures`/`FeaturesDocument`）被 GFF/typing/cluster 多处使用，建议拆到 `feature_models.py`，让 GFF/GBK 解析器只依赖模型。
6. **`mcp_tools.py` 的约 9 个参数 helper + 9 个 handler** 接近 300 行；如果 MCP 继续扩展，建议按 tool 分组，或把参数校验抽成通用 `Arguments` 类。

---

## 7. 测试体系评价

**强项**

- 1,183 个测试项 / 23.7k 测试 LOC，测试:源码 ≈1.8:1；对 0.5.x 项目非常充足。
- 纯函数、CLI、provider、MCP、golden、parity、streaming、concurrency 都有专门文件；“provider 真实性”测试只断言 URL/transform，不发网络请求，离线可跑。
- 61 个 golden + 8 个 parity 资产，说明“字节契约”不是口头承诺。
- 4 个失败是断言与实现语义不一致，而不是大范围失效；修复成本低。

**缺口**

- 没有覆盖率数据/门禁；无法回答“哪些错误分支没测”。
- 缺少对“输出字段含分隔符/控制字符”的对抗测试（P1-2/P2-3 的直接原因）。
- 缺少“强制重建失败”的 fault injection（P1-1）。
- 缺少“click 解析错误”的 envelope 测试（P2-1）。
- 缺少“畸形 GBK/GFF/BLAST/PAF” 的 fuzz/property 测试（P2-2）。
- progress 测试依赖 rich 内部帧行为，属于实现细节耦合；应改为断言可观察行为。
- 9 个 skip 依赖本地大样本/真实基因组；CI 永远覆盖不到这些路径，建议在 nightly 或 release 前跑一次。

---

## 8. 安全与供应链评价

**已做对的**

- 无 `shell=True`、无 `eval/exec/pickle/yaml.load`；所有外部命令 argv list。
- `minimap2` 路径强制 absolute，避免 option injection；`db build` 名称拒绝路径分隔符。
- `db install` SHA256 + 原子替换；`download_file` `.part` + `os.replace`；`generate_sequences` 临时文件 + replace。
- ZIP/tar 解压没有直接路径穿越 RCE（stdlib 清理 + 只取指定 tar 成员）。
- 许可证边界清晰：GPL/非商业内容不入 wheel；bundle 有 30 MiB 预算测试。

**主要风险**

1. **明文 HTTP + 无 checksum**（P2-6）：三个 provider 的内容可被 MITM 替换；这是最重要的供应链问题。
2. **无下载/解压大小上限**：可被恶意/故障上游用于磁盘耗尽。
3. **输出注入**（P1-2/P2-3）：不是 RCE，但会让下游 CSV/YAML/表格解析器错误解释数据；对 agent 管道是“数据完整性”问题。
4. **写路径 OSError 未统一**：只读 datadir、磁盘满、权限拒绝等可能以 `UNEXPECTED` (exit 1) 出现，而不是 `DB_*`/`FILE_IO_ERROR`；建议在 `write_records`、`db_dir.mkdir`、`generate_sequences` 等写入点包 `OSError`。
5. **无 SECURITY.md**：对“下载并解析第三方内容”的工具，建议明确信任边界、漏洞报告渠道、`db fetch` 的 MITM 风险。

---

## 9. 性能与可扩展性评价

| 场景 | 现状 | 风险 |
|---|---|---|
| CLI 启动 | ≈0.42 s，急切导入 19 providers + 所有 cmd + mcp | agent 高频调用时累积明显；可 lazy import |
| 小 contig 筛查 | BLAST 子进程 + 全量 stdout | 小数据很快；大数据 O(output) 内存 |
| reads 筛查 | PAF 全量 `list[PafRecord]`，100k 行 ≈151 MB | 百万/千万行会 OOM；与“streaming”声明不符 |
| `--jobs` | 线程池 + 每 worker BLAST/minimap2 | 内存/CPU 成倍；有 oversubscription 警告，无内存预算 |
| DB header 读取 | 每个 sample 重读整个 `sequences` header | 大 DB × 多样本重复 I/O |
| `db list` | registry 行读 manifest；本地目录读 FASTA/索引后缀 | 多个 DB 时仍较快，但不验证索引是否可读 |
| 构建 | `generate_sequences` 流式写、`verify_sequences` 再读、`count_records` 再读 | 3 次 I/O；大 DB 构建慢，但正确性优先 |

**建议**：优先把 reads PAF 改成增量聚合；其次 lazy import / 缓存 DB header；最后再考虑 BLAST 输出流式解析。对“mass screening”定位，PAF 内存是第一优先级。

---

## 10. 文档、契约与发布工程评价

**强项**

- `SPEC.md` 是少见的“可执行契约”文档：公式、边界、quirk、错误消息都有出处。
- `AGENTS.md` 把环境、目录、约定、测试策略、许可证政策写清楚；对 agent/新 contributor 都友好。
- `docs/` 双语 + MkDocs strict + CI docs workflow；README 有安装、快速开始、数据库表、输出契约、开发任务表。
- `CHANGELOG.md` 详细到“为什么改”，不是只列 commit。

**问题**

- **计数漂移**（P3-9）：7 个 bundle/23 个 catalog/1183 tests/12 schemas 的最新事实没有同步到 README/AGENTS/docs/recipe。
- **契约版本纪律**（P2-7）：`cluster/2` 在 patch 中破坏性变更，与“schema version bump”规则冲突。
- **错误码文档**：`snake_case` 描述错误；错误码没有机器可读 enum。
- **sdist 内容与 recipe 描述不符**：recipe 说 sdist 不含 tests，实际含 279 个测试文件；stale SHA 会在 bioconda 提交时失败。
- **发布流程**：tag 直接发布，没有测试/构建 smoke；没有 `twine check`/wheel 安装验证。
- **社区文件**：缺 `SECURITY.md`/`CONTRIBUTING.md`/`CITATION.cff`/issue 模板；单作者 bus factor 1。

---

## 11. 改进路线图（建议）

### 立即（0–3 天，P1）

1. 修 4 个 progress 测试：断言当前 description + 计数，或让 `describe` 保留 `"Screening "` 前缀；在干净 env 重跑全量。
2. 用 `csv.writer` 替换 5 个 CSV 渲染器；输入侧 `summary`/`typing_input` 改用 `csv.reader`；补逗号/引号/换行/paired-end 测试。
3. 强制重建事务化：staging 目录 + 原子切换，或旧 manifest tombstone + 恢复；`make_blast_db` 先 preflight 再删索引；加 fault-injection 测试。
4. click usage error 统一 envelope（未知 option/缺 required/非法值），补测试。

### 近期（1–2 周，P2）

5. 修 GBK/GFF 解析器边界（去 assert、长度检查、`ValidationError` 转换）；把 fuzz 最小样例加入 regression。
6. 在 `db build`/provider 入库边界拒绝 gene/accession/function 中的 `\t\r\n`/控制字符；Markdown 处理换行/frontmatter quoting。
7. reads PAF 增量聚合（或至少限制/分批），修正“no full-PAF materialization”文档；加 1M 行 memory 测试。
8. 三个 HTTP provider：HTTPS 优先、SHA256 pin、下载/解压上限、stderr 警告。
9. `cluster/2` 契约决策：发布 `/3` 或明确标为实验性；补 schema snapshot 测试。
10. 修复文档计数（7 bundles/23 DB/1183 tests/12 schemas）、error code 大小写、AGENTS typer、recipe SHA/sdist 描述；删除 `stderr.txt`，单一化 `vp_ok` 资产。

### 1.0 前（架构/工程）

11. `ScreenRequest` + `BatchRunner` 重构；拆分 `cmd_screen`/`gbfeatures`/`mcp_tools`；消除三份批处理重复。
12. 把 `ToolDocument`/`utc_timestamp` 移到叶子模块，消除 `records → formats.json` 倒置和 `db.py` 延迟导入。
13. 统一所有 stdout JSON 的 schema（含 receipt/search），注册 `dblist`/`dboutdated`；错误码 enum 进 `gapit schema error`。
14. 覆盖率门禁（先记录 baseline，再逐步提高）；nightly clean-env job；release workflow 串 gates + wheel smoke test。
15. 补 `SECURITY.md`/`CONTRIBUTING.md`/`CITATION.cff`/issue 模板；明确 bus factor/维护预期。
16. 可选：lazy provider imports、DB header 缓存、cluster build 的 BLAST 索引可配置、gene-only 不强依赖 minimap2。

---

## 12. 附录：方法与局限

### 12.1 本次实际执行的操作

- 阅读 `pyproject.toml`、`pixi.toml`、CI/release/docs workflow、README/AGENTS/SPEC/PLAN/CHANGELOG 与全部 `src/gapit` 模块（96 个）。
- 运行 `pixi run lint`、`pixi run typecheck`、`pixi run fmt --check`、`pixi run fmt -- --check`、`pixi run test`、`pytest --collect-only`。
- 在 py311 环境单独运行 `tests/test_progress.py` 复现 4 个失败。
- AST 扫描：模块/函数/LOC/复杂度/注解/导入图；`ruff --select S/C901/TRY/BLE` 补充扫描。
- 定向复现：CSV 列错位（screen/summary/paired reads）、TSV 控制字符注入、force 重建后 manifest/sequences 不一致、GBK/GFF 畸形输入异常、`make_blast_db` 先删索引。
- 模糊测试：BLAST/PAF/cs/header/GBK/GFF/FASTA 七类各 300 个变异输入。
- 微基准：PAF 100k 行内存、BLAST/PAF 解析吞吐、CLI 启动时间。

### 12.2 局限

- 未执行 `pixi run -e parity`（parity env 未安装，需要额外 bioconda solve），因此未独立复核 abricate 字节 parity；本报告对 parity 的评价基于代码、测试资产与文档。
- 未在 macOS/Windows 上验证；项目 classifier 含 macOS，但 CI 只跑 Linux。
- 未做专业安全审计/渗透测试；HTTP/供应链结论基于静态代码与威胁建模。
- 无覆盖率工具（`pytest-cov` 不在环境），因此没有量化覆盖率；测试评价基于测试结构、断言密度与抽样阅读。
- 复杂度为 AST 决策点近似（非 mccabe 精确值），用于相对排序而非绝对阈值。
- 内存基准为合成 500bp `cs` 行，实际 PAF 行长度分布会不同；数量级结论成立，绝对数字仅作参考。

---

**最终结论**：gapit 的“工程质量底盘”很扎实，类型、测试、文档、许可证与确定性设计都值得保留；当前最需要补的是**失败原子性**（manifest/staging）、**输出层转义**（CSV/TSV/Markdown）与**干净环境门禁**。这三类问题修完后，它会从“高质量 0.5.x 项目”进入“可以放心交给 agent/流水线长期运行”的区间。
