# 脑膜炎奈瑟菌血清群分型（`nme_sero`）

`gapit` 提供一个按需下载的 **`nme_sero`** 数据库，只覆盖脑膜炎奈瑟菌
（*Neisseria meningitidis*）的**基因存在**层面：

- 血清群 A / B / C / W / X / Y（W/Y 双扩增时输出 W/Y）
- `ctrA` 胶囊位点标记
- `sodC` 单拷贝深度内参

这符合 gapit 的定位：**判断已知基因在不在**。MLST/CC 交给 gMLST。
PorA VR1/VR2、FetA、PorB、BAST/MenDeVAR 属于**等位基因级**序列分型，
需要蛋白等位基因库和 BAST 组合查表，刻意**不放进**这个数据库；菌株分型交给
gMLST 数据库，抗原等位基因 ID 交给外部分型流程。

## 安装

```bash
gapit db fetch -d nme_sero
gapit db list | grep nme_sero
gapit db info -d nme_sero
```

数据库在 fetch 时下载（GPL-3.0 / PubMLST 条款），不随 MIT wheel 分发。
本次验证使用的上游快照为 commit
`35c4a5bff24fbc3625f92ce575d7c2e99ecd4378`（2026-02-19）。

## 从 contig 做血清群

```bash
gapit screen isolate.fna -d nme_sero
gapit typing isolate.tsv
```

`nme_sero/typing.json` 提供 `men_serogroup` 方案，控制基因是 `sodC`；
没有 `sodC` 的样本会判为 untypeable，而不是静默给出一个血清群。

## 从 FASTQ 做血清群 + sodC 深度内参

`nme_sero/screen.json` 设置：

```json
"reads": {
  "min_identity": 80.0,
  "min_mapq": 20,
  "min_depth": 5.0,
  "min_relative_depth": 0.10,
  "relative_to": "sodC"
}
```

因此直接跑：

```bash
gapit screen -d nme_sero -1 R1.fastq.gz -2 R2.fastq.gz
```

会按“基因深度 / sodC 深度 ≥10%”判定存在。`sodC` 是保守的单拷贝
Neisseria Cu-Zn 超氧化物歧化酶，作为相对深度分母；它不是严格的种特异性
标记。可用 `--min-relative-depth` 或 `--relative-to` 临时覆盖。

## 为什么没有 `gapit meningotype` 子命令？

`gapit typing` 是基于**基因存在/缺失**规则的判定引擎；PorA/FetA/PorB/BAST
是**等位基因级**调用，需要蛋白等位基因库、精确全长肽匹配和 BAST 组合查表，
超出通用 `screen -> typing` 契约，因此不放进这个数据库。菌株分型交给
gMLST 数据库；抗原等位基因 ID 交给外部分型流程。

## `gapit typing` 能否判断 EX7E motif？

不能直接判断。`gapit typing` 读取的是筛查表（基因、覆盖度、一致性），看不到
比对上的密码子，也不能翻译 synG 的 EX7E 窗口；通用规则引擎没有 motif/等位基因
原语。

可以用**等位基因特异探针**近似：把 P、G、S 三种带侧翼的探针放进数据库，把规则的
`identity_floor` 提到约 99.5，用存在/缺失选择表型（即 `typing` 手册里的
“allele probe trick”）。但这种做法只能识别精确的探针序列，对侧翼的沉默变异敏感，
也无法报告新等位基因；而且当较长的血清群扩增子与探针重叠时，abricate 的
`-culling_limit 1` 可能遮蔽探针命中。

忠实实现应沿用上游 `menwy`：比对 `synG`，取出 919–945 窗口（坐标相对于 3114 bp
的 synG 参考序列，**不是**基因组坐标），翻译后读取第 3 位（P = W，G = Y，
S = W/Y）。这属于等位基因级序列解释，应放在专门的 finetyping 工具中，而不是这个
基因存在数据库。

因此 W/Y 判定不会对这一步沉默：规则的 `notes`（JSON）、TSV 的 `NOTES` 列和
Markdown 表都会带一条明确的 `NEXT STEP`，说明 EX7E 窗口、P/G/S 对应关系和
`menwy` 后续工具，并注明坐标是 synG 相对坐标，不能直接当作基因组坐标使用。

LLM 能理解这条指令并把样本路由到 `menwy`，但不应在没有比对后的 synG 序列或工具
的情况下自行“翻译”窗口：这一步需要比对、链方向处理和细菌遗传密码。没有 `menwy`
（或等价等位基因 caller）时，agent 的正确行为是停下并报告 W vs Y 未解决，而不是
根据 W/Y 共检出猜一个结果。通过 MCP 的 `typing` 工具也能拿到同样的信息：把筛查表
文本作为 `table` 传入，读取 `files[].phenotypes.men_serogroup.notes`。

## Caveats

- **W/Y**：当前用 W 和 Y amplicon 共检出判断；上游 synG EX7E motif 细分
  未实现，罕见 W/Y 变异株建议用 PCR 或上游 `menwy` 复核。原因见上一节。
- `ctrA` 是胶囊位点标记，不是普适毒力基因；无荚膜菌株可能缺失。
- `sodC` 是保守单拷贝内参，不是严格的种特异性标记。
- 数据库内容遵循 GPL-3.0（meningotype）+ PubMLST Neisseria 条款；按需
  下载，不随 wheel 分发。

## 验证

在上游 A/B/C/W/X/Y 测试基因组上，`gapit screen -d nme_sero` +
`gapit typing` 分别返回 A、B、C、W、X、Y；`sodC` 和 `ctrA` 在六个基因组中
均检出。
