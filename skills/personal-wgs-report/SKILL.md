---
name: personal-wgs-report
description: 用自己的人类全基因组测序（WGS，约 30× 双端 FASTQ）数据，在一台装了 Docker 的 Linux 机器上从头分析到出报告：质控、比对与变异检测（有 NVIDIA 显卡用 Parabricks，没有就走 CPU）、注释与临床相关筛查（ACMG 次要发现、隐性携带、药物基因组、重复扩增、SMA / 地贫 / CYP2D6 等难测基因、HLA）、祖源与多基因评分，可选结构变异和扩展分析（局部祖源、古 DNA、尼安德特人片段、纯合片段、血型、衰老指标、血液里的微生物 DNA），最后写一份大白话中文报告。用户说“我测了全基因组 / 拿到了 WGS 原始数据 / 帮我分析我的基因组”时使用。
---

# 个人全基因组分析与报告

这个 skill 把一套已经在真实数据上跑通的分析流程（脚本 + 说明）交给你（AI）执行。数据从头到尾都留在用户自己的机器上；
你负责：和用户确认需求 → 准备环境 → 跑流程 → 检查结果 → 按模板把结果写成用户看得懂的报告。

- 默认参数针对**东亚人群**（POP=EAS）。其他人群核心流程照样能跑；祖源相关的扩展分析要先按 `references/customizing.md` 改参数。
- 这是**研究 / 自我了解性质**的分析，不是临床诊断。报告里任何可能影响就医、用药、生育决定的发现，都要写明“需在有资质的临床实验室验证，并咨询医生 / 遗传咨询师”。

## 第 0 步：先和用户确认（不要跳过）

一次问清，能从机器上查到的自己查（`nproc`、`free -g`、`df -h`、`nvidia-smi`、`docker info`）：

1. **数据**：FASTQ 在哪、几对文件（多条泳道就有多对）、测序公司给的 MD5 文件在不在；是不是全基因组（不是外显子组 / 低深度）；样本类型（血 / 唾液）。
2. **几个人**：只分析本人，还是还有伴侣 / 家人（每个人都要本人同意；两个人一起分析会多出“同一隐性基因都携带”的生育相关结果）。
3. **机器**：建议 ≥16 核、≥64 GB 内存、≥1 TB 空闲磁盘（核心流程约 600 GB，加扩展分析再加 ~100 GB）；有没有 NVIDIA 显卡（≥16 GB 显存可用 Parabricks，30× 数据比对 + 变异检测约 1 小时；纯 CPU 约 8–15 小时）。
4. **网络**：要下载约 60–150 GB 参考数据和容器镜像；国内访问 GitHub / Google / Ensembl 慢时需要代理（写进 `config.env` 的 HTTPS_PROXY），quay.io 慢可以换镜像站（BIOC）。
5. **想看什么**：核心流程（必做）／结构变异（+3–5 小时）／扩展分析（祖源细节、古 DNA、尼安德特人、ROH、血型、衰老指标、血液微生物，+半天）。
6. **知情选择（“不知情权”）**：报告里可能出现一些人不想知道的结果。逐项问用户要不要写进报告——
   ACMG 次要发现（遗传性肿瘤、心源性猝死等可干预疾病基因）、APOE（阿尔茨海默病风险）、隐性携带、多基因风险评分、神经肌肉病基因的意义未明变异。用户说不要的，报告里只写“按你的要求未列出”。
7. **隐私告知**：原始数据和完整结果只在本机；但你（AI）读取的结果摘要会发送给你所运行的模型服务。要明确告诉用户这一点并得到同意；用户不同意就只帮他跑流程，让他自己读 `summary.md`。

## 第 1 步：建工作目录、写配置

```bash
bash <本仓库>/skills/personal-wgs-report/scripts/setup/init_workdir.sh /data/wgs      # 任意有空间的目录
cd /data/wgs && $EDITOR config.env      # 必改：SAMPLE（只用字母数字下划线，不要用真名）、FASTQ_R1/FASTQ_R2（逗号分隔，一一对应）、THREADS、MEM_GB
```

之后**所有命令都在工作目录里执行**（脚本靠当前目录找 `config.env` 和数据）。工作目录结构：`raw/ ref/ work/ qc/ logs/ results/ tools/ tmp/`，`scripts` 是指向本仓库的软链接。
第二个人：复制一份 `config.<名字>.env` 改 SAMPLE 和 FASTQ，运行任何脚本时加 `WGS_CONFIG=config.<名字>.env`。

## 第 2 步：准备环境（一次性，约 2–6 小时，大部分是下载）

按顺序执行，每一步看到 `... DONE` 再走下一步。长任务放进 `tmux` / `nohup`，日志写到 `logs/`：

```bash
bash scripts/setup/check_env.sh        # 体检：Docker、磁盘、内存、GPU、FASTQ 可读；有 [不行] 先解决
bash scripts/setup/install_tools.sh    # fastp / seqkit / bwa-mem2、Python 环境（SpliceAI、hmmix、CNVpytor…）、Cyrius / Gauchian / Yleaf、中文字体
bash scripts/setup/pull_images.sh      # 拉容器镜像（版本固定在 scripts/lib/env.sh）
bash scripts/setup/get_reference.sh    # GRCh38 分析集（no-alt + decoy + 屏蔽假重复区）
bash scripts/setup/index_reference.sh  # BWA 索引（GPU 用）或 bwa-mem2 索引（CPU 用），按 USE_GPU 自动选
bash scripts/setup/get_resources.sh    # 质控 / 注释 / 基因 - 疾病表 / HLA·KIR 索引（约 40 GB）
bash scripts/setup/prep_annotation.sh  # 建索引、注释范围、基因表
bash scripts/setup/get_resources_sv.sh   # 可选：结构变异（约 20 GB）
bash scripts/setup/get_resources_ext.sh  # 可选：扩展分析（约 45 GB）
```

REVEL、PanelApp 在部分地区下载会被拦（脚本会提示怎么手动下载）；缺了也能跑，只是少一个错义打分 / 少神经肌肉面板。

## 第 3 步：跑分析

```bash
tmux new -s wgs
bash scripts/run_all.sh [--with-sv] [--with-ext] 2>&1 | tee logs/run_all.log
```

- 每步成功后写 `logs/done/<步骤>.<样本>`，中断后重跑会自动接着做；想重做某一步就删掉对应 done 文件，或 `--only <步骤名>`。
- 每步的详细日志在 `logs/<步骤>.<样本>.log`，出错先看这里，再看 `references/pitfalls.md`。
- 时间（30×，16–24 核）：GPU 机器核心流程约 4–6 小时；纯 CPU 约 15–25 小时；`--with-sv` +3–5 小时；`--with-ext` +6–12 小时。
- 每一步做什么、输出在哪：`references/pipeline.md`。

## 第 4 步：检查结果靠不靠谱（写报告前必须做）

先读 `results/<样本>/summary.md` 第 1 节（数据质量），对照下面的标准。不合格要先告诉用户，并在报告开头说明影响：

| 指标 | 正常范围 | 不正常意味着 |
|---|---|---|
| 常染色体有效平均深度 | ≥ 20×（15–20× 可用但灵敏度下降） | 太低：杂合插入缺失、嵌合、结构变异、重复扩增容易漏 |
| ≥10× 覆盖比例 | ≥ 90% | 覆盖不均（建库 / GC 偏好） |
| 比对率 | ≥ 99%（血液样本） | 偏低：污染、样本降解、接头问题 |
| 重复率 | ≤ 30%（NovaSeq X 光学重复偏高属常见） | 文库复杂度低，有效深度打折 |
| 污染 FREEMIX | < 2% | 高：样本混了别人的 DNA，杂合判断不可信——**停下来**和用户确认 |
| 推断性别 | 与用户自述一致 | 不一致：样本搞错 / 性染色体异常，先确认再继续 |
| Ti/Tv（全基因组 SNV） | 约 2.0–2.1 | 明显偏低：假阳性多 |
| 杂合 / 纯合比 | 非洲裔约 1.8–2.0，其他约 1.3–1.7（东亚常见 1.2–1.5） | 极低：近亲婚配或污染判断错误；极高：污染 |

再快速扫一遍其他各节有没有“未运行 / 失败”，有的话看对应日志补跑，或者在报告里说明缺了哪部分。

## 第 5 步：写报告

1. 读 `results/<样本>/summary.md`（数据底稿，所有数字都从这里取）、`assets/report_template.md`（结构与每节写法）、`references/interpretation.md`（每类结果怎么解读、阈值、措辞）。
2. 需要细节时按底稿里的“来源”去读原始文件（如 `work/findings/<样本>.carrier.tsv`、PharmCAT 报告、`results/<样本>/figs/*.csv`）。
3. 把报告写成 `results/<样本>/report.md`，图用相对路径引用 `figs/*.png`（`results/<样本>/figs/index.json` 列出了已生成的图）。
4. 转成单文件网页（图片内嵌，可离线打开）：`tools/venv-ext/bin/python scripts/report/render_html.py results/<样本>/report.md`

写作要求：

- **大白话先行**：每节第一句话就是结论（“对你意味着什么”），然后才是数字和依据；术语第一次出现时用一句话解释（参考 `references/primer.md`）。用第二人称“你”。
- **只写有数据支撑的话**：数字一律来自底稿；底稿里没有的不要编；分析没跑就写“未做”。
- **分清确定程度**：“检出且 ClinVar 多星致病” ≠ “软件预测可能有害” ≠ “意义未明”。意义未明变异（VUS）不要写得像坏消息。
- **不制造焦虑，也不轻描淡写**：携带一个隐性致病等位人人都有（平均 1–5 个），主要影响生育；多基因评分只是相对高低，不是患病概率；ACMG 次要发现阳性则明确建议临床验证和就诊。
- **尊重第 0 步的选择**：用户不想看的类别不要写结果。
- **重要发现要给读段证据**：ACMG / 携带 / 双等位等变异在 `results/<样本>/evidence_igv.html` 里有读段图，报告里提一句可以打开看，并提醒低质量（qc_note 非空）的要谨慎。
- **不要给诊疗指令**：用药部分写“告知医生 / 药师，由医生决定”，不要让用户自行停药换药。

## 第 6 步：交付

告诉用户：报告在哪（`report.md` / `report.html`）、最重要的 3–5 条结论、哪些需要线下验证、哪些分析没跑或失败了、原始数据和中间文件占了多少空间（`du -sh work ref`）以及哪些可以删（`tmp/`、`work/ext/panel/chr/` 可删；BAM、gVCF 建议保留备份）。

## 隐私与伦理（必须遵守）

- 原始数据（FASTQ / BAM / VCF）和完整结果**只留在用户机器上**：不要上传到任何在线工具、网盘、第三方解读网站、公开仓库，除非用户明确要求并知道后果。
- 联网步骤只下载公共数据库；`bloodgroups.py` / `trait_snps.py` 查 Ensembl 时只发 rsID，不发个人基因型。
- 分析别人的数据必须有本人同意；未成年人的成人发病疾病（如 APOE、遗传性肿瘤）不做预测性解读。
- 报告是给用户本人的；用户想公开分享（博客、社交媒体）时，提醒他：基因组信息会暴露亲属的信息，且无法撤回；
  去掉样本编号、测序公司订单号、具体致病变异（尤其涉及家人的）后再发。
- 不要把用户的结果写进本仓库或任何会被提交的文件。

## 常用小工具（回答用户追问时用）

```bash
tools/venv-annot/bin/python scripts/interpret/quick_sites.py <样本> rs671:ALDH2 rs1229984     # 任意 rsID 的基因型（从 gVCF 取，区分“纯合参考”和“没测到”）
tools/venv-annot/bin/python scripts/interpret/pileup_check.py chr12 111803962 <样本>        # 某个位点（例：ALDH2 rs671）的读段细节：支持读段数、正反链、比对 / 碱基质量、是否重复
tools/venv-annot/bin/python scripts/interpret/pgs_search.py "Type 2 diabetes"                # 在 PGS Catalog 找适合本人群的评分
FIG_PERSON=我 tools/venv-ext/bin/python scripts/report/make_figures.py <样本>                 # 重画图（图里称呼改成“我”，适合写成第一人称文章）
```

## 出问题时

`references/pitfalls.md` 收集了实际跑过遇到的坑（Parabricks 显存 / 驱动、DeepVariant CPU 很慢、VerifyBamID 失败、PharmCAT 不认识的 CYP2D6 等位、RFC1 被误报扩增、线粒体 poly-C 伪影、hmmix 需要 bcftools 等）。
解决不了的，把出错步骤的日志最后 50 行给用户看，说明卡在哪一步、可以怎么绕过（大多数分析步骤互相独立，失败了不影响其他结果）。

## 文件地图

| 位置 | 内容 |
|---|---|
| `config.example.env` | 配置模板（`init_workdir.sh` 复制成 `config.env`） |
| `scripts/setup/` | 一次性准备：工具、镜像、参考、数据库 |
| `scripts/01_qc.sh` … `15_sv.sh` | 主流程各步骤（`run_all.sh` 按顺序调用） |
| `scripts/interpret/` | 解读：变异分级、携带、PRS、PCA、性状位点、质控汇总… |
| `scripts/ext/` | 扩展分析：局部祖源、古 DNA、古人类片段、ROH、衰老、血型、拷贝数、病毒 |
| `scripts/report/` | 画图、汇总成底稿、Markdown → HTML |
| `assets/report_template.md` | 报告结构与每节写法 |
| `assets/population/EAS.json` | 东亚人群的扩展分析参数（其他人群照此另写） |
| `assets/genelists/`、`assets/pgs_selected_EAS.tsv` | ACMG SF v3.3 基因与规则、性状位点、东亚 PRS 评分清单 |
| `references/` | 流程说明、解读指南、常见问题、定制方法、科普素材 |
