# personal-wgs-report：让 AI 帮你读懂自己的全基因组

测了全基因组（WGS），拿到几十 GB 的原始数据，然后呢？这个仓库是一个 **Agent Skill**：把它交给 Claude Code 之类的 AI 编程助手，
AI 会在你自己的电脑 / 服务器上，从原始 FASTQ 一路跑到一份**大白话的中文报告**。数据从头到尾不离开你的机器。

流程最初是作者分析自己的 30× 全基因组时写的，已在真实数据上完整跑通；这里去掉了所有个人数据，参数化成任何人都能用的形式。

## 能得到什么

**核心流程**（GPU 机器约 4–6 小时，纯 CPU 约 15–25 小时）

- 数据质量：深度、覆盖、重复率、污染、性别核对
- 比对与变异检测：有 NVIDIA 显卡用 Parabricks（BWA-MEM + DeepVariant，约 1 小时），没有就用 bwa-mem2 + DeepVariant CPU
- 和健康直接相关：ACMG SF v3.3 的 84 个可干预疾病基因、隐性病携带（两个人一起分析时给出“同一基因都携带”）、
  神经肌肉病基因、重复扩增（51 个致病位点，含 RFC1 重复单元判别）、SMA（SMN1/2）、α-地贫缺失、GBA1、结构变异（可选）
- 用药：PharmCAT 药物基因组（CYP2C19、CYP2D6（Cyrius）、华法林相关、他汀相关、巯嘌呤类…）+ HLA 用药风险等位 + 线粒体氨基糖苷致聋位点
- 多基因风险评分：默认 16 个用东亚 / 中国人群开发的 PGS，对照 1000 Genomes 东亚人群给百分位
- 性状：酒精代谢、乳糖耐受、耳垢、APOE（可选择不看）等
- 祖源：大洲判断、1000 Genomes 主成分、母系（线粒体）和父系（Y）单倍群；HLA / KIR 分型

**扩展分析**（可选，+半天）：染色体祖源涂色（局部祖源）、和古 DNA（AADR）比较（PCA、f3、qpAdm）、尼安德特人 / 丹尼索瓦人片段、
父母血缘（ROH）、红细胞血型扩展分型、拷贝数会变的基因（AMY1、LPA、UGT2B17…）、衰老指标（线粒体拷贝数、Y 丢失、克隆性造血、端粒）、血液里的微生物 DNA。

最后得到：`report.md` / `report.html`（单文件，图片内嵌）、14 张左右的图、关键变异的读段证据页（igv-reports）、PharmCAT 用药报告，以及所有中间结果。

## 需要什么

- **数据**：人类全基因组双端测序 FASTQ（Illumina 或华大，约 30×；多条泳道的多对文件都可以）
- **机器**：Linux + Docker；建议 ≥16 核、≥64 GB 内存、≥1 TB 空闲磁盘；NVIDIA 显卡（≥16 GB 显存）可选但快很多
- **网络**：需要下载约 60–150 GB 的公共参考数据和容器镜像（国内可能需要代理）
- **一个 AI 编程助手**：Claude Code、Codex CLI 等能在终端里执行命令的工具

## 用法

### 交给 AI（推荐）

```bash
git clone https://github.com/Hirriririir/personal-wgs-report.git
mkdir -p ~/.claude/skills && ln -s "$PWD/personal-wgs-report/skills/personal-wgs-report" ~/.claude/skills/personal-wgs-report
```

然后在 Claude Code 里说：

> 用 personal-wgs-report 分析我的全基因组，FASTQ 在 /data/mywgs/，机器有一块 4090。

AI 会先问清楚你的数据、机器和想看哪些结果（包括哪些结果你**不想知道**），再一步步准备环境、跑流程、检查质量、写报告。
其他支持 Agent Skills 的工具，把 `skills/personal-wgs-report` 放进它的 skills 目录；不支持的，直接让 AI 读 `skills/personal-wgs-report/SKILL.md`。

### 自己跑

```bash
S=personal-wgs-report/skills/personal-wgs-report
bash $S/scripts/setup/init_workdir.sh /data/wgs && cd /data/wgs
vi config.env                                   # 样本名、FASTQ 路径、线程、是否用 GPU
for x in check_env install_tools pull_images get_reference index_reference get_resources prep_annotation; do
  bash scripts/setup/$x.sh || break
done
bash scripts/run_all.sh                         # 可加 --with-sv / --with-ext
less results/<样本名>/summary.md                # 结果底稿；报告结构见 assets/report_template.md
```

详细说明：[SKILL.md](skills/personal-wgs-report/SKILL.md)（流程总览）、[references/pipeline.md](skills/personal-wgs-report/references/pipeline.md)（每一步）、
[references/interpretation.md](skills/personal-wgs-report/references/interpretation.md)（怎么解读）、[references/pitfalls.md](skills/personal-wgs-report/references/pitfalls.md)（常见问题）。

## 隐私

- 所有计算都在你自己的机器上；脚本只从公共数据库下载数据，查 Ensembl 时只发送 rsID。
- 让云端 AI 写报告时，AI 会读取结果摘要（不是原始数据），这部分内容会发送给 AI 服务商。介意的话可以只让 AI 跑流程，自己读 `summary.md`。
- 不要把自己的结果提交到这个仓库或任何公开的地方。基因组信息也关乎你的亲属，公开后无法撤回。

## 免责声明

本项目用于个人研究和自我了解，**不是医学诊断**，没有经过临床验证。任何可能影响就医、用药、生育决定的结果，都需要在有资质的临床实验室验证，并咨询医生或遗传咨询师。

## 默认针对东亚人群

PRS 评分、局部祖源参考组、古 DNA 地区和 qpAdm 模型默认按东亚人群设置。其他人群的改法见 [references/customizing.md](skills/personal-wgs-report/references/customizing.md)。

## 用到的主要工具与数据

NVIDIA Parabricks、BWA / bwa-mem2、DeepVariant、GLnexus、samtools / bcftools、fastp、mosdepth、VerifyBamID2、somalier、Ensembl VEP、gnomAD、ClinVar、AlphaMissense、REVEL、SpliceAI、
GenCC、ClinGen、PanelApp Australia、ExpansionHunter / stranger / REViewer、SMNCopyNumberCaller、Cyrius、Gauchian、T1K（IPD-IMGT/HLA、IPD-KIR）、mutserve、GATK、HaploGrep3、Yleaf、
PharmCAT（CPIC / DPWG）、PGS Catalog、plink2、1000 Genomes、Manta、Delly、CNVpytor、AnnotSV、Beagle、RFMix、AADR、ADMIXTOOLS、EIGENSOFT、hmmix、Kraken2、TelSeq、igv-reports。
请在使用结果时引用相应的工具和数据库；各数据库有自己的使用条款（部分仅限非商业用途）。

## License

代码以 MIT 许可发布（见 [LICENSE](LICENSE)）。下载的第三方数据和容器遵循各自的许可。

---

## English

**personal-wgs-report** is an Agent Skill that lets an AI coding assistant (e.g. Claude Code) analyse your own human whole-genome sequencing data end to end on your own Linux machine. It covers QC, GPU (Parabricks) or CPU alignment and DeepVariant calling, annotation (VEP, gnomAD, ClinVar, AlphaMissense, SpliceAI), ACMG secondary findings, carrier screening, pharmacogenomics (PharmCAT), repeat expansions, hard-to-call genes (SMN1/2, CYP2D6, GBA1, HBA), HLA/KIR typing, ancestry, polygenic scores, and optional SV calling and ancestry/ancient-DNA/archaic-introgression extras. It then writes a plain-language report (Chinese by default). Defaults target East Asian ancestry; see `references/customizing.md` to adapt. Research use only, not a medical diagnosis. Your data never leaves your machine.
