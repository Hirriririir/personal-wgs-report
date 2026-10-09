# 流程说明：每一步做什么、输出在哪

所有路径相对工作目录。`<s>` = 样本名（config 的 SAMPLE）。时间按 30× 数据、16–24 核估计。

## 准备（scripts/setup/，一次性）

| 脚本 | 内容 | 体积 / 时间 |
|---|---|---|
| `init_workdir.sh` | 建目录、`scripts` 软链接、`config.env` | — |
| `check_env.sh` | Docker、磁盘、内存、GPU、FASTQ 可读性体检 | — |
| `install_tools.sh` | `tools/bin`（fastp 1.4.0、seqkit 2.14.0、bwa-mem2 2.2.1）；Python 环境 `venv-annot`（SpliceAI 1.3.1 + TensorFlow 2.15、cyvcf2、pysam、igv-reports）、`venv-ext`（hmmix 0.9.2、matplotlib、CrossMap…）、`venv-illumina`（Cyrius、Gauchian）、`venv-yleaf`（Yleaf 4.1.5）、`venv-cnv`（CNVpytor 1.3.2）；中文字体 | ~5 GB |
| `pull_images.sh` | 容器镜像（版本见 `scripts/lib/env.sh`；GPU 用 Parabricks 4.7.1，CPU 用 DeepVariant 1.9.0） | ~25 GB |
| `get_reference.sh` | NCBI GRCh38 no-alt + hs38d1 decoy analysis set（chrY PAR 已屏蔽）+ 按 GRC_exclusions 屏蔽假重复区（GIAB v3 做法，坐标不变） | 3 GB |
| `index_reference.sh` | GPU：BWA 索引（`ref/GRCh38/bwa/`）；CPU：bwa-mem2 索引（建索引要 ~80 GB 内存） | 5–15 GB，~1 h |
| `get_resources.sh` | somalier 位点与 1000G 参考、VerifyBamID2 资源、STR 目录、ClinVar、gnomAD v4.1 约束、AlphaMissense、GENCODE v50、REVEL（可选）、GenCC、ClinGen 剂量、PanelApp（可选）、VEP 116 merged 缓存、T1K 的 HLA / KIR 索引、chrM | ~40 GB |
| `prep_annotation.sh` | ClinVar 改 chr 命名、建索引；注释范围（外显子 ±50 bp + ClinVar 位点）；基因表；蛋白编码基因 BED、CHIP 基因编码区 | — |
| `get_resources_sv.sh` | AnnotSV 3.5 注释包、Delly 排除区 | ~20 GB |
| `get_resources_ext.sh` | 1000 Genomes 3202 人（PLINK2）、AADR v66.p1 1240K、hmmix 资源 + 4 个古人类基因组、Beagle 遗传图谱、hg19→hg38 链、cytoBand、Natural Earth 底图、Kraken2 standard-8GB | ~45 GB |

## 主流程（scripts/run_all.sh 按顺序调用）

| 步骤 | 做什么 | 主要输出 | 时间 |
|---|---|---|---|
| `01_qc.sh` | MD5 核对、seqkit 统计、fastp 质控报告（不改写读段）、读段名统计仪器 / 泳道 | `qc/fastp/<s>.*.json/html`、`qc/<s>.lanes.txt` | 0.5–1 h |
| `02_align_call.sh` | GPU：Parabricks fq2bam（BWA-MEM、排序、标重复，光学重复距离 2500）+ DeepVariant（WGS 模型，gVCF）。CPU：bwa-mem2 → samtools fixmate / sort / markdup → DeepVariant CPU | `work/<s>/<s>.bam`、`<s>.dv.vcf.gz`、`<s>.dv.g.vcf.gz` | GPU ~1 h；CPU 8–15 h |
| `03_bam_qc.sh` | mosdepth 覆盖度、samtools stats、VerifyBamID2 污染、somalier 指纹、按 X/Y 深度推断性别 | `work/<s>/qc/` | ~1 h |
| `05_str.sh` | ExpansionHunter 5 + stranger（51 个致病重复位点）、REViewer 读段图、RFC1 重复单元组成 | `work/<s>/str/` | ~0.5 h |
| `06_special_loci.sh` | SMNCopyNumberCaller、Cyrius（CYP2D6）、Gauchian（GBA1）、α-珠蛋白簇读深 | `work/<s>/special/` | ~0.5 h |
| `07_hla_kir.sh` | T1K：HLA（IPD-IMGT/HLA）和 KIR（IPD-KIR）分型 | `work/<s>/hla/`、`work/<s>/kir/` | ~0.5 h |
| `08_mito.sh` | mutserve 2（同质 / 异质性）、GATK Mutect2 线粒体模式交叉验证、HaploGrep3 单倍群 | `work/<s>/mito/` | ~0.3 h |
| `09_ychr.sh` | 男性：chrY 小 BAM → Yleaf（YFull / FTDNA / ISOGG 三棵树） | `work/<s>/ychr/` | ~0.3 h |
| `04_cohort_genotype.sh` | GLnexus 联合分型（工作目录里所有人）→ 拆多等位、左对齐 | `work/joint/cohort.dv.norm.vcf.gz` | ~0.5 h |
| `10_annotate.sh` | 只注释外显子 ±50 bp + ClinVar 位点：VEP 116（merged 缓存、gnomAD v4.1、MANE）+ AlphaMissense + REVEL + NMD + ClinVar | `work/annot/cohort.vep.vcf.gz` | ~1 h |
| `10b_spliceai.sh` | 解读范围内约 5000 个基因上的罕见变异跑 SpliceAI，写回 INFO/SpliceAI | `work/annot/cohort.vep.spliceai.vcf.gz` | GPU ~1 h；CPU 数小时 |
| `11_findings.sh` | 分级（ACMG SF / 携带 / 双等位 / 神经肌肉 / 风险等位）、纯合功能缺失基因、性状位点、质控汇总 | `work/findings/` | 数分钟 |
| `11b_evidence.sh` | 关键变异的读段证据（igv-reports 单文件 HTML） | `results/<s>/evidence_igv.html` | 数分钟 |
| `12_pgx.sh` | PharmCAT 3.4（gVCF 展开参考位点；CYP2D6 用 Cyrius、HLA-A/B 用 T1K 作外部分型） | `work/<s>/pgx/` | 数分钟 |
| `13_ancestry.sh` | somalier relate / ancestry；plink2 PCA（全球 + POP 内部，1000 Genomes LD 修剪 SNP），样本投影 | `work/ancestry/` | ~1 h（首次建面板） |
| `14_prs.sh` | PGS Catalog 评分（默认东亚 16 个）：参考人群打分 + 样本从 gVCF 取基因型打分 → 百分位 | `work/prs/summary.json` | 1–2 h |
| `15_sv.sh`（可选） | Manta、Delly、CNVpytor → AnnotSV → sv_findings.py 严格过滤 | `work/<s>/sv/`、`work/findings/<s>.sv.tsv` | 3–5 h |
| `report` | make_figures.py（图）+ collect_results.py（结果底稿） | `results/<s>/figs/`、`summary.md/json` | 数分钟 |

GPU / CPU 实测（24 核 EPYC 7F72 + RTX 4090，30× 原始数据）：GPU 比对 16 分钟（fq2bam 全程 27.5 分钟）、DeepVariant 27 分钟；
CPU 24 线程：bwa 0.7.19 约 3.9 小时、bwa-mem2 约 2.0 小时、DeepVariant CPU 约 5.3 小时。

## 扩展分析（scripts/ext/，`run_all.sh --with-ext`）

参数在 `assets/population/<POP>.json`。

| 步骤 | 做什么 | 主要输出 |
|---|---|---|
| `x00_panel.sh` | 参考人群 MAF≥1% 的常染色体 SNP（东亚约 790 万）；1000 Genomes 已定相基因型按染色体导出；样本在这些位点的基因型 | `work/ext/panel/` |
| `x01_lai.sh` | Beagle 5.5 定相 → RFMix v2 局部祖源；留出的参考个体 + 对照人群同流程 | `work/ext/lai/` |
| `x02_aadr.sh` | AADR 1240K 位点转 GRCh38、样本编码；挑地区古今个体 + 外群；smartpca 投影；outgroup f3；qpAdm（含个体级对照） | `work/ext/aadr/` |
| `x03_roh.sh` | bcftools roh（参考人群频率 + 遗传图谱），样本与参考人群同参数 | `work/ext/roh/` |
| `x04_archaic.sh` | hmmix 二倍体流程 + 4 个古人类基因组归类；参考人群对照 | `work/ext/archaic/` |
| `x05_aging.sh` | TelSeq 端粒、线粒体拷贝数、Y 丢失、克隆性造血筛查 | `work/ext/aging/` |
| `x06_bloodgroup_cnv_stats.sh` | 红细胞血型扩展分型、拷贝数会变的基因、基因组概览 | `work/ext/` |
| `x07_virome.sh` | 未比对 + EBV 诱饵读段 → Kraken2 | `work/ext/virome/` |

## 磁盘

单个样本：FASTQ 30–50 GB（gz）、BAM 60–80 GB、gVCF ~5 GB、其他 ~20 GB；参考与数据库 ~70 GB（核心）+ 20（SV）+ 60（扩展，含面板 VCF）；容器镜像 ~25 GB；临时文件峰值 ~100 GB。
跑完可删：`tmp/`、`work/ext/panel/chr/`、`work/ext/lai/chr/`（只留汇总）、`work/<s>/str/*.bam`。建议备份：BAM（或 CRAM）、gVCF、`results/`。

## 多人

每人一份配置（`config.<名字>.env`），运行时加 `WGS_CONFIG=`。单样本步骤（01–03、05–09、11b、12、15、x01、x05–x07）按配置里的样本跑；
多人步骤（04、10、10b、11、13、14、x00、x02–x04）自动包含工作目录里所有已有 gVCF 的人，有新人加入会自动重算（`run_all.sh` 记录了样本集合）。
两个人时 `11_findings.sh` 额外输出 `pair.*.tsv`（同一隐性基因都携带）。
