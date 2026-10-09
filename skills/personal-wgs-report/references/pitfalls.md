# 常见问题与坑（都是实际跑过遇到的）

## 环境

- **Docker 权限**：当前用户要在 `docker` 组里（`sudo usermod -aG docker $USER` 后重新登录）。容器一律以当前用户 uid 运行（`-u $(id -u):$(id -g)`），输出文件不会变成 root 的。
- **路径挂载**：容器按原路径挂载工作目录；工作目录顶层的软链接（`scripts`、共享的 `ref` 等）和 config 里 FASTQ 所在目录会自动挂上（`scripts/lib/mounts.sh`）。
  FASTQ 放在其他位置时，确保在 config 里写的是绝对路径。
- **GPU**：Parabricks 需要 NVIDIA 驱动 ≥ 525、显存 ≥ 16 GB（24 GB 更稳）。装了 nvidia-container-toolkit 用 `GPU_MODE=toolkit`；没装时 `GPU_MODE=manual` 会手动挂 `/dev/nvidia*` 和驱动库（`scripts/gpu_docker.sh`），
  驱动库不在 `/usr/lib/x86_64-linux-gnu` 的发行版要改脚本里的 `L=`。显存不足时 fq2bam 会报 out of memory，可加 `--low-memory`（改 02 脚本）。
- **国内网络**：quay.io / Docker Hub 慢 → `BIOC=` 换镜像站，或在别处 `docker pull` + `docker save/load`。NCBI / Ensembl / GitHub / Google Cloud 慢 → `HTTPS_PROXY=`。
  `fetch.py` 支持断点续传和 MD5 核对，中断后重跑即可。
- **REVEL / PanelApp 下载 403**：部分地区被拦，按脚本提示手动下载；不放也能跑。
- **bwa-mem2 建索引** 要 ~80 GB 内存；内存不够就设 `USE_GPU=no` 也改用经典 bwa（把 02 脚本里的 bwa-mem2 换成 `$D $IMG_BWA bwa mem`，速度慢约一倍）。

## 比对与变异检测

- **Parabricks 加 `--gvcf`** 会同时写出 `<s>.dv.vcf.gz`（只含变异）和 `<s>.dv.g.vcf.gz`；后面多个步骤两个都要用。
- **DeepVariant CPU 很慢**（24 线程约 5 小时）：`--num_shards` = THREADS；中间文件在 `tmp/<s>/dv`，空间要几十 GB。
- **光学重复多**（NovaSeq X 常见 15–25%）：属于平台特性，不是样本问题；`OPTICAL_DUP_PIXEL_DISTANCE` 图案化芯片用 2500。
- **多对 FASTQ**：同一文库多条泳道就是多对，按顺序全写进 FASTQ_R1 / FASTQ_R2（逗号分隔）；不要先 cat 成一个文件（会丢失读段组信息，但其实不影响结果——只是标重复按文库进行）。

## 质控

- **VerifyBamID2 失败**：常见于参考资源下载不全或 BAM 太大内存不够；不影响后续，报告写“污染未评估”，或用 somalier 的杂合比例粗看。
- **性别 ambiguous**：X/Y 深度比不在典型范围（如 XXY、嵌合、低深度）；`sample_sex()` 会按 female 处理（对 STR 和 X 连锁更保守），并在报告里说明。

## 注释与解读

- **VEP 很吃内存**：`--fork` 上限 16、`--buffer_size 20000`；内存 <64 GB 时把 fork 调小。
- **SpliceAI**：CPU 也能跑，但很慢；失败不影响其他步骤（解读少一类剪接预测）。TensorFlow 2.15 + numpy<2 已固定版本，别升级。
- **RFC1 被标成扩增**：ExpansionHunter 的简并单元把良性 AAAAG / AAAGG 扩增也算进去，一定要看 `rfc1_motif.json`（见 interpretation.md）。
- **线粒体 8271–8289 缺失、310 / 16189 poly-C**：mutserve 会报成“异质性”，多为比对伪影或单倍群标志性缺失（如 B 类群的 9 bp 缺失），汇总时已排除。
- **PharmCAT 不认识的 CYP2D6 等位**（如 `*36+*10`）：日志里 “Undocumented CYP2D6 named allele in outside call”，PharmCAT 对 CYP2D6 不给推荐；按 Cyrius 结果手算活性分（interpretation.md）。
- **PharmCAT 报 RYR1 “Uncertain Susceptibility”**：意思是没查到已知易感变异，不是“有不确定风险”。
- **ClinVar 冲突解读**（conflicting）：不进 P/LP；在 all_candidates.tsv 里能看到，需要专业判断。
- **低 GQ / 低 VAF 的“致病”变异**：qc_note 非空的多为比对伪影（同源区、重复区）。打开 `evidence_igv.html` 看读段：支持读段都在一条链上、都在读段末端、比对质量低 → 多半是假的。

## PRS

- 评分文件只给效应等位、没给另一等位（部分 TPMI 评分）：`prs_weights.py` 自动按 1000 Genomes 位点补全 REF/ALT，并去掉回文位点。
- 评分在参考人群和样本上必须用**同一套位点**：样本的基因型从 gVCF 取（参考区块 GQ≥20 记 0/0），所以缺失不多；`union.genotype.log` 里能看到取到了多少。
- 换人群（POP≠EAS）必须换评分：用 `pgs_search.py` 找该人群开发 / 验证过的评分。

## 扩展分析

- **hmmix 需要 bcftools**：`x04_archaic.sh` 会在 `tools/bin-docker/` 里生成一个调用容器的 bcftools 包装并放进 PATH。
- **Beagle 内存**：每条染色体 ~24 GB（`BEAGLE_MEM`），并行数 `JOBS` 按内存调。
- **AADR 个体标签**：带 `_o`（离群）、`Ignore`、`_dup`、`_rel` 的已排除；qpAdm 的来源人群和外群名单是针对东亚写的，其他人群要重选（customizing.md）。
- **smartpca 投影**：古代样本位点少（低覆盖），投影会向原点收缩，看相对位置别看绝对坐标。
- **Kraken2 结果**：“未能归类”占大头是正常的；standard-8GB 库是精简版，罕见物种会漏。

## 其他

- **中文字体**：图里中文显示成方块 → `install_tools.sh` 会下载 Noto Sans SC 到 `tools/fonts/`，或系统装 `fonts-noto-cjk`。
- **pgrep / pkill 匹配到自己**：用 `pgrep -f "Expansio[n]Hunter"` 这类写法，避免匹配到监控命令本身。
- **重跑某一步**：删 `logs/done/<步骤>.<样本>`（多人步骤是 `.cohort`）后再 `run_all.sh`，或 `run_all.sh --only <步骤>`。
