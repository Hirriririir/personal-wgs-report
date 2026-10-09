# 所有 shell 脚本的公共入口（用 source 引入，不要直接执行）：
#   进入工作目录、读取配置、给出默认值、定义容器镜像和小工具。
# 约定：每个脚本都在工作目录里运行，工作目录下有 scripts/（指向本仓库的软链接）、raw/ ref/ work/ qc/ logs/ results/。
# 同一工作目录分析第二个人（比如伴侣、家人）：给他写一份 config.<名字>.env，运行时加 WGS_CONFIG=config.<名字>.env；
# 参考、数据库、工具共用，每人的结果在 work/<样本名>/，联合分型等多人步骤会自动把工作目录里所有人一起算。
set -euo pipefail

WGS_WORKDIR=${WGS_WORKDIR:-$PWD}
cd "$WGS_WORKDIR"
export WGS_CONFIG=${WGS_CONFIG:-config.env}
if [ ! -f "$WGS_CONFIG" ]; then
  echo "找不到 $WGS_WORKDIR/$WGS_CONFIG：先运行 scripts/setup/init_workdir.sh <工作目录>，或在工作目录里执行本脚本" >&2
  exit 1
fi
# shellcheck disable=SC1090
. "./$WGS_CONFIG"

: "${SAMPLE:?$WGS_CONFIG 里要写 SAMPLE}"
SEX=${SEX:-auto}
THREADS=${THREADS:-$(nproc)}
MEM_GB=${MEM_GB:-32}
USE_GPU=${USE_GPU:-auto}
GPU_MODE=${GPU_MODE:-auto}
POP=${POP:-EAS}
LIBRARY=${LIBRARY:-lib1}
PLATFORM_MODEL=${PLATFORM_MODEL:-NovaSeqX}
OPTICAL_DUP_PIXEL_DISTANCE=${OPTICAL_DUP_PIXEL_DISTANCE:-2500}

# 容器镜像（版本固定，结果可复现；国内可把 BIOC 换成镜像站前缀）
BIOC=${BIOC:-quay.io/biocontainers}
IMG_SAMTOOLS=$BIOC/samtools:1.24--h9dcdb79_1
IMG_BCFTOOLS=$BIOC/bcftools:1.24--h118bc1c_4
IMG_BWA=$BIOC/bwa:0.7.19--h577a1d6_1
IMG_MOSDEPTH=$BIOC/mosdepth:0.3.14--h87be163_2
IMG_VERIFYBAMID=$BIOC/verifybamid2:2.0.3--hc004090_0
IMG_SOMALIER=$BIOC/somalier:0.3.5--h5205c93_0
IMG_GLNEXUS=$BIOC/glnexus:1.4.1--h17e8430_5
IMG_EH=$BIOC/expansionhunter:5.0.0--hc26b3af_5
IMG_REVIEWER=$BIOC/reviewer:0.2.7--h48da230_0
IMG_STRANGER=$BIOC/stranger:0.10.2--pyhdfd78af_0
IMG_SMN=$BIOC/smncopynumbercaller:1.1.2--py312h7e72e81_1
IMG_T1K=$BIOC/t1k:1.0.10--h5814d7d_0
IMG_MUTSERVE=$BIOC/mutserve:2.0.3--hdfd78af_0
IMG_HAPLOGREP=$BIOC/haplogrep3:3.2.2--hdfd78af_1
IMG_PLINK2=$BIOC/plink2:2.0.0a.7.8--hea6131b_0
IMG_KRAKEN2=$BIOC/kraken2:2.17.2--pl5321h3be2455_0
IMG_TELSEQ=$BIOC/telseq:0.0.2--h06902ac_8
IMG_BEAGLE=$BIOC/beagle:5.5_27Feb25.75f--hdfd78af_0
IMG_RFMIX=$BIOC/rfmix:2.03.r0.9505bfa--h503566f_8
IMG_ADMIXTOOLS=$BIOC/admixtools:8.0.2--h75d7a4a_0
IMG_EIGENSOFT=$BIOC/eigensoft:9.0.0--h519d32d_0
IMG_MANTA=$BIOC/manta:1.6.0--py27h9948957_6
IMG_DELLY=$BIOC/delly:2.7.0--h3752d28_0
IMG_CNVPYTOR=$BIOC/cnvpytor:1.3.2--pyhdfd78af_0
IMG_SURVIVOR=$BIOC/survivor:1.0.7--h077b44d_7
IMG_ANNOTSV=$BIOC/annotsv:3.5.10--hdfd78af_0
IMG_GATK=${IMG_GATK:-broadinstitute/gatk:4.7.0.0}
IMG_VEP=${IMG_VEP:-ensemblorg/ensembl-vep:release_116.2}
IMG_PHARMCAT=${IMG_PHARMCAT:-pgkb/pharmcat:3.4.0}
IMG_DEEPVARIANT=${IMG_DEEPVARIANT:-google/deepvariant:1.9.0}
IMG_PARABRICKS=${IMG_PARABRICKS:-nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1}

REF=ref/GRCh38/GRCh38.fa
BIN=tools/bin                      # fastp、seqkit、bwa-mem2 等静态程序（setup/install_tools.sh 下载）
D=scripts/dr.sh                     # CPU 容器入口
PY=${PY:-tools/venv-annot/bin/python}   # 注释 / 解读用的 Python
PYB=${PYB:-tools/venv-ext/bin/python}   # 扩展分析和画图用的 Python

[ -n "${HTTPS_PROXY:-}" ] && export https_proxy=$HTTPS_PROXY http_proxy=$HTTPS_PROXY HTTPS_PROXY HTTP_PROXY=$HTTPS_PROXY

log() { echo "[$(date "+%F %T")] $*"; }

# 是否用 GPU：USE_GPU=yes/no，auto 时看有没有 nvidia-smi
use_gpu() {
  case "$USE_GPU" in
    yes) return 0 ;;
    no) return 1 ;;
    *) command -v nvidia-smi > /dev/null 2>&1 && nvidia-smi -L > /dev/null 2>&1 ;;
  esac
}

# 性别：SEX=auto 时用 03_bam_qc 推断的结果（ambiguous 按 female 处理，对重复扩增和 X 染色体更保守）
sample_sex() {
  local s=${SEX}
  if [ "$s" = auto ]; then
    s=$(cut -f1 "work/$SAMPLE/qc/$SAMPLE.sex.txt" 2>/dev/null || echo ambiguous)
    [ "$s" = male ] || s=female
  fi
  echo "$s"
}

# 工作目录里已有变异检测结果的所有样本（多人步骤用）
all_samples() { for g in work/*/*.dv.g.vcf.gz; do [ -e "$g" ] && basename "$g" .dv.g.vcf.gz; done; }

# FASTQ 列表（逗号分隔 → 数组）
fastq_lists() {
  IFS=',' read -r -a FQ1 <<< "${FASTQ_R1:?$WGS_CONFIG 里要写 FASTQ_R1}"
  IFS=',' read -r -a FQ2 <<< "${FASTQ_R2:?$WGS_CONFIG 里要写 FASTQ_R2}"
  [ "${#FQ1[@]}" = "${#FQ2[@]}" ] || { echo "FASTQ_R1 与 FASTQ_R2 的文件数不一样" >&2; exit 1; }
}
