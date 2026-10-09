#!/usr/bin/env bash
# 装宿主机上的小工具和 Python 环境（其余工具都在容器里）：
#   tools/bin：fastp、seqkit、bwa-mem2（静态二进制）
#   tools/venv-annot：SpliceAI（TensorFlow）、cyvcf2、pysam、pandas —— 注释与解读
#   tools/venv-ext：matplotlib、hmmix、pysam、cyvcf2、shapely —— 扩展分析与画图
#   tools/venv-illumina + Cyrius / Gauchian：CYP2D6、GBA1
#   tools/venv-yleaf + Yleaf：Y 染色体单倍群
#   tools/venv-cnv：CNVpytor（读深 CNV）
# 装了 uv 会快很多（https://docs.astral.sh/uv/），没有就用 python3 -m venv + pip。
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
mkdir -p tools/bin
cd tools
dl() { [ -s "$2" ] || curl -fSL --retry 5 -o "$2" "$1"; }
# --- 静态程序
if [ ! -x bin/fastp ]; then dl http://opengene.org/fastp/fastp.1.4.0 bin/fastp && chmod +x bin/fastp; fi
if [ ! -x bin/seqkit ]; then
  dl https://github.com/shenwei356/seqkit/releases/download/v2.14.0/seqkit_linux_amd64.tar.gz seqkit.tar.gz
  tar -xzf seqkit.tar.gz -C bin && rm seqkit.tar.gz
fi
if [ ! -x bin/bwa-mem2 ]; then
  dl https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 bwa-mem2.tar.bz2
  tar -xjf bwa-mem2.tar.bz2 && cp bwa-mem2-2.2.1_x64-linux/bwa-mem2* bin/ && rm -rf bwa-mem2.tar.bz2 bwa-mem2-2.2.1_x64-linux
fi
# --- Python 环境
mkvenv() {  # mkvenv <名字> <python 版本> <包...>
  local v=$1 py=$2; shift 2
  [ -x venv-$v/bin/python ] && { echo "venv-$v 已存在"; return; }
  if command -v uv > /dev/null; then
    uv venv -p "$py" venv-$v && VIRTUAL_ENV=$PWD/venv-$v uv pip install "$@"
  else
    python3 -m venv venv-$v && venv-$v/bin/pip install -q --upgrade pip && venv-$v/bin/pip install -q "$@"
  fi
}
mkvenv annot 3.11 "spliceai==1.3.1" "tensorflow==2.15.1" "numpy<2" cyvcf2 pysam pandas openpyxl intervaltree "igv-reports==1.17.0"
mkvenv ext 3.12 "hmmix==0.9.2" matplotlib pandas numpy scipy pysam cyvcf2 shapely scikit-learn zstandard CrossMap markdown
mkvenv illumina 3.12 pysam scipy statsmodels pandas numpy
mkvenv yleaf 3.12 pandas numpy networkx graphviz
if [ ! -x venv-cnv/bin/cnvpytor ]; then   # CNVpytor（读深 CNV，15_sv.sh 用）+ 它自带的 GRCh38 GC / 可比对性资源
  mkvenv cnv 3.11 "cnvpytor==1.3.2" && venv-cnv/bin/cnvpytor -download
fi
[ -d Cyrius ] || git clone -q https://github.com/Illumina/Cyrius.git && git -C Cyrius checkout -q a2ee08a
[ -d Gauchian ] || git clone -q https://github.com/Illumina/Gauchian.git && git -C Gauchian checkout -q e69ceee
if [ ! -d Yleaf ]; then
  git clone -q https://github.com/genid/Yleaf.git && git -C Yleaf checkout -q v4.1.5
  venv-yleaf/bin/pip install -q ./Yleaf 2> /dev/null || VIRTUAL_ENV=$PWD/venv-yleaf uv pip install ./Yleaf
fi
# 中文字体（画图用）：系统里没有中文字体时下载思源黑体 / Noto Sans SC（SIL OFL 开源字体）
if ! fc-list :lang=zh 2>/dev/null | grep -q . && [ ! -s fonts/NotoSansSC-Regular.otf ]; then
  mkdir -p fonts
  for w in Regular Medium Bold; do
    dl https://github.com/notofonts/noto-cjk/raw/main/Sans/SubsetOTF/SC/NotoSansSC-$w.otf fonts/NotoSansSC-$w.otf || true
  done
fi
echo "TOOLS DONE"
