#!/usr/bin/env bash
# 扩展分析（祖源、古人类、病毒）用的资源，约 45 GB：
#   1000 Genomes 高深度 3,202 人（PLINK 2 格式，已定相）       → 主成分、局部祖源、纯合片段对照、PRS 百分位
#   AADR v66 古 DNA 数据集（1240K）                            → 主成分投影、f3、qpAdm
#   hmmix 辅助文件 + 4 个高覆盖古人类基因组                     → 尼安德特 / 丹尼索瓦片段
#   Beagle GRCh38 遗传图谱、UCSC hg19→hg38 坐标转换链           → 定相、AADR 坐标换算
#   UCSC cytoBand、Natural Earth 底图                          → 画染色体图、古 DNA 地图
#   Kraken2 standard-8GB 数据库                                → 比不上人类的读段分类（找病毒）
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
F="python3 scripts/setup/fetch.py"
mkdir -p ref/pca ref/aadr ref/hmmix/archaic ref/maps ref/liftover ref/kraken2 ref/geo ref/annot
# --- 1000 Genomes（PLINK 2 官方资源页的 3,202 人版本）
[ -s ref/pca/all_hg38.pgen ] || {
  $F url "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1" ref/pca/all_hg38.pgen.zst
  $D $IMG_PLINK2 plink2 --zst-decompress ref/pca/all_hg38.pgen.zst ref/pca/all_hg38.pgen && rm -f ref/pca/all_hg38.pgen.zst*; }
$F url "https://www.dropbox.com/s/ngbo2xm5ojw9koy/all_hg38_noannot.pvar.zst?dl=1" ref/pca/all_hg38.pvar.zst
$F url "https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1" ref/pca/all_hg38.psam
$F url "https://www.dropbox.com/s/129gx0gl2v7ndg6/deg2_hg38.king.cutoff.out.id?dl=1" ref/pca/deg2_hg38.king.cutoff.out.id
$F url https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt ref/pca/high_ld_hg38.txt
# --- AADR v66.p1（Harvard Dataverse）
$F dataverse doi:10.7910/DVN/FFIDCW ref/aadr 'v66.p1_1240K.aadr*' 'v66.p1__files.md5sum'
# --- hmmix
$F zenodo 11212339 ref/hmmix hg38_Outgroup_1000g_HGDP.txt hg38_mutationrate.bed hg38_strick_callability_mask.bed hg38_ancestral.tar.gz
[ -d ref/hmmix/hg38_ancestral ] || tar -xzf ref/hmmix/hg38_ancestral.tar.gz -C ref/hmmix
$F zenodo 13368126 ref/hmmix/archaic 'individuals_highcov.*'
# --- 遗传图谱、坐标转换
$F url https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip ref/maps/plink.GRCh38.map.zip
if [ ! -d ref/maps/chr_in_chrom_field ]; then
  python3 - <<'PY'
import os, zipfile
z = zipfile.ZipFile("ref/maps/plink.GRCh38.map.zip")
for d in ("chr_in_chrom_field", "no_chr_in_chrom_field"):
    os.makedirs(f"ref/maps/{d}", exist_ok=True)
for n in z.namelist():
    if not n.endswith(".map"):
        continue
    lines = z.read(n).decode().splitlines()
    base = os.path.basename(n)
    open(f"ref/maps/no_chr_in_chrom_field/{base}", "w").write("\n".join(lines) + "\n")
    open(f"ref/maps/chr_in_chrom_field/{base}", "w").write("\n".join("chr" + l if not l.startswith("chr") else l for l in lines) + "\n")
PY
fi
$F url https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz ref/liftover/hg19ToHg38.over.chain.gz
# --- 画图：染色体长度与着丝粒（UCSC cytoBand）、地图底图（Natural Earth 1:50m 陆地 / 湖泊 / 河流，公有领域）
$F url https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz ref/annot/cytoBand.hg38.txt.gz
for x in land lakes rivers_lake_centerlines; do
  $F url https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_$x.geojson ref/geo/ne_50m_$x.geojson
done
# --- Kraken2 standard-8GB（取 S3 上最新的一版）
if [ ! -s ref/kraken2/std8/hash.k2d ]; then
  k=$(curl -fsSL "https://genome-idx.s3.amazonaws.com/?list-type=2&prefix=kraken/k2_standard_08_GB_" | grep -oE "kraken/k2_standard_08_GB_[0-9]{8}\.tar\.gz" | sort | tail -1)
  $F url "https://genome-idx.s3.amazonaws.com/$k" "ref/kraken2/$(basename "$k")"
  mkdir -p ref/kraken2/std8 && tar -xzf "ref/kraken2/$(basename "$k")" -C ref/kraken2/std8
fi
echo EXT RESOURCES DONE
