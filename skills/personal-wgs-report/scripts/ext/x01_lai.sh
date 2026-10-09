#!/usr/bin/env bash
# 局部祖源（染色体“涂色”）：把样本每段染色体判成参考组之一（东亚默认：北方 = CHB，南方 = CDX + KHV）
#   定相：Beagle 5.5，参考 = 1000 Genomes 3202 人已定相单倍型；局部祖源：RFMix v2
#   参考组各留若干人不进参考，和其他对照人群一起跑同一流程，检验结果可信度
#   参数在 assets/population/<POP>.json 的 lai 段；先跑 ext/x00_panel.sh
# 用法：scripts/ext/x01_lai.sh [prep|phase|rfmix|summary|all]
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
s=$SAMPLE; PN=work/ext/panel; O=work/ext/lai; MAP=ref/maps/chr_in_chrom_field
mkdir -p $O/chr
BEAGLE="$D $IMG_BEAGLE beagle"; RFMIX="$D $IMG_RFMIX rfmix"; BCF="$D $IMG_BCFTOOLS bcftools"
JOBS=${JOBS:-3}; T=$(( THREADS / JOBS > 1 ? THREADS / JOBS : 1 )); STEP=${1:-all}
G=$($PYB -c 'import sys; sys.path.insert(0,"scripts/ext"); import popcfg; print(popcfg.load()["lai"].get("generations", 20))')
prep() {  # 参考组 / 对照组（固定随机种子，可复现）
  $PYB - <<'PY'
import random, sys
sys.path.insert(0, "scripts/ext")
import popcfg
cfg = popcfg.load()["lai"]
kg = popcfg.kg_samples()
by_pop = {}
for i, (sp, p) in kg.items():
    by_pop.setdefault(p, []).append(i)
random.seed(20261009)
q, ref = [], []
ref_pops = {p: g for g, ps in cfg["ref"].items() for p in ps}
for p in sorted(set(cfg["controls"]) | set(ref_pops)):
    ids = sorted(by_pop.get(p, []))
    random.shuffle(ids)
    k = cfg.get("holdout", 10)
    if p in cfg["controls"]:
        q += [(i, p) for i in ids[:k]]
    if p in ref_pops:
        ref += [(i, ref_pops[p]) for i in (ids[k:] if p in cfg["controls"] else ids)]
O = "work/ext/lai"
open(f"{O}/query_1kg.tsv", "w").write("".join(f"{i}\t{p}\n" for i, p in q))
open(f"{O}/query_1kg.ids", "w").write("".join(f"{i}\n" for i, _ in q))
open(f"{O}/sample_map.tsv", "w").write("".join(f"{i}\t{g}\n" for i, g in ref))
open(f"{O}/ref.ids", "w").write("".join(f"{i}\n" for i, _ in ref))
print(len(q), "个对照，", len(ref), "个参考个体，参考组", sorted(cfg["ref"]))
PY
}
phase_one() {
  set -euo pipefail
  c=$1; d=$O/chr
  if [ ! -s $d/$s.chr$c.phased.vcf.gz.csi ]; then
    $BEAGLE -Xmx${BEAGLE_MEM:-24}g gt=$PN/$s.sites.vcf.gz ref=$PN/chr/kg.chr$c.vcf.gz map=$MAP/plink.chrchr$c.GRCh38.map chrom=chr$c \
        out=$d/$s.chr$c.phased impute=false nthreads=$T seed=20261009 > $d/beagle.$s.chr$c.stdout 2>&1
    $BCF index -f $d/$s.chr$c.phased.vcf.gz
  fi
  if [ ! -s $d/ref.chr$c.vcf.gz.csi ]; then
    $BCF view -S $O/ref.ids -Oz -o $d/ref.chr$c.vcf.gz $PN/chr/kg.chr$c.vcf.gz && $BCF index -f $d/ref.chr$c.vcf.gz
    $BCF view -S $O/query_1kg.ids -Oz -o $d/q1kg.chr$c.vcf.gz $PN/chr/kg.chr$c.vcf.gz && $BCF index -f $d/q1kg.chr$c.vcf.gz
  fi
  if [ ! -s $d/query.$s.chr$c.vcf.gz.csi ]; then
    $BCF merge -m none -Oz -o $d/query.$s.chr$c.vcf.gz $d/$s.chr$c.phased.vcf.gz $d/q1kg.chr$c.vcf.gz && $BCF index -f $d/query.$s.chr$c.vcf.gz
  fi
}
rfmix_one() {
  set -euo pipefail
  c=$1; d=$O/chr
  [ -s $d/rfmix.$s.chr$c.msp.tsv ] && return 0
  $RFMIX -f $d/query.$s.chr$c.vcf.gz -r $d/ref.chr$c.vcf.gz -m $O/sample_map.tsv -g $PN/gmap.tsv -o $d/rfmix.$s.chr$c \
      --chromosome=chr$c -G $G -e 0 --n-threads=$T --random-seed=20261009 > $d/rfmix.$s.chr$c.log 2>&1
}
export -f phase_one rfmix_one
export O PN MAP BEAGLE RFMIX BCF T G s
case $STEP in prep|all) prep ;; esac
case $STEP in phase|all) seq 1 22 | xargs -P $JOBS -I{} bash -c 'phase_one {}' ;; esac
case $STEP in rfmix|all) seq 1 22 | xargs -P $JOBS -I{} bash -c 'rfmix_one {}' ;; esac
case $STEP in summary|all) $PYB scripts/ext/lai_summary.py $s > $O/summary.$s.log 2>&1 && cat $O/summary.$s.log ;; esac
log "LAI $STEP DONE"
