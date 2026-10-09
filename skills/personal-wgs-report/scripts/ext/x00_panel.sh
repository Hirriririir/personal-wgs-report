#!/usr/bin/env bash
# 扩展分析的公共参考面板（只做一次，约 1 小时、15 GB）：
#   work/ext/panel/sites.tsv      参考人群（config 的 POP）无亲缘个体中 MAF≥1% 的常染色体双等位 SNP（东亚约 790 万个）
#   work/ext/panel/chr/kg.chr*.vcf.gz   1000 Genomes 3202 人在这些位点上的已定相基因型（定相参考、局部祖源、ROH 对照）
#   work/ext/panel/<样本>.sites.vcf.gz   样本在这些位点上的基因型（从 gVCF 取，参考区块 GQ≥20 记 0/0）
#   work/ext/panel/gmap.tsv       RFMix 用的遗传图谱
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
O=work/ext/panel; P=ref/pca; MAP=ref/maps/chr_in_chrom_field
mkdir -p $O/chr
PL="$D $IMG_PLINK2 plink2"; BCF="$D $IMG_BCFTOOLS bcftools"
JOBS=${JOBS:-4}; T=$(( THREADS / JOBS > 1 ? THREADS / JOBS : 1 ))
$PYB - <<'PY'
import os, sys
sys.path.insert(0, "scripts/ext")
import popcfg, wgs
kg = popcfg.kg_samples()
ids = sorted(i for i, (sp, p) in kg.items() if sp == wgs.POP)
open("work/ext/panel/pop_unrel.ids", "w").write("".join(f"{i}\n" for i in ids))
print(f"{wgs.POP} 无亲缘个体 {len(ids)} 人")
PY
if [ ! -s $O/sites.tsv ]; then
  $PL --pfile $P/all_hg38 vzs --keep $O/pop_unrel.ids --chr 1-22 --snps-only just-acgt --max-alleles 2 --maf 0.01 \
      --set-all-var-ids '@:#:$r:$a' --rm-dup exclude-all --make-just-pvar cols= --out $O/sites --threads $THREADS > $O/sites.log 2>&1
  awk '!/^#/{print $3; print "chr"$3 > "'$O/sites.chrids'"}' $O/sites.pvar > $O/sites.ids
  awk '!/^#/{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5}' $O/sites.pvar > $O/sites.tsv
fi
log "参考面板位点：$(wc -l < $O/sites.tsv)"
[ -s $O/gmap.tsv ] || for c in $(seq 1 22); do awk '{print $1"\t"$4"\t"$3}' $MAP/plink.chrchr$c.GRCh38.map; done > $O/gmap.tsv
export O P PL BCF T
one() {
  set -euo pipefail
  c=$1; d=$O/chr
  [ -s $d/kg.chr$c.vcf.gz.csi ] && return 0
  $PL --pfile $P/all_hg38 vzs --chr $c --snps-only just-acgt --max-alleles 2 --set-all-var-ids '@:#:$r:$a' --rm-dup exclude-all \
      --extract $O/sites.chrids --export vcf bgz id-paste=iid --output-chr chrM --out $d/kg.chr$c --threads $T > $d/kg.chr$c.log 2>&1
  $BCF index -f $d/kg.chr$c.vcf.gz
}
export -f one
seq 1 22 | xargs -P $JOBS -I{} bash -c 'one {}'
for s in $(all_samples); do
  [ -s $O/$s.sites.vcf.gz.csi ] && continue
  $PY scripts/interpret/gvcf_genotype_stream_par.py $O/sites.tsv $O/$s.sites.vcf.gz $s=work/$s/$s.dv.g.vcf.gz > $O/$s.sites.log
  $BCF index -f $O/$s.sites.vcf.gz
done
log "PANEL DONE"
