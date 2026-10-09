#!/usr/bin/env bash
# 纯合片段（ROH）：父母有没有血缘关系、祖上人群大小。样本与 1000 Genomes 参考人群无亲缘个体用同一套位点、同一参数：
#   bcftools roh --GTs-only 30，等位基因频率 = 参考人群，遗传图谱 = Beagle GRCh38 图谱。先跑 ext/x00_panel.sh
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
O=work/ext/roh; PN=work/ext/panel; MAP=ref/maps/chr_in_chrom_field
mkdir -p $O/map $O/chr
PL="$D $IMG_PLINK2 plink2"; BCF="$D $IMG_BCFTOOLS bcftools"
JOBS=${JOBS:-6}
# 1) 参考人群等位基因频率 → bcftools AF 文件
if [ ! -s $O/pop.af.tab.gz.tbi ]; then
  $PL --pfile ref/pca/all_hg38 vzs --keep $PN/pop_unrel.ids --chr 1-22 --snps-only just-acgt --max-alleles 2 \
      --set-all-var-ids '@:#:$r:$a' --rm-dup exclude-all --extract $PN/sites.ids --freq --out $O/pop --threads $THREADS > $O/pop.freq.log 2>&1
  awk 'NR>1{split($2,a,":"); print "chr"$1"\t"a[2]"\t"$3","$4"\t"$5}' $O/pop.afreq | sort -k1,1V -k2,2n > $O/pop.af.tab
  $D $IMG_SAMTOOLS sh -c "bgzip -f $O/pop.af.tab && tabix -f -s1 -b2 -e2 $O/pop.af.tab.gz"
fi
# 2) 遗传图谱 → bcftools 格式（位置  速率 cM/Mb  累计 cM）
for c in $(seq 1 22); do
  [ -s $O/map/chr$c.txt ] || awk 'BEGIN{print "position COMBINED_rate(cM/Mb) Genetic_Map(cM)"} {r=(NR>1 && $4>p)?($3-g)/($4-p)*1e6:0; print $4, r, $3; p=$4; g=$3}' \
    $MAP/plink.chrchr$c.GRCh38.map > $O/map/chr$c.txt
done
# 3) 样本
for s in $(all_samples); do
  $BCF roh --GTs-only 30 --AF-file $O/pop.af.tab.gz -m "$O/map/{CHROM}.txt" -O r -o $O/$s.roh.txt $PN/$s.sites.vcf.gz > $O/$s.roh.log 2>&1
  grep "^RG" $O/$s.roh.txt > $O/$s.roh.RG.txt || true
done
# 4) 参考人群（按染色体并行）
one() {
  c=$1
  [ -s $O/chr/pop.chr$c.roh.txt ] && return 0
  $BCF roh --GTs-only 30 --AF-file $O/pop.af.tab.gz -m "$O/map/{CHROM}.txt" -S $PN/pop_unrel.ids -O r -o $O/chr/pop.chr$c.roh.txt \
    $PN/chr/kg.chr$c.vcf.gz > $O/chr/pop.chr$c.log 2>&1
}
export -f one; export O PN BCF
seq 1 22 | xargs -P $JOBS -I{} bash -c 'one {}'
cat $O/chr/pop.chr*.roh.txt | grep "^RG" > $O/pop.roh.RG.txt
$PYB scripts/ext/roh_summary.py > $O/summary.log 2>&1 && cat $O/summary.log
log "ROH DONE"
