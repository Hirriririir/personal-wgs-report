#!/usr/bin/env bash
# 比对质控：mosdepth 覆盖度、samtools stats / flagstat、VerifyBamID2 污染、somalier 指纹；再按 X / Y 深度推断性别。
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
s=$SAMPLE; B=work/$s/$s.bam; Q=work/$s/qc; mkdir -p $Q work/somalier
$D $IMG_MOSDEPTH mosdepth -t 4 -n --fast-mode --by 1000 -T 1,10,15,20,30 $Q/$s $B
$D $IMG_SAMTOOLS sh -c "samtools flagstat -@ 8 $B > $Q/$s.flagstat.txt && samtools stats -@ 8 -r $REF $B > $Q/$s.samtools_stats.txt"
$D $IMG_VERIFYBAMID verifybamid2 --SVDPrefix ref/qc/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat \
  --Reference $REF --BamFile $B --NumThread 8 --Output $Q/$s.verifybamid2 > $Q/$s.verifybamid2.log 2>&1 || log "VerifyBamID2 失败（见 $Q/$s.verifybamid2.log），不影响后续"
$D $IMG_SOMALIER somalier extract -d work/somalier --sites ref/qc/somalier/sites.hg38.vcf.gz -f $REF $B
$PY scripts/interpret/infer_sex.py $s
log "BAM QC DONE $s"
