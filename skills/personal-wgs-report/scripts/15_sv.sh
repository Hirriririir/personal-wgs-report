#!/usr/bin/env bash
# 结构变异 / CNV（可选，约 3–5 小时）：Manta（断点）+ Delly（断点 + 读对）+ CNVpytor（读深）→ AnnotSV 注释 → sv_findings.py 过滤
# 健康人的大片段致病 CNV 很少，SV 软件的假阳性很多，所以解读规则偏严（见 sv_findings.py 开头）。
# 用法：scripts/15_sv.sh [manta] [delly] [cnvpytor] [annotate]   默认全做
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
s=$SAMPLE; O=work/$s/sv; steps=${*:-manta delly cnvpytor annotate}
mkdir -p $O/annot ref/sv
if [ ! -s ref/sv/primary.bed.gz.tbi ]; then   # 只在 1–22、X、Y 上找（不含 decoy / 未定位片段 / chrM）
  awk -v OFS="\t" '$1~/^chr([0-9]+|X|Y)$/{print $1,0,$2}' $REF.fai > ref/sv/primary.bed
  $D $IMG_BCFTOOLS sh -c "bgzip -f ref/sv/primary.bed && tabix -f -p bed ref/sv/primary.bed.gz"
fi
for step in $steps; do case $step in
manta)
  rm -rf $O/manta
  $D $IMG_MANTA sh -c "configManta.py --bam work/$s/$s.bam --referenceFasta $REF --callRegions ref/sv/primary.bed.gz --runDir $O/manta \
    && $O/manta/runWorkflow.py -j $THREADS" > $O/manta.log 2>&1
  log "Manta 完成" ;;
delly)
  $D $IMG_DELLY delly sr -h $(( THREADS > 8 ? 8 : THREADS )) -g $REF -x ref/sv/delly_human.hg38.excl.tsv -o $O/delly.bcf work/$s/$s.bam > $O/delly.log 2>&1
  $D $IMG_BCFTOOLS sh -c "bcftools view $O/delly.bcf -Oz -o $O/delly.vcf.gz && bcftools index -f -t $O/delly.vcf.gz"
  log "Delly 完成" ;;
cnvpytor)
  P=$O/$s.pytor; C=tools/venv-cnv/bin/cnvpytor; rm -f $P
  ( $C -root $P -rd work/$s/$s.bam -chrom $(echo chr{1..22} chrX chrY) -j 8
    $C -root $P -gc $REF -j 8 || true
    $C -root $P -his 1000 10000 100000 -j 8
    $C -root $P -partition 1000 10000 100000 -j 8
    $C -root $P -call 1000 > $O/cnvpytor.calls.1k.tsv
    $C -root $P -call 10000 > $O/cnvpytor.calls.10k.tsv ) > $O/cnvpytor.log 2>&1
  log "CNVpytor 完成" ;;
annotate)
  A="$D $IMG_ANNOTSV AnnotSV"; B="$D $IMG_BCFTOOLS"
  # Manta：PASS、≥50 bp；Delly：PASS
  $B sh -c "bcftools view -f PASS -i 'SVLEN>=50 || SVLEN<=-50 || SVTYPE=\"BND\" || SVTYPE=\"INV\"' $O/manta/results/variants/diploidSV.vcf.gz -Oz -o $O/annot/manta.pass.vcf.gz \
    && bcftools index -f -t $O/annot/manta.pass.vcf.gz && bcftools view -f PASS $O/delly.vcf.gz -Oz -o $O/annot/delly.pass.vcf.gz && bcftools index -f -t $O/annot/delly.pass.vcf.gz"
  # CNVpytor：e-val1<1e-4、q0<0.5（MAPQ0 读段比例）、pN<0.5（N 比例）、≥5 kb
  awk -v s=$s -F'\t' 'BEGIN{OFS="\t"} {split($2,a,/[:-]/); if ($5<1e-4 && $9<0.5 && $10<0.5 && $3>=5000) print a[1],a[2],a[3],($1=="deletion"?"DEL":"DUP"),s,$4}' \
    $O/cnvpytor.calls.1k.tsv | sort -k1,1V -k2,2n > $O/annot/cnvpytor.bed
  for x in manta delly; do
    $A -SVinputFile $O/annot/$x.pass.vcf.gz -SVinputInfo 1 -annotationsDir ref/sv/annotsv -genomeBuild GRCh38 -SVminSize 50 \
       -outputDir $O/annot -outputFile $O/annot/$x.annotsv.tsv > $O/annot/$x.annotsv.log 2>&1 || log "AnnotSV 失败：$x"
  done
  $A -SVinputFile $O/annot/cnvpytor.bed -svtBEDcol 4 -samplesidBEDcol 5 -annotationsDir ref/sv/annotsv -genomeBuild GRCh38 \
     -outputDir $O/annot -outputFile $O/annot/cnvpytor.annotsv.tsv > $O/annot/cnvpytor.annotsv.log 2>&1 || log "AnnotSV 失败：cnvpytor"
  $PY scripts/interpret/sv_findings.py > logs/sv_findings_$s.log 2>&1 || log "sv_findings 失败（见 logs/sv_findings_$s.log）"
  log "SV 解读完成：work/findings/$s.sv.tsv" ;;
*) echo "未知步骤 $step"; exit 1;;
esac; done
log "SV DONE $s"
