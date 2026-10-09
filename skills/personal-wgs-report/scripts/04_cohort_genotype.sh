#!/usr/bin/env bash
# DeepVariant gVCF → GLnexus（DeepVariantWGS 配置）分型 → 拆多等位、左对齐的规范化 VCF：work/joint/cohort.dv.norm.vcf.gz
# 默认只有本样本；工作目录里如果还有别的样本（work/<名字>/<名字>.dv.g.vcf.gz，比如家人），会一起联合分型。
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
J=work/joint; mkdir -p $J; rm -rf $J/glnexus.DB
gvcfs=$(ls work/*/*.dv.g.vcf.gz)
log "联合分型：$(echo $gvcfs | wc -w) 个样本"
$D $IMG_GLNEXUS sh -c "glnexus_cli --config DeepVariantWGS --dir $J/glnexus.DB --threads $THREADS --mem-gbytes $MEM_GB $gvcfs > $J/cohort.dv.glnexus.bcf"
rm -rf $J/glnexus.DB
$D $IMG_BCFTOOLS sh -c "bcftools norm -m -any -f $REF $J/cohort.dv.glnexus.bcf -Oz -o $J/cohort.dv.norm.vcf.gz --threads 8 \
  && bcftools index -t $J/cohort.dv.norm.vcf.gz && bcftools stats -s - $J/cohort.dv.norm.vcf.gz > $J/cohort.dv.norm.stats.txt"
log "JOINT DONE"
