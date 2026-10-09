#!/usr/bin/env bash
# 解读：小变异分级（ACMG SF 次要发现 / 隐性携带 / 神经肌肉病基因 / 风险等位）、纯合功能缺失基因、性状位点、数据质量汇总
# 输出 work/findings/（每人 <样本>.*.tsv；工作目录里恰好两个人时多出 pair.*.tsv：同一隐性基因两人都携带）
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
mkdir -p work/findings logs
VCF=work/annot/cohort.vep.spliceai.vcf.gz
[ -s $VCF.tbi ] || { VCF=work/annot/cohort.vep.vcf.gz; log "没有 SpliceAI 结果，只用 VEP 注释"; }
$PY scripts/interpret/findings.py $VCF > logs/findings.log 2>&1 || { echo "findings.py 失败，见 logs/findings.log" >&2; exit 1; }
for s in $(all_samples); do
  $PY scripts/interpret/knockouts.py $s > work/findings/$s.knockouts.txt 2> logs/knockouts_$s.log || log "knockouts $s 失败"
done
$PY scripts/interpret/trait_snps.py > logs/trait_snps.log 2>&1 || log "trait_snps 失败（见 logs/trait_snps.log）"
$PY scripts/interpret/qc_summary.py > logs/qc_summary.log 2>&1 || log "qc_summary 失败（见 logs/qc_summary.log）"
log "FINDINGS DONE：$(tail -n +1 work/findings/summary.json | tr -d '\n ' | cut -c1-200)"
