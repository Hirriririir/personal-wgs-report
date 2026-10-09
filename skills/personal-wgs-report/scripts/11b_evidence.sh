#!/usr/bin/env bash
# 关键变异的读段证据：igv-reports 生成单个 HTML（内嵌读段 + 参考序列 + igv.js，可离线打开），给人复核“是不是真的”。
# 变异来源：work/findings 里本样本的 ACMG SF / 携带 / 双等位 / 风险等位 / 神经肌肉（P/LP 和倾向致病）以及 pair.* 文件。
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
s=$SAMPLE; O=work/findings; R=results/$s; mkdir -p $R
n=$($PY - "$s" <<'PY'
import csv, glob, os, sys
s = sys.argv[1]
keys = set()
files = [f"work/findings/{s}.{k}.tsv" for k in ("acmg_sf", "carrier", "biallelic", "nmd", "risk_alleles")] + glob.glob("work/findings/pair.*.tsv")
for f in files:
    if not os.path.exists(f):
        continue
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r.get("sample") not in (s, None):
            continue
        if f.endswith("nmd.tsv") and not (r["tier"].startswith(("P/LP", "LP")) or "倾向" in r["tier"]):
            continue
        keys.add((r["chrom"], int(r["pos"])))
with open(f"work/findings/{s}.evidence_sites.tsv", "w") as o:
    for c, p in sorted(keys):
        o.write(f"{c}\t{p}\n")
print(len(keys))
PY
)
[ "$n" -gt 0 ] || { log "没有需要看读段证据的变异"; exit 0; }
$D $IMG_BCFTOOLS sh -c "bcftools view -s $s -T $O/$s.evidence_sites.tsv work/annot/cohort.vep.vcf.gz -Ou \
  | bcftools annotate -x INFO/CSQ -Oz -o $O/$s.evidence.vcf.gz && bcftools index -f -t $O/$s.evidence.vcf.gz"
tools/venv-annot/bin/create_report $O/$s.evidence.vcf.gz --fasta $REF --tracks work/$s/$s.bam --flanking 60 --standalone \
  --samples $s --sample-columns GT DP AD GQ --title "关键变异读段证据（$s）" --output $R/evidence_igv.html > logs/evidence_$s.log 2>&1
log "EVIDENCE DONE：$n 个位点 → $R/evidence_igv.html"
