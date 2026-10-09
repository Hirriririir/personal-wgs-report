#!/usr/bin/env bash
# 药物基因组：PharmCAT 3.4（CPIC / DPWG / FDA 指南）
#   输入：从 DeepVariant gVCF 展开 PharmCAT 位点（高质量参考区块 → 明确的 0/0，避免把“没测到”当成“野生型”）
#   外部分型：CYP2D6 用 Cyrius（结构变异 / 杂合基因），HLA-A / HLA-B 用 T1K
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
set +e
s=$SAMPLE; O=work/$s/pgx; mkdir -p $O ref/pgx
[ -s ref/pgx/pharmcat_positions.vcf.bgz.csi ] || for f in pharmcat_regions.bed pharmcat_positions.vcf.bgz pharmcat_positions.vcf.bgz.csi; do
  docker run --rm $IMG_PHARMCAT cat /pharmcat/$f > ref/pgx/$f; done
$PY scripts/interpret/pgx_input.py $s ref/pgx/pharmcat_positions.vcf.bgz $O/$s.pgx_input.vcf || { echo "pgx_input 失败" >&2; exit 1; }
$D $IMG_BCFTOOLS sh -c "bgzip -f $O/$s.pgx_input.vcf && bcftools index -f -t $O/$s.pgx_input.vcf.gz"
# outside calls（基因 TAB 双倍型）
$PY - "$s" > $O/$s.outside.tsv <<'PY'
import csv, sys
s = sys.argv[1]
try:   # CYP2D6（Cyrius）：有多个候选（用 ; 分隔）时不交给 PharmCAT
    for r in csv.DictReader(open(f"work/{s}/special/{s}.cyp2d6.tsv"), delimiter="\t"):
        g = r.get("Genotype", "None")
        if g and g != "None" and ";" not in g:
            print(f"CYP2D6\t{g}")
except FileNotFoundError:
    pass
def two_field(a):   # T1K 等位基因取两位字段；同一位置多个候选只取第一个
    a = a.split(",")[0].split("*", 1)[1]
    return "*" + ":".join(a.split(":")[:2])
try:
    for line in open(f"work/{s}/hla/{s}_genotype.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] in ("HLA-A", "HLA-B") and f[2] != ".":
            a1 = two_field(f[2]); a2 = two_field(f[5]) if len(f) > 5 and f[5] != "." else a1
            print(f"{f[0]}\t{a1}/{a2}")
except FileNotFoundError:
    pass
PY
# pharmcat_pipeline 会自动读取输入 VCF 同目录下的 <-bf 名>.outside.tsv
$D $IMG_PHARMCAT pharmcat_pipeline $O/$s.pgx_input.vcf.gz -o $O -bf $s \
  -reporterJson -reporterHtml -reporterCallsOnlyTsv -re -matcherHtml > $O/pharmcat.log 2>&1 || echo "PharmCAT 失败（见 $O/pharmcat.log）"
log "PGX DONE $s：$(ls $O/*.report.html 2>/dev/null)"
