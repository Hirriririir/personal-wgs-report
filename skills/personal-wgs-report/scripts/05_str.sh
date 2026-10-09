#!/usr/bin/env bash
# 重复扩增：ExpansionHunter v5 + stranger 目录（51 个致病位点，含 DMPK/CNBP/PABPN1/NOTCH2NLC/LRP12/GIPC1/RILPL1/RFC1/FGF14 等）
#           → stranger 按正常 / 前突变 / 致病阈值标注 → REViewer 逐位点画读段图
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
s=$SAMPLE; sex=$(sample_sex)
O=work/$s/str; mkdir -p $O/reviewer
CAT=ref/str/stranger_catalog_grch38_chr.json
if [ ! -s $CAT ]; then  # stranger 目录用 "19:…" 命名，参考是 "chr19"
python3 - "$CAT" <<'PY'
import json, sys
d = json.load(open("ref/str/stranger_variant_catalog_grch38.json"))
fix = lambda r: r if r.startswith("chr") else "chr" + r
for x in d:
    for k in ("ReferenceRegion", "OfftargetRegions"):
        if k in x:
            x[k] = [fix(r) for r in x[k]] if isinstance(x[k], list) else fix(x[k])
json.dump(d, open(sys.argv[1], "w"), indent=1)
PY
fi
$D $IMG_EH ExpansionHunter --reads work/$s/$s.bam \
  --reference $REF --variant-catalog $CAT --output-prefix $O/$s.eh --sex $sex --threads 8 > $O/$s.eh.log 2>&1
$D $IMG_SAMTOOLS sh -c "samtools sort -o $O/$s.eh_realigned.sorted.bam $O/$s.eh_realigned.bam && samtools index $O/$s.eh_realigned.sorted.bam"
$D $IMG_STRANGER sh -c "stranger -f ref/str/stranger_variant_catalog_grch38.json $O/$s.eh.vcf > $O/$s.eh.stranger.vcf" 2> $O/stranger.log
# 逐个位点画图：某个位点没分型（深度不够）时 REViewer 会报错退出，不能连累其他位点
loci=$(python3 -c 'import json,sys; print(" ".join(x["LocusId"] for x in json.load(open(sys.argv[1]))))' $CAT)
cat > $O/reviewer.sh <<SH
for l in $loci; do
  REViewer --reads $O/$s.eh_realigned.sorted.bam --vcf $O/$s.eh.vcf --reference ref/GRCh38/GRCh38.fa \
    --catalog $CAT --locus \$l --output-prefix $O/reviewer/$s || echo "REViewer failed: \$l"
done
SH
$D $IMG_REVIEWER sh $O/reviewer.sh > $O/reviewer.log 2>&1
$PY scripts/interpret/str_table.py $s > $O/$s.str_table.log 2>&1 || true
$PY scripts/interpret/rfc1_motif.py $s > $O/$s.rfc1.txt 2>&1 || true
echo "STR DONE $s ($(ls $O/reviewer/*.svg 2>/dev/null | wc -l) plots, $(grep -c 'REViewer failed' $O/reviewer.log) loci without plot)"
