#!/usr/bin/env bash
# SpliceAI：只算“解读范围内基因（GenCC 致病基因 + ACMG SF + 神经肌肉面板，约 5000 个）上的罕见变异（人群频率 <1%）”，
# 结果以 INFO/SpliceAI 合并回注释后的 VCF → work/annot/cohort.vep.spliceai.vcf.gz。
# 官方预计算分数要 Illumina 账号，这里自己算：有 GPU 约 1 小时，纯 CPU 要几个小时（可以跳过，后面的解读会少一类剪接预测）。
# 用法：scripts/10b_spliceai.sh [并行份数，默认 GPU 3 份 / CPU 按线程数]
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
if use_gpu; then N=${1:-3}; else N=${1:-$(( THREADS / 4 > 1 ? THREADS / 4 : 1 ))}; export CUDA_VISIBLE_DEVICES=""; fi
O=work/annot/spliceai; mkdir -p $O
BCF="$D $IMG_BCFTOOLS"
# 1) 基因区间（GENCODE v50，两侧各 50 bp）
if [ ! -s $O/genes.bed ]; then
  cut -f1 ref/genelists/gene_table.tsv | tail -n +2 | sort -u > $O/genes.txt
  zcat ref/annot/gencode/gencode.v50.annotation.gtf.gz | awk -F'\t' '$3=="gene"' \
    | $PY -c '
import sys, re
keep = set(open(sys.argv[1]).read().split())
for line in sys.stdin:
    f = line.split("\t")
    g = re.search(r"gene_name \"([^\"]+)\"", f[8]).group(1)
    if g in keep and f[0] != "chrM":
        print(f"{f[0]}\t{max(0, int(f[3]) - 51)}\t{int(f[4]) + 50}\t{g}")
' $O/genes.txt | sort -k1,1V -k2,2n > $O/genes.bed
fi
# 2) 罕见变异（只保留位点，不要样本）
$PY - <<'PY'
from cyvcf2 import VCF
import re
v = VCF("work/annot/cohort.vep.vcf.gz")
fmt = re.search(r"Format: ([^\"]+)", v.get_header_type("CSQ")["Description"]).group(1).split("|")
i_af = fmt.index("MAX_AF")
regions = [l.split()[:3] for l in open("work/annot/spliceai/genes.bed")]
out = open("work/annot/spliceai/rare.sites.vcf", "w")
out.write("##fileformat=VCFv4.2\n")
for c in [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]:
    out.write(f"##contig=<ID={c}>\n")
out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
n = 0
for chrom, s, e in regions:
    for r in v(f"{chrom}:{int(s)+1}-{e}"):
        csq = r.INFO.get("CSQ") or ""
        afs = [float(x.split("|")[i_af]) for x in csq.split(",") if x.split("|")[i_af] not in ("",)]
        if afs and max(afs) >= 0.01:
            continue
        out.write(f"{r.CHROM}\t{r.POS}\t.\t{r.REF}\t{r.ALT[0]}\t.\t.\t.\n")
        n += 1
print("基因区间内的罕见变异（重叠基因会重复，下一步去重）：", n)
PY
$BCF sh -c "bcftools sort $O/rare.sites.vcf -Ou | bcftools norm -d exact -Oz -o $O/rare.sites.vcf.gz && bcftools index -f -t $O/rare.sites.vcf.gz"
# 3) 拆成 N 份并行（GPU 上每个进程按需占显存）
rm -f $O/chunk.*
$BCF sh -c "bcftools view -h $O/rare.sites.vcf.gz > $O/header.txt && bcftools view -H $O/rare.sites.vcf.gz > $O/body.txt"
log "SpliceAI：$(wc -l < $O/body.txt) 个罕见位点，分 $N 份"
split -n l/$N -d -a 2 $O/body.txt $O/body.
for b in $O/body.[0-9][0-9]; do cat $O/header.txt $b > $O/chunk.${b##*.}.vcf; done
rm -f $O/body.txt $O/body.[0-9][0-9]
for c in $O/chunk.[0-9][0-9].vcf; do
  TF_FORCE_GPU_ALLOW_GROWTH=true TF_CPP_MIN_LOG_LEVEL=2 tools/venv-annot/bin/spliceai -I $c -O ${c%.vcf}.out.vcf \
    -R $REF -A grch38 -D 50 > ${c%.vcf}.log 2>&1 &
done
wait
# 4) 合并并写回
$BCF sh -c "for f in $O/chunk.*.out.vcf; do bcftools sort \$f -Oz -o \$f.gz && bcftools index -f -t \$f.gz; done && \
  bcftools concat -a $O/chunk.*.out.vcf.gz -Oz -o $O/rare.spliceai.vcf.gz && bcftools index -f -t $O/rare.spliceai.vcf.gz && \
  bcftools annotate -a $O/rare.spliceai.vcf.gz -c INFO/SpliceAI work/annot/cohort.vep.vcf.gz -Oz -o work/annot/cohort.vep.spliceai.vcf.gz && \
  bcftools index -f -t work/annot/cohort.vep.spliceai.vcf.gz"
log "SPLICEAI DONE"
