#!/usr/bin/env bash
# 没比对上人类参考的读段是什么（约 0.1–0.5%）：未比对读段 + 比对到 EBV 诱饵序列（chrEBV）的读段 → Kraken2（standard-8GB 库）逐条分类
# 能看到：血液里常见的细环病毒（TTV）、潜伏的 EBV / 疱疹病毒、试剂盒污染菌、参考基因组没收录的人类序列等。
# 只能说明“测到了什么”，测不到不等于没感染（血浆 / 血细胞里的病毒 DNA 量通常很低）。
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
s=$SAMPLE; O=work/ext/virome; mkdir -p $O
if [ ! -s $O/$s.k2_1.fq.gz ]; then
  $D $IMG_SAMTOOLS sh -c "samtools view -@4 -b -f 4 work/$s/$s.bam > $O/$s.unmapped.bam \
    && samtools view -@4 -b work/$s/$s.bam chrEBV > $O/$s.chrEBV.bam \
    && samtools merge -f -@4 $O/$s.k2in.bam $O/$s.unmapped.bam $O/$s.chrEBV.bam \
    && samtools fastq -@4 -1 $O/$s.k2_1.fq.gz -2 $O/$s.k2_2.fq.gz -s $O/$s.k2_s.fq.gz -0 /dev/null -n $O/$s.k2in.bam" > $O/$s.prep.log 2>&1
fi
$D $IMG_KRAKEN2 kraken2 --db ref/kraken2/std8 --threads $(( THREADS > 8 ? 8 : THREADS )) --report $O/$s.k2.report --output $O/$s.k2.out \
  --use-names $O/$s.k2_1.fq.gz $O/$s.k2_2.fq.gz $O/$s.k2_s.fq.gz > $O/$s.k2.log 2>&1
$PYB - "$s" <<'PY'
import json, sys
s = sys.argv[1]
rep = [l.rstrip("\n").split("\t") for l in open(f"work/ext/virome/{s}.k2.report")]
clade = {r[5].strip(): int(r[1]) for r in rep}
rank = {r[5].strip(): r[3] for r in rep}
total = clade.get("unclassified", 0) + clade.get("root", 0)
viral = sorted(([n, r[3], int(r[1])] for r in rep for n in [r[5].strip()] if r[3] in ("S", "S1", "G", "F") and int(r[1]) > 0
                and any(k in n.lower() for k in ("virus", "viridae", "phage"))), key=lambda x: -x[2])   # [名称, 分类级别, 读段数]
top_bact = sorted(((r[5].strip(), int(r[1])) for r in rep if r[3] == "G"), key=lambda x: -x[1])[:10]
out = {"total_reads": total, "unclassified": clade.get("unclassified", 0), "bacteria": clade.get("Bacteria", 0),
       "archaea": clade.get("Archaea", 0), "viruses": clade.get("Viruses", 0), "human": clade.get("Homo sapiens", 0),
       "viral_taxa": viral[:30], "top_genera": top_bact}
json.dump(out, open(f"work/ext/virome/{s}.summary.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(out, ensure_ascii=False, indent=1)[:2000])
PY
log "VIROME DONE $s"
