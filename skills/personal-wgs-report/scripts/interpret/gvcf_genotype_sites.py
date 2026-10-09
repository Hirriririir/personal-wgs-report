#!/usr/bin/env python3
"""在指定的双等位 SNP 位点上，从 DeepVariant gVCF 取每个人的基因型，输出 REF/ALT 与位点表一致的 VCF。

gVCF 里没写出来的位点要么在参考区块（纯合参考）里，要么没覆盖；只有区块 GQ≥20 才判 0/0，否则记缺失。
这样 PCA 投影、PRS 打分时纯合参考是明确的 0/0，而不是缺失或等位基因对不上。

用法：gvcf_genotype_sites.py sites.tsv(chrom pos id ref alt，chrom 可带可不带 chr) out.vcf.gz name1=g1.vcf.gz [name2=g2.vcf.gz ...]
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import sys
import pysam

sites_path, out_path, *pairs = sys.argv[1:]
samples = [p.split("=", 1) for p in pairs]
gv = [pysam.VariantFile(path) for _, path in samples]

def call(vf, chrom, pos, ref, alt):
    for rec in vf.fetch(chrom, pos - 1, pos):
        smp = rec.samples[0]
        gt = smp.get("GT")
        if rec.alts is None or rec.alts == ("<*>",) or (rec.alts and rec.alts[0] == "<*>" and len(rec.alts) == 1):
            gq = smp.get("GQ") or 0
            return "0/0" if gq >= 20 else "./."
        if rec.pos != pos or rec.ref[0] != ref:
            continue  # 上游 indel 覆盖到这里，跳过
        if gt is None or None in gt:
            return "./."
        alleles = [rec.alleles[i] for i in gt]
        if any(a not in (ref, alt) for a in alleles if len(a) == 1) or any(len(a) != 1 for a in alleles if a != "<*>"):
            return "./."  # 多等位 / 不同的 ALT / indel
        n_alt = sum(a == alt for a in alleles)
        return {0: "0/0", 1: "0/1", 2: "1/1"}[n_alt]
    return "./."

header = ["##fileformat=VCFv4.2", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
chroms = []
rows = []
with open(sites_path) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        chrom, pos, vid, ref, alt = line.split()[:5]
        c = chrom if chrom.startswith("chr") else "chr" + chrom
        pos = int(pos)
        gts = [call(v, c, pos, ref, alt) for v in gv]
        rows.append(f"{c}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts))
        if c not in chroms:
            chroms.append(c)
header += [f"##contig=<ID={c}>" for c in chroms]
header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(n for n, _ in samples))
tmp = out_path[:-3] if out_path.endswith(".gz") else out_path
with open(tmp, "w") as o:
    o.write("\n".join(header) + "\n" + "\n".join(rows) + "\n")
if out_path.endswith(".gz"):
    pysam.tabix_compress(tmp, out_path, force=True)
    pysam.tabix_index(out_path, preset="vcf", force=True)
miss = [sum(r.split("\t")[9 + i] == "./." for r in rows) for i in range(len(samples))]
print(f"{len(rows)} sites; missing per sample: " + ", ".join(f"{n}={m}" for (n, _), m in zip(samples, miss)))
