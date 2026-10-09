#!/usr/bin/env python3
"""古人类渗入：生成 hmmix 的 ingroup 观测文件（等价于 hmmix create_ingroup，省掉 vcftools 依赖）。
规则同 hmmix 0.9.2：只取 SNP；在严格可检测区（strict mask）内；去掉外群（1000G+HGDP 非洲人）里出现过的位点；
有祖先等位基因信息时，只保留携带衍生等位基因（非祖先）的位点。输出 chrom pos ancestral_base genotype。
用法：hmmix_obs.py <in.vcf.gz> <out.txt（可含 {sample}）>"""
import bisect
import collections
import sys

import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录
import pysam  # noqa: E402

H = "ref/hmmix"
vcf_path, out_path = sys.argv[1], sys.argv[2]  # out_path 可含 {sample}，多样本 VCF 时每人一个文件

mask = collections.defaultdict(lambda: ([], []))
for line in open(f"{H}/hg38_strick_callability_mask.bed"):
    c, s, e = line.split()[:3]
    mask[c][0].append(int(s))
    mask[c][1].append(int(e))


def in_mask(c, pos):  # pos 1-based；BED 0-based 半开
    st, en = mask[c]
    i = bisect.bisect_right(st, pos - 1) - 1
    return i >= 0 and pos - 1 < en[i]


out_pos = collections.defaultdict(set)
with open(f"{H}/hg38_Outgroup_1000g_HGDP.txt") as f:
    f.readline()
    for line in f:
        c, p = line.split("\t", 2)[:2]
        out_pos[c].add(int(p))
print("outgroup positions loaded", sum(len(v) for v in out_pos.values()), file=sys.stderr)

vf = pysam.VariantFile(vcf_path)
samples = list(vf.header.samples)
outs = {sm: open(out_path.replace("{sample}", sm) if "{sample}" in out_path else out_path, "w") for sm in samples}
for o in outs.values():
    print("chrom", "pos", "ancestral_base", "genotype", sep="\t", file=o)
n = collections.Counter()
for chrom_n in range(1, 23):
    c = f"chr{chrom_n}"
    anc = "".join(l.strip().upper() for l in open(f"{H}/hg38_ancestral/homo_sapiens_ancestor_{chrom_n}.fa") if not l.startswith(">"))
    for rec in vf.fetch(c):
        if rec.filter.keys() not in ([], ["PASS"]):
            continue
        alleles = rec.alleles
        if len(alleles[0]) != 1 or alleles[0] not in "ACGT":
            continue
        if not in_mask(c, rec.pos):
            continue
        if rec.pos in out_pos[c]:
            continue
        a = anc[rec.pos - 1] if rec.pos - 1 < len(anc) else "N"
        if a not in alleles:
            continue
        for sm in samples:
            gt = rec.samples[sm]["GT"]
            if gt is None or None in gt:
                continue
            bases = [alleles[g] for g in gt]
            if not all(len(b) == 1 and b in "ACGT" for b in bases):
                continue
            g = "".join(bases)
            if g.count(a) != 2:
                print(c, rec.pos, a, g, sep="\t", file=outs[sm])
                n[sm] += 1
    print(c, "done", file=sys.stderr)
for o in outs.values():
    o.close()
print("written per sample", dict(n), file=sys.stderr)
