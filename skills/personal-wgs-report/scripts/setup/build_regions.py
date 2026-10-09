#!/usr/bin/env python3
"""注释范围：GENCODE 全部外显子两侧各 50 bp + ClinVar 收录的所有位点，合并成 ref/annot/regions/exons50_clinvar.bed。
VEP 只注释这些区域（全基因组注释太慢，也用不上）。"""
import gzip
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401

iv = {}
with gzip.open("ref/annot/gencode/gencode.v50.annotation.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t", 5)
        if p[2] != "exon":
            continue
        iv.setdefault(p[0], []).append((max(0, int(p[3]) - 51), int(p[4]) + 50))
with gzip.open("ref/annot/clinvar/clinvar.chr.vcf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c, pos, _, ref = line.split("\t", 4)[:4]
        iv.setdefault(c, []).append((int(pos) - 1, int(pos) - 1 + max(1, len(ref))))
os.makedirs("ref/annot/regions", exist_ok=True)
order = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y", "M"]]
n = 0
with open("ref/annot/regions/exons50_clinvar.bed", "w") as o:
    for c in sorted(iv, key=lambda x: order.index(x) if x in order else 99):
        cur = None
        for s, e in sorted(iv[c]):
            if cur and s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                if cur:
                    o.write(f"{c}\t{cur[0]}\t{cur[1]}\n"); n += 1
                cur = [s, e]
        if cur:
            o.write(f"{c}\t{cur[0]}\t{cur[1]}\n"); n += 1
print(f"{n} merged regions")
