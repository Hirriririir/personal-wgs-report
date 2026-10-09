#!/usr/bin/env python3
"""顺序扫描版（按染色体并行）：在大量位点（上千万）上从 gVCF 取基因型。位点表与 gVCF 都按坐标排序，双指针一遍扫完。
规则同 gvcf_genotype_sites.py：参考区块 GQ≥20 记 0/0；变异记录按 REF/ALT 匹配；否则缺失。
用法：gvcf_genotype_stream_par.py sites.tsv(chrom pos id ref alt；已按染色体 1..22 / 位置排序) out.vcf.gz name=g.vcf.gz [...]
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import sys
from cyvcf2 import VCF
import pysam

sites_path, out_path, *pairs = sys.argv[1:]
samples = [p.split("=", 1) for p in pairs]

# 读位点（按染色体分组）
by_chr = {}
order = []
with open(sites_path) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        c, p, i, r, a = line.split()[:5]
        c = c if c.startswith("chr") else "chr" + c
        if c not in by_chr:
            by_chr[c] = []
            order.append(c)
        by_chr[c].append((int(p), i, r, a))
for c in by_chr:
    by_chr[c].sort()

def genotype_chrom(args):
    path, c = args
    res = {c: ["./."] * len(by_chr[c])}
    vf = VCF(path)
    for c in [c]:
        sites = by_chr[c]
        n = len(sites)
        j = 0
        try:
            it = vf(c)
        except Exception:  # noqa: BLE001
            continue
        for rec in it:
            if j >= n:
                break
            start = rec.POS            # 1-based
            alts = rec.ALT
            is_block = (not alts) or alts == ["<*>"]
            end = rec.INFO.get("END") if is_block else rec.POS + len(rec.REF) - 1
            end = int(end) if end else start
            while j < n and sites[j][0] < start:
                j += 1
            if j >= n:
                break
            gq = rec.format("GQ")
            gqv = int(gq[0][0]) if gq is not None else 0
            if is_block:
                if gqv >= 20:
                    k = j
                    while k < n and sites[k][0] <= end:
                        res[c][k] = "0/0"
                        k += 1
                continue
            # 变异记录：只处理位点正好落在 POS 且是 SNV 的情况
            k = j
            while k < n and sites[k][0] == start:
                pos, sid, ref, alt = sites[k]
                gt = rec.genotypes[0]
                if len(rec.REF) == 1 and rec.REF == ref and gt[0] >= 0 and gt[1] >= 0:
                    al = [rec.REF] + alts
                    ga = [al[x] if x < len(al) else "?" for x in gt[:2]]
                    if all(g in (ref, alt) for g in ga):
                        na = sum(g == alt for g in ga)
                        res[c][k] = ("0/0", "0/1", "1/1")[na]
                    elif any(g == "<*>" for g in ga) and all(g in (ref, "<*>") for g in ga):
                        res[c][k] = "0/0"
                k += 1
            # 覆盖到后面位点的缺失：这些位点设缺失（保持 ./.）
    return (path, c, res[c])

from multiprocessing import Pool
with Pool(16) as pool:
    parts = pool.map(genotype_chrom, [(p, c) for _, p in samples for c in order], chunksize=1)
gts = []
for _, p in samples:
    gts.append({c: r for (pp, c, r) in parts if pp == p})
tmp = out_path[:-3]
with open(tmp, "w") as o:
    o.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    for c in order:
        o.write(f"##contig=<ID={c}>\n")
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(n for n, _ in samples) + "\n")
    for c in order:
        for k, (p, i, r, a) in enumerate(by_chr[c]):
            o.write(f"{c}\t{p}\t{i}\t{r}\t{a}\t.\tPASS\t.\tGT\t" + "\t".join(g[c][k] for g in gts) + "\n")
pysam.tabix_compress(tmp, out_path, force=True)
import os
os.remove(tmp)
tot = sum(len(v) for v in by_chr.values())
miss = [sum(x == "./." for c in order for x in g[c]) for g in gts]
print(f"{tot} sites; missing: " + ", ".join(f"{n}={m} ({m / max(1, tot):.1%})" for (n, _), m in zip(samples, miss)))
