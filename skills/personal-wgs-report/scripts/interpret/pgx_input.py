#!/usr/bin/env python3
"""给 PharmCAT 生成输入 VCF：
  - DeepVariant 在药物基因区域里检出的真实变异原样保留；
  - PharmCAT 位点表里、本人没有变异且落在高质量参考区块（GQ≥20）的位点，按 PharmCAT 的 REF/ALT 写一条 0/0
    （插入缺失位点必须带完整 REF，单碱基 REF + ALT=. 会被 PharmCAT 判为无效而当缺失）；
  - 覆盖不足 / 低质量的位点不写，PharmCAT 会如实报“缺失”。
用法：pgx_input.py <sample> <pharmcat_positions.vcf.bgz> <out.vcf>
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import sys

import pysam

s, pos_path, out_path = sys.argv[1:4]
dv = pysam.VariantFile(f"work/{s}/{s}.dv.vcf.gz")
gv = pysam.VariantFile(f"work/{s}/{s}.dv.g.vcf.gz")
pc = pysam.VariantFile(pos_path)

def confident_ref(chrom, start, end):
    """[start,end) 全部落在 GQ≥20 的纯合参考区块里，且没有任何非参考记录。"""
    covered = start
    for rec in gv.fetch(chrom, start, end):
        alts = rec.alts or ()
        smp = rec.samples[0]
        if alts not in ((), ("<*>",)):
            gt = smp.get("GT")
            if gt and any(a not in (0, None) for a in gt):
                return False
        gq = smp.get("GQ") or 0
        if gq < 20:
            return False
        if rec.start > covered:
            return False
        covered = max(covered, rec.stop)
    return covered >= end

lines = []
header = str(dv.header).rstrip("\n").split("\n")
seen = set()
for p in pc.fetch():
    chrom = p.chrom if p.chrom.startswith("chr") else "chr" + p.chrom
    start, end = p.start, p.stop
    # 本人在这个位点附近的真实变异
    var = [r for r in dv.fetch(chrom, max(0, start - 1), end + 1) if r.filter.keys() in (["PASS"], [])]
    for r in var:
        k = (r.chrom, r.pos, r.ref, r.alts)
        if k not in seen:
            seen.add(k)
            lines.append((r.chrom, r.pos, str(r).rstrip("\n")))
    if var:
        continue
    if confident_ref(chrom, start, end):
        alts = ",".join(a for a in (p.alts or ()) if a)
        ref = p.ref
        fmt_cols = "GT"
        lines.append((chrom, p.pos, f"{chrom}\t{p.pos}\t{p.id or '.'}\t{ref}\t{alts or '.'}\t.\tPASS\t.\t{fmt_cols}\t0/0"))
order = {c: i for i, c in enumerate(dv.header.contigs)}
lines.sort(key=lambda x: (order.get(x[0], 10**6), x[1]))
with open(out_path, "w") as o:
    for h in header:
        o.write(h + "\n")
    for _, _, l in lines:
        o.write(l + "\n")
print(f"{s}: wrote {len(lines)} records ({len(seen)} variant records)")
