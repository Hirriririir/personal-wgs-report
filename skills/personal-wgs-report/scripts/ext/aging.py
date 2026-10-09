#!/usr/bin/env python3
"""衰老相关指标（只读 BAM 的局部区域）→ work/ext/aging/aging.<s>.json
  1) 线粒体 DNA 拷贝数 = 2 × chrM 平均深度 / 常染色体平均深度（MAPQ≥20、碱基质量≥20、去重复）
  2) 男性：Y 染色体相对深度（Yleaf 标记位点：单拷贝、可比对）；没有嵌合性 Y 丢失（LOY）时期望 ≈0.5
  3) 克隆性造血（CHIP）：25 个常见基因的编码区逐位点计数，找 VAF 3–40%、≥3 条支持读段的“非胚系样”变异
     30× 全基因组只能看到 VAF 较高的克隆，只作筛查提示
  4) 端粒长度：TelSeq 估计值（x05_aging.sh 先跑），读这里汇总
用法：aging.py <sample>"""
import glob
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import wgs  # noqa: E402
import depth_lib as dl  # noqa: E402
import numpy as np  # noqa: E402
import pysam  # noqa: E402

s = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
O = "work/ext/aging"
os.makedirs(O, exist_ok=True)
bam = dl.open_bam(s)
fa = pysam.FastaFile(wgs.REF)
auto = dl.autosomal_depth(bam, fa, minq=20, baseq=20)
mt = float(np.mean(dl.depth(bam, "chrM", 0, 16569, 20, 20)))
res = {"sample": s, "auto_mean_depth_q20": round(auto, 2), "chrM_mean_depth_q20": round(mt, 1), "mtDNA_copies_per_cell": round(2 * mt / auto, 1)}
if wgs.sex_of(s) == "male":
    marks = []
    for f in sorted(glob.glob(f"work/{s}/ychr/yleaf/*/*.out"), key=lambda f: "ftdna" not in f):
        for line in open(f):
            x = line.split("\t")
            if x[0].lower() == "chry" and len(x) > 7 and x[7].isdigit():
                marks.append(int(x[1]))
        if marks:
            break
    if not marks:   # 没有 Yleaf 结果：用 chrY 雄性特异区（X 退化区）里随机点
        rnd = random.Random(3)
        marks = [rnd.randrange(2_800_000, 26_600_000) for _ in range(3000)]
    pts = random.Random(11).sample(marks, min(3000, len(marks)))
    yd = float(np.mean(np.concatenate([dl.depth(bam, "chrY", p - 100, p + 100, 20, 20) for p in pts])))
    res.update(chrY_marker_mean_depth_q20=round(yd, 2), Y_to_auto_ratio=round(yd / auto, 3), LOY_fraction_est=round(max(0.0, 1 - 2 * yd / auto), 3))
# CHIP
germline = set()
vf = pysam.VariantFile(f"work/{s}/{s}.dv.vcf.gz")
hits, scanned = [], 0
for line in open("ref/annot/gencode/chip_cds.bed"):
    c, st, en, g = line.split()
    st, en = int(st), int(en)
    scanned += en - st
    for rec in vf.fetch(c, st, en):
        germline.add((c, rec.pos))
    cov = bam.count_coverage(c, st, en, quality_threshold=20, read_callback=dl.read_filter(20))
    ref = fa.fetch(c, st, en).upper()
    for i in range(en - st):
        counts = {b: cov[k][i] for k, b in enumerate("ACGT")}
        dp = sum(counts.values())
        if dp < 8:
            continue
        for b, n in counts.items():
            if b == ref[i] or n < 3:
                continue
            vaf = n / dp
            if 0.03 <= vaf <= 0.40:
                hits.append(dict(gene=g, chrom=c, pos=st + i + 1, ref=ref[i], alt=b, alt_reads=n, depth=dp, vaf=round(vaf, 3),
                                 germline_call=(c, st + i + 1) in germline))
res["chip_scanned_bp"] = scanned
res["chip_candidates"] = hits
t = f"{O}/telseq.{s}.txt"
if os.path.exists(t):
    rows = [x.rstrip("\n").split("\t") for x in open(t) if x.strip()]
    if len(rows) > 1:
        d = dict(zip(rows[0], rows[1]))
        res["telseq_length_kb"] = round(float(d["LENGTH_ESTIMATE"]), 2)
json.dump(res, open(f"{O}/aging.{s}.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "chip_candidates"}, ensure_ascii=False))
print(f"CHIP 候选 {len(hits)} 个（germline_call=True 的多为胚系杂合或比对伪影）")
for h in hits:
    print(h)
