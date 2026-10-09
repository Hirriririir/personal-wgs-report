"""BAM 局部深度小工具（只读索引覆盖的区域，很快）：常染色体基准深度、区域平均深度。"""
import os
import random

import numpy as np
import pysam


def read_filter(minq):
    def ok(r):
        return not (r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_unmapped) and r.mapping_quality >= minq
    return ok


def depth(bam, chrom, start, end, minq=20, baseq=0):
    cov = bam.count_coverage(chrom, start, end, quality_threshold=baseq, read_callback=read_filter(minq))
    return np.sum(cov, axis=0)


def autosomal_points(fa, n=3000, seed=11):
    """常染色体上 n 个随机点：优先用参考面板的常见 SNP（可比对、非 N），没有就在非 N 区随机取。"""
    rnd = random.Random(seed)
    p = "work/ext/panel/sites.tsv"
    if os.path.exists(p):
        sites = [x.split()[:2] for x in open(p)]
        return [(c, int(q)) for c, q in rnd.sample(sites, n)]
    chroms = [f"chr{i}" for i in range(1, 23)]
    lens = [fa.get_reference_length(c) for c in chroms]
    out = []
    while len(out) < n:
        c = rnd.choices(chroms, weights=lens)[0]
        q = rnd.randrange(1_000_000, fa.get_reference_length(c) - 1_000_000)
        seq = fa.fetch(c, q - 100, q + 100).upper()
        if "N" not in seq:
            out.append((c, q))
    return out


def autosomal_depth(bam, fa, minq=20, baseq=0, n=3000):
    pts = autosomal_points(fa, n)
    return float(np.mean(np.concatenate([depth(bam, c, q - 100, q + 100, minq, baseq) for c, q in pts])))


def open_bam(s):
    return pysam.AlignmentFile(f"work/{s}/{s}.bam")
