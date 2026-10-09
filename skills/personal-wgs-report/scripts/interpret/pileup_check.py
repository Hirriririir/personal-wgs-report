#!/usr/bin/env python3
"""单个位点的读段核查：每个样本列出支持参考 / 变异的读段数、正反链、比对质量、碱基质量、读段内位置、是否重复。
用法：pileup_check.py chr pos [sample ...]"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import collections
import statistics
import sys

import pysam

chrom, pos = sys.argv[1], int(sys.argv[2])
samples = sys.argv[3:] or [wgs.SAMPLE]
fa = pysam.FastaFile("ref/GRCh38/GRCh38.fa")
print(f"{chrom}:{pos} ref={fa.fetch(chrom, pos - 1, pos)}  context={fa.fetch(chrom, pos - 16, pos - 1)}[{fa.fetch(chrom, pos - 1, pos)}]{fa.fetch(chrom, pos, pos + 15)}")
for s in samples:
    bam = pysam.AlignmentFile(f"work/{s}/{s}.bam")
    stat = collections.defaultdict(lambda: {"n": 0, "fwd": 0, "mapq": [], "bq": [], "cyc": [], "dup": 0, "clip": 0})
    for col in bam.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=0, stepper="nofilter", ignore_orphans=False):
        for pr in col.pileups:
            r = pr.alignment
            if pr.is_del or pr.is_refskip or r.is_secondary or r.is_supplementary:
                continue
            b = r.query_sequence[pr.query_position]
            d = stat[b]
            d["n"] += 1
            d["fwd"] += 0 if r.is_reverse else 1
            d["mapq"].append(r.mapping_quality)
            d["bq"].append(r.query_qualities[pr.query_position])
            d["cyc"].append(pr.query_position if not r.is_reverse else r.query_length - 1 - pr.query_position)
            d["dup"] += r.is_duplicate
            d["clip"] += any(op == 4 for op, _ in (r.cigartuples or []))
    print(f"  {s}:")
    for b, d in sorted(stat.items(), key=lambda x: -x[1]["n"]):
        print(f"    {b}: n={d['n']} (dup {d['dup']}) fwd/rev={d['fwd']}/{d['n'] - d['fwd']} "
              f"MAPQ med={statistics.median(d['mapq']):.0f} min={min(d['mapq'])} BQ med={statistics.median(d['bq']):.0f} "
              f"readpos med={statistics.median(d['cyc']):.0f} softclipped={d['clip']}")
