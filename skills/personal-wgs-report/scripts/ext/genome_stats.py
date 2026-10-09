#!/usr/bin/env python3
"""基因组概览：PASS 变异统计、与 1000 Genomes 3,202 人对比的“新”变异（只数 1000G 严格可检测区里 GQ≥20 的）、
1 Mb 窗口杂合 SNV 密度 → work/ext/genome_stats.<s>.json + work/ext/het_density.<s>.tsv
用法：genome_stats.py <sample>"""
import collections
import io
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402
import pysam  # noqa: E402
import zstandard  # noqa: E402

s_ = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
vcf = pysam.VariantFile(f"work/{s_}/{s_}.dv.vcf.gz")
fa = pysam.FastaFile(wgs.REF)
TI = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def pvar_by_chrom():
    """按染色体依次产出 (chrom, set((pos, ref, alt)))；pvar 只有 5 列"""
    with open("ref/pca/all_hg38.pvar.zst", "rb") as fh:
        rd = io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(fh), encoding="utf-8")
        cur, s = None, set()
        for line in rd:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != cur:
                if cur is not None:
                    yield cur, s
                cur, s = f[0], set()
            p = int(f[1])
            for a in f[4].split(","):
                s.add((p, f[3], a))
        yield cur, s


# 1000G 严格可检测区（hmmix 提供的 hg38 strict mask）
import bisect
MASK = collections.defaultdict(lambda: ([], []))
for line in open("ref/hmmix/hg38_strick_callability_mask.bed"):
    c, a, b = line.split()[:3]
    MASK[c][0].append(int(a))
    MASK[c][1].append(int(b))


def in_mask(c, pos):
    st, en = MASK[c]
    i = bisect.bisect_right(st, pos - 1) - 1
    return i >= 0 and pos - 1 < en[i]


chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
st = collections.Counter()
het_bins = collections.Counter()
gen = pvar_by_chrom()
pending = {}
for c in chroms:
    kc = c[3:]
    while kc not in pending:
        k, sset = next(gen)
        pending[k] = sset
    s = pending.pop(kc)
    for rec in vcf.fetch(c):
        if rec.filter.keys() not in ([], ["PASS"]):
            continue
        smp = rec.samples[0]
        gt = smp["GT"]
        if gt is None or None in gt or max(gt) == 0:
            continue
        gq = smp.get("GQ") or 0
        hiconf = gq >= 20 and in_mask(c, rec.pos)
        for a in [rec.alleles[i] for i in set(gt) if i > 0]:
            kind = "SNV" if len(rec.ref) == 1 and len(a) == 1 else "indel"
            het = gt[0] != gt[1]
            st[kind] += 1
            st[f"{kind}_{'het' if het else 'hom'}"] += 1
            if kind == "SNV":
                st["ti" if (rec.ref, a) in TI else "tv"] += 1
                if het and gq >= 20 and c != "chrX":
                    het_bins[(c, rec.pos // 1_000_000)] += 1
            if hiconf:
                st[f"{kind}_hiconf"] += 1
                if (rec.pos, rec.ref, a) not in s:
                    st[f"{kind}_hiconf_not_in_1kg"] += 1
                    st[f"{kind}_hiconf_not_in_1kg_{'het' if het else 'hom'}"] += 1
    print(c, dict(st), flush=True)
    del s
st["titv"] = round(st["ti"] / max(1, st["tv"]), 3)
os.makedirs("work/ext", exist_ok=True)
json.dump(dict(st), open(f"work/ext/genome_stats.{s_}.json", "w"), indent=1)
with open(f"work/ext/het_density.{s_}.tsv", "w") as o:
    o.write("chrom\tbin_mb\thet_snv\tnonN_bases\n")
    for c in chroms[:-1]:
        L = fa.get_reference_length(c)
        for b in range(L // 1_000_000 + 1):
            seq = fa.fetch(c, b * 1_000_000, min(L, (b + 1) * 1_000_000))
            nonN = len(seq) - seq.upper().count("N")
            o.write(f"{c}\t{b}\t{het_bins.get((c, b), 0)}\t{nonN}\n")
print("DONE", json.dumps(st))
