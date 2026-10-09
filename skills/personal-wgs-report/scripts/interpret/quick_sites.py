#!/usr/bin/env python3
"""临时查几个 rsID 在某人 gVCF 里的基因型（Ensembl 查 GRCh38 坐标，结果缓存到 tmp/）。
用法：quick_sites.py <sample> rsid[:标签] ..."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import json
import os
import sys
import time
import urllib.request

import pysam

s, *items = sys.argv[1:]
cache_p = "tmp/quick_sites_cache.json"
cache = json.load(open(cache_p)) if os.path.exists(cache_p) else {}

def ensembl(rs):
    if rs in cache:
        return cache[rs]
    url = f"https://rest.ensembl.org/variation/human/{rs}?content-type=application/json"
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.load(r)
            for m in d.get("mappings", []):
                if m.get("assembly_name") == "GRCh38" and m.get("seq_region_name") in [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]:
                    cache[rs] = ["chr" + m["seq_region_name"] if m["seq_region_name"] != "MT" else "chrM", int(m["start"]), int(m["end"]), m["allele_string"]]
                    return cache[rs]
        except Exception:  # noqa: BLE001
            time.sleep(2)
    cache[rs] = None
    return None

vf = pysam.VariantFile(f"work/{s}/{s}.dv.g.vcf.gz")
fa = pysam.FastaFile("ref/GRCh38/GRCh38.fa")
for it in items:
    rs, _, label = it.partition(":")
    m = ensembl(rs)
    if not m:
        print(f"{label or rs}\t{rs}\t坐标未知")
        continue
    chrom, start, end, alleles = m
    ref = fa.fetch(chrom, start - 1, end).upper() if start <= end else "-"
    got = "未覆盖"
    for rec in vf.fetch(chrom, max(0, start - 2), end + 1):
        smp = rec.samples[0]
        alts = rec.alts or ()
        if alts in ((), ("<*>",)):
            if rec.start <= start - 1 < (rec.stop or rec.start + 1):
                got = f"{ref}/{ref}（参考纯合，GQ {smp.get('GQ')}）"
            continue
        gt = smp.get("GT")
        if gt is None or None in gt or rec.pos != start:
            continue
        al = [rec.alleles[i] for i in gt]
        got = f"{'/'.join(al)}（变异，GQ {smp.get('GQ')}，深度 {smp.get('DP')}）"
        break
    print(f"{label or rs}\t{rs}\t{chrom}:{start}\t{alleles}\t{got}")
json.dump(cache, open(cache_p, "w"))
