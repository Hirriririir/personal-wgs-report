#!/usr/bin/env python3
"""某人在注释区域内的变异统计 + 纯合功能缺失（“人类敲除”）基因清单。用法：knockouts.py <sample>"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import collections
import re
import sys

from cyvcf2 import VCF

s = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
v = VCF(wgs.annot_vcf())
fmt = re.search(r"Format: ([^\"]+)", v.get_header_type("CSQ")["Description"]).group(1).split("|")
I = {k: i for i, k in enumerate(fmt)}
si = v.samples.index(s)
LOF = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}
n = het = hom = norsid = 0
ko = collections.defaultdict(list)
hetlof = set()
for r in v:
    if r.FILTER not in (None, "PASS"):
        continue
    g = r.genotypes[si]
    if 1 not in g[:2]:
        continue
    n += 1
    if g[0] == g[1] == 1:
        hom += 1
    else:
        het += 1
    csq = r.INFO.get("CSQ") or ""
    ex = [e.split("|") for e in csq.split(",")]
    if not any(x[I["Existing_variation"]] for x in ex):
        norsid += 1
    gq = r.format("GQ")[si][0]
    for x in ex:
        if x[I["SOURCE"]] != "Ensembl" or not x[I["CANONICAL"]]:
            continue
        cons = set(x[I["Consequence"]].split("&"))
        if cons & LOF and x[I["BIOTYPE"]] == "protein_coding" and gq >= 20:
            lab = x[I["HGVSp"]].split(":")[-1] or x[I["HGVSc"]].split(":")[-1]
            if g[0] == g[1] == 1:
                ko[x[I["SYMBOL"]]].append(lab)
            else:
                hetlof.add(x[I["SYMBOL"]])
print(f"{s} 注释区域（外显子±50bp+ClinVar位点）内：变异 {n}，杂合 {het}，纯合 {hom}，无 rsID（可能是新变异）{norsid}")
print(f"纯合功能缺失基因（敲除）{len(ko)} 个，杂合功能缺失基因 {len(hetlof)} 个")
ors = sorted(g for g in ko if g.startswith(("OR", "TAS2R")))
print(f"  其中嗅觉/味觉受体 {len(ors)} 个: " + ", ".join(ors[:30]))
rest = [f"{g}({ko[g][0]})" for g in sorted(ko) if g not in ors]
print("  非受体类 (%d): " % len(rest) + ", ".join(rest))
