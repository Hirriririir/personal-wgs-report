#!/usr/bin/env python3
"""古人类片段汇总（hmmix 二倍体解码，后验 ≥0.8）→ work/ext/archaic/summary.json + segments.<s>.tsv
归类：片段内与 3 个尼安德特人（Altai / Vindija / Chagyrskaya）共享的衍生变异数取最大值 N，与丹尼索瓦人共享数 D；
N>D → 尼安德特型，D>N → 丹尼索瓦型，相等 → 难分，均为 0 → 未能归类。每段标注重叠的蛋白编码基因，并标出文献里的著名渗入基因。"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

D = "work/ext/archaic"
NEA = ["AltaiNeandertal", "Chagyrskaya-Phalanx", "Vindija33.19"]
KNOWN = ["OAS1", "OAS2", "OAS3", "BNC2", "SLC16A11", "SLC16A13", "LZTFL1", "EPAS1", "TLR1", "TLR6", "TLR10", "POU2F3", "HYAL2",
         "TBX15", "WARS2", "MUC19", "STAT2", "HMGCR", "CERT1", "KRT71", "KRT74", "KRT5", "HLA-A", "HLA-B", "HLA-C", "PNMA1",
         "ZNF365", "SIPA1L2", "GLI3", "CHMP1A", "DSTYK"]


def segments(path):
    d = pd.read_csv(path, sep="\t")
    a = d[(d["state"] == "Archaic") & (d["mean_prob"] >= 0.8)].copy()
    N = a[[c for c in NEA if c in a]].max(axis=1)
    Dn = a["Denisova"]
    a["cls"] = np.where((N == 0) & (Dn == 0), "unclassified", np.where(N > Dn, "Neanderthal", np.where(Dn > N, "Denisovan", "ambiguous")))
    return a


def totals(a):
    return dict(total_Mb=round(a["length"].sum() / 1e6, 2), nea_Mb=round(a.loc[a["cls"] == "Neanderthal", "length"].sum() / 1e6, 2),
                den_Mb=round(a.loc[a["cls"] == "Denisovan", "length"].sum() / 1e6, 2), n_segments=int(len(a)),
                median_kb=round(float(a["length"].median()) / 1e3, 1) if len(a) else 0)


genes = collections.defaultdict(list)
for line in open("ref/annot/gencode/pc_genes.bed"):
    c, s_, e, g = line.split()
    genes[c].append((int(s_), int(e), g))
ctrl = {}
pop = dict(x.split() for x in open(f"{D}/kg/ids.tsv")) if os.path.exists(f"{D}/kg/ids.tsv") else {}
for i, p in pop.items():
    f = f"{D}/kg/decoded.{i}.diploid.txt"
    if os.path.exists(f):
        ctrl[i] = dict(pop=p, **totals(segments(f)))
out = {"controls": ctrl, "samples": {}}
for s in wgs.all_samples():
    f = f"{D}/decoded.{s}.diploid.txt"
    if not os.path.exists(f):
        continue
    a = segments(f)
    a["genes"] = [",".join(sorted({g for gs, ge, g in genes[c] if gs < e and ge > st})) for c, st, e in zip(a["chrom"], a["start"], a["end"])]
    a.drop(columns=["variants", "DAV_variants"], errors="ignore").to_csv(f"{D}/segments.{s}.tsv", sep="\t", index=False)
    tt = totals(a)
    hits = a[a["genes"].apply(lambda x: any(g in x.split(",") for g in KNOWN))]
    tt["known_gene_hits"] = [dict(chrom=r.chrom, start=int(r.start), end=int(r.end), cls=r.cls,
                                  genes=[g for g in r.genes.split(",") if g in KNOWN]) for r in hits.itertuples()]
    tt["longest"] = [dict(chrom=r.chrom, start=int(r.start), end=int(r.end), kb=round(r.length / 1e3), cls=r.cls, genes=r.genes[:80])
                     for r in a.sort_values("length", ascending=False).head(10).itertuples()]
    if ctrl:
        c = pd.DataFrame(ctrl.values())
        tt["percentile_vs_controls"] = {k: round(float((c[k] < tt[k]).mean() * 100), 1) for k in ("total_Mb", "nea_Mb", "den_Mb")}
    out["samples"][s] = tt
json.dump(out, open(f"{D}/summary.json", "w"), ensure_ascii=False, indent=1)
print(json.dumps(out["samples"], ensure_ascii=False, indent=1)[:3000])
if ctrl:
    print(pd.DataFrame(ctrl.values()).groupby("pop")[["total_Mb", "nea_Mb", "den_Mb"]].median())
