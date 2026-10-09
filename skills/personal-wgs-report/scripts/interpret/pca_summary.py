#!/usr/bin/env python3
"""把样本的 PCA 投影分数校准到参考 eigenvec 尺度（参考样本自身投影做线性回归），
再算到各人群中心的距离 → work/ancestry/pca_summary.json"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import csv
import json
import math
import statistics

def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))

psam = {r["#IID"]: r for r in tsv("ref/pca/all_hg38.psam")}
out = {}
for setn, label_col in (("kg_pruned", "SuperPop"), (f"kg_{wgs.POP.lower()}_pruned", "Population")):
    ev = {r["#IID"]: r for r in tsv(f"ref/pca/{setn}.pca.eigenvec")}
    sp = {r["#IID"]: r for r in tsv(f"ref/pca/{setn}.selfproj.sscore")}
    common = [i for i in ev if i in sp]
    fits = {}
    for k in ("PC1", "PC2", "PC3"):
        xs = [float(sp[i][f"{k}_AVG"]) for i in common]
        ys = [float(ev[i][k]) for i in common]
        mx, my = statistics.fmean(xs), statistics.fmean(ys)
        b = sum((a - mx) * (c - my) for a, c in zip(xs, ys)) / sum((a - mx) ** 2 for a in xs)
        r = statistics.correlation(xs, ys)
        fits[k] = (mx, my, b, r)
    proj = {r["#IID"]: r for r in tsv(f"work/ancestry/cohort.{setn}.proj.sscore")}
    cal = {s: {k: fits[k][1] + fits[k][2] * (float(proj[s][f"{k}_AVG"]) - fits[k][0]) for k in fits} for s in proj}
    groups = {}
    for i, r in ev.items():
        g = psam.get(i, {}).get(label_col, "?")
        groups.setdefault(g, []).append((float(r["PC1"]), float(r["PC2"]), float(r["PC3"])))
    cents = {g: tuple(statistics.fmean(p[j] for p in v) for j in range(3)) for g, v in groups.items()}
    sds = {g: tuple(statistics.pstdev(p[j] for p in v) for j in range(3)) for g, v in groups.items()}
    res = {"fit_r": {k: round(v[3], 4) for k, v in fits.items()}, "centroids": cents, "samples": {}}
    for s, c in cal.items():
        d = {g: math.dist((c["PC1"], c["PC2"]), cents[g][:2]) for g in cents}
        near = sorted(d, key=d.get)
        res["samples"][s] = {"PC": c, "nearest": near[:3], "dist": {g: round(d[g], 5) for g in near[:5]}}
    out[setn] = res
json.dump(out, open("work/ancestry/pca_summary.json", "w"), ensure_ascii=False, indent=1)
for setn, r in out.items():
    print(setn, "fit r:", r["fit_r"])
    for s, v in r["samples"].items():
        print(" ", s, {k: round(x, 4) for k, x in v["PC"].items()}, "nearest:", v["nearest"], v["dist"])
    print("  centroids:", {g: tuple(round(x, 4) for x in c[:2]) for g, c in r["centroids"].items()})
