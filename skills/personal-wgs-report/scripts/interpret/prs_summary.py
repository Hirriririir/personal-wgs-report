#!/usr/bin/env python3
"""PRS 百分位：样本评分在 1000 Genomes 参考人群（config 里的 POP）分布中的位置 → work/prs/summary.json。
用 SCORE1_AVG（按实际有基因型的等位基因数归一），避免个别位点缺失把总分拉低。"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import csv
import json
import os

SEX_SPECIFIC = {"乳腺癌": "female", "前列腺癌": "male"}
out = []
for line in open(wgs.CONFIG.get("PGS_LIST") or "ref/pgs/selected.tsv"):
    if line.startswith("#") or not line.strip():
        continue
    trait, pid, note = line.rstrip("\n").split("\t")
    ref_p, cpl_p = f"work/prs/{pid}.ref.sscore", f"work/prs/{pid}.cohort.sscore"
    if not (os.path.exists(ref_p) and os.path.exists(cpl_p)):
        continue
    ref = [float(r["SCORE1_AVG"]) for r in csv.DictReader(open(ref_p), delimiter="\t")]
    cpl = {(r.get("#IID") or r.get("IID")): r for r in csv.DictReader(open(cpl_p), delimiter="\t")}
    n_used = sum(1 for _ in open(f"work/prs/{pid}.ref.sscore.vars")) if os.path.exists(f"work/prs/{pid}.ref.sscore.vars") else 0
    wp = f"work/prs/{pid}.weights.tsv"
    n_w = len({tuple(x.split(":")[:2]) for x in open(wp).read().split("\n")[1:] if x}) if os.path.exists(wp) else 0
    pct, z = {}, {}
    mu = sum(ref) / len(ref); sd = (sum((x - mu) ** 2 for x in ref) / (len(ref) - 1)) ** 0.5
    for s, r in cpl.items():
        v = float(r["SCORE1_AVG"])
        pct[s] = round(100 * sum(x < v for x in ref) / len(ref), 1)
        z[s] = round((v - mu) / sd, 2) if sd else None
    gl = open("work/prs/union.genotype.log").read().strip() if os.path.exists("work/prs/union.genotype.log") else ""
    out.append({"trait": trait, "pgs_id": pid, "note": note, "n_used": n_used, "n_weights": n_w, "pct": pct, "z": z,
                "ref_n": len(ref), "ref_pop": wgs.POP, "genotype_log": gl, "sex_specific": SEX_SPECIFIC.get(trait)})
json.dump(out, open("work/prs/summary.json", "w"), ensure_ascii=False, indent=1)
for r in out:
    print(r["trait"], r["pgs_id"], r["n_used"], r["pct"], r["z"], r["genotype_log"])
