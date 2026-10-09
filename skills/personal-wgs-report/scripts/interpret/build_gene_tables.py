#!/usr/bin/env python3
"""汇总解读用的基因表 → ref/genelists/gene_table.tsv

每个基因一行：
  gencc_moi       GenCC（Definitive/Strong/Moderate）给出的遗传方式集合：AR / AD / XL / SD / MT
  gencc_diseases  对应疾病名（去重，截断）
  acmg_sf         ACMG SF v3.3 报告规则（all / biallelic / HFE / TTN），否则空
  nmd_myopathy    PanelApp AUS Myopathy Superpanel 置信度（3=green, 2=amber）；没下载到面板时为空
  nmd_superpanel  PanelApp AUS Neuromuscular Superpanel 置信度
  nmd_moi         PanelApp 给的遗传方式
  loeuf, pli      gnomAD v4.1 约束（MANE / canonical 转录本）
"""
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录

G = "ref/genelists"
genes = collections.defaultdict(lambda: {"moi": set(), "dis": [], "level": {}})
RANK = {"Definitive": 3, "Strong": 2, "Moderate": 1}
MOI = {"Autosomal recessive": "AR", "Autosomal dominant": "AD", "X-linked": "XL", "X-linked recessive": "XL",
       "Semidominant": "SD", "Mitochondrial": "MT"}
for r in csv.DictReader(open(f"{G}/gencc_submissions.tsv"), delimiter="\t"):
    if r["classification_title"] not in ("Definitive", "Strong", "Moderate"):
        continue
    g = r["gene_symbol"]
    m = MOI.get(r["moi_title"])
    if m:
        genes[g]["moi"].add(m)
        lv = RANK[r["classification_title"]]
        genes[g]["level"][m] = max(genes[g]["level"].get(m, 0), lv)
    d = r["disease_title"]
    if d and d not in genes[g]["dis"]:
        genes[g]["dis"].append(d)

acmg = {}
for line in open(f"{G}/acmg_sf_v3.3.tsv"):
    if line.startswith(("#", "gene\t")):
        continue
    f = line.rstrip("\n").split("\t")
    acmg[f[0]] = (f[3], f[1], f[2])

def panel(path):
    """超级面板由多个子面板合并，同一基因可出现多次：置信度取最高，遗传方式 / 表型合并。
    有 .tsv 用 .tsv；只有 PanelApp API 的 .json 时先转换；都没有返回空表（面板是可选的）。"""
    out = {}
    tsv = path
    js = path[:-4] + ".json"
    if not os.path.exists(tsv) and os.path.exists(js):
        d = json.load(open(js))
        with open(tsv, "w") as o:
            o.write("gene\tconfidence\tmoi\tphenotypes\n")
            for g in d.get("genes", []):
                sym = (g.get("gene_data") or {}).get("gene_symbol") or g.get("entity_name")
                phe = "; ".join(g.get("phenotypes") or []).replace("\t", " ")
                o.write(f"{sym}\t{g.get('confidence_level', '')}\t{g.get('mode_of_inheritance', '')}\t{phe}\n")
    if not os.path.exists(tsv):
        print(f"[可选] 没有 {tsv}，神经肌肉面板列留空")
        return out
    with open(tsv) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            g, conf, moi, phe = p[0], int(p[1] or 0), p[2], " ".join(p[3:])
            if g in out:
                c0, m0, ph0 = out[g]
                moi = m0 if moi in m0 else (m0 + " / " + moi if m0 else moi)
                phe = ph0 if phe in ph0 else (ph0 + "; " + phe)[:400]
                conf = max(conf, int(c0))
            out[g] = (str(conf), moi, phe)
    return out


myo = panel(f"{G}/panelapp_aus_3101.tsv")
nmd = panel(f"{G}/panelapp_aus_4092.tsv")

cons = {}
with open("ref/annot/gnomad/gnomad.v4.1.constraint_metrics.tsv") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        g = r.get("gene")
        if not g:
            continue
        mane = r.get("mane_select", "").lower() == "true"
        canon = r.get("canonical", "").lower() == "true"
        if g in cons and not mane:
            continue
        if mane or canon or g not in cons:
            cons[g] = (r.get("lof.oe_ci.upper", ""), r.get("lof.pLI", ""))

allg = sorted(set(genes) | set(acmg) | set(myo) | set(nmd))
with open(f"{G}/gene_table.tsv", "w") as o:
    o.write("gene\tgencc_moi\tmoi_strong\tgencc_diseases\tacmg_sf\tacmg_category\tacmg_condition\tnmd_myopathy\tnmd_superpanel\tnmd_moi\tnmd_phenotypes\tloeuf\tpli\n")
    for g in allg:
        a = acmg.get(g, ("", "", ""))
        m = myo.get(g, ("", "", ""))
        n = nmd.get(g, ("", "", ""))
        moi = ",".join(sorted(genes[g]["moi"])) if g in genes else ""
        dis = " | ".join(genes[g]["dis"][:6]) if g in genes else ""
        lo, pli = cons.get(g, ("", ""))
        strong = ",".join(sorted(k for k, v in (genes[g]["level"].items() if g in genes else []) if v >= 2))
        o.write("\t".join([g, moi, strong, dis, a[0], a[1], a[2], m[0], n[0], n[1] or m[1], (n[2] or m[2])[:300], lo, pli]) + "\n")
c = collections.Counter()
for g in allg:
    mo = genes[g]["moi"] if g in genes else set()
    if "AR" in mo: c["AR"] += 1
    if "XL" in mo: c["XL"] += 1
    if "AD" in mo: c["AD"] += 1
print(f"{len(allg)} genes; GenCC AR={c['AR']} XL={c['XL']} AD={c['AD']}; ACMG SF={len(acmg)}; myopathy panel={len(myo)}; NMD superpanel={len(nmd)}; constraint={len(cons)}")
