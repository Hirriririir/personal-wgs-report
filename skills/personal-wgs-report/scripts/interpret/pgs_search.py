#!/usr/bin/env python3
"""查 PGS Catalog：每个性状列出 GWAS 来源含目标人群（config 的 POP，默认 EAS）/ 在该人群验证过的评分，供挑选。
用法：pgs_search.py [性状名 ...]；结果挑好后写进 ref/pgs/selected.tsv（格式见 assets/pgs_selected_EAS.tsv）。"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import json
import sys
import urllib.request

TRAITS = {
    "CAD": "EFO_0001645", "T2D": "MONDO_0005148", "Breast cancer": "MONDO_0007254",
    "Prostate cancer": "MONDO_0008315", "Colorectal cancer": "MONDO_0005575",
    "Atrial fibrillation": "EFO_0000275", "Stroke": "EFO_0000712", "Hypertension": "EFO_0000537",
    "LDL cholesterol": "EFO_0004611", "BMI": "EFO_0004340", "Height": "OBA_VT0001253",
    "Gastric cancer": "MONDO_0001056", "Alzheimer": "MONDO_0004975", "Myopia": "HP_0000545",
    "Lung cancer": "MONDO_0008903", "Liver cancer (HCC)": "EFO_0000182", "Esophageal cancer": "MONDO_0007576",
    "Nasopharyngeal carcinoma": "MONDO_0015459", "Gout": "EFO_0004274",
}

def get(u):
    with urllib.request.urlopen(u, timeout=90) as r:
        return json.load(r)

only = sys.argv[1:] or list(TRAITS)
for name in only:
    efo = TRAITS[name]
    try:
        d = get(f"https://www.pgscatalog.org/rest/score/search?trait_id={efo}&include_children=0&limit=250")
    except Exception as e:  # noqa: BLE001
        print(name, "ERR", e)
        continue
    rows = []
    for s in d.get("results", []):
        dist = s.get("ancestry_distribution") or {}
        gwas = (dist.get("gwas") or {}).get("dist", {})
        ev = (dist.get("eval") or {}).get("dist", {})
        pub = s.get("publication") or {}
        rows.append((round(gwas.get(wgs.POP, 0), 1), round(ev.get(wgs.POP, 0), 1), s["id"], s["variants_number"],
                     pub.get("firstauthor", ""), (pub.get("date_publication") or "")[:4], s.get("name", "")))
    rows.sort(key=lambda r: (-(r[0] > 0), -r[1], -int(r[5] or 0)))
    print(f"== {name} ({efo}) n={d.get('count')}")
    for r in rows[:7]:
        print(f"   gwas{wgs.POP}%=" + "{} eval%={} {} nvar={} {} {} {}".format(*r))
