#!/usr/bin/env python3
"""打印 / 导出 ExpansionHunter + stranger 结果：每个主位点一行。用法：str_table.py <sample>"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import csv
import sys

s = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
rows = []
for line in open(f"work/{s}/str/{s}.eh.stranger.vcf"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
    fmt = dict(zip(f[8].split(":"), f[9].split(":")))
    if info.get("VARID") != info.get("REPID"):
        continue
    rows.append({"locus": info.get("REPID"), "disease": info.get("Disease", ""), "moi": info.get("InheritanceMode", ""),
                 "ru": info.get("DisplayRU") or info.get("RU", ""), "repcn": fmt.get("REPCN", ""), "ci": fmt.get("REPCI", ""),
                 "normal_max": info.get("STR_NORMAL_MAX", ""), "patho_min": info.get("STR_PATHOLOGIC_MIN", ""),
                 "status": info.get("STR_STATUS", ""), "filter": f[6], "support": fmt.get("ADSP", "") + "|" + fmt.get("ADFL", "") + "|" + fmt.get("ADIR", ""),
                 "lc": fmt.get("LC", "")})
with open(f"work/{s}/str/{s}.str_summary.tsv", "w") as o:
    w = csv.DictWriter(o, list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow(r)
for r in rows:
    print(f'{r["locus"]:<12}{r["repcn"]:<10}{r["ci"]:<16}≤{r["normal_max"]:<5}≥{r["patho_min"]:<6}{r["status"]:<22}{r["filter"]:<10}{r["disease"][:40]}')
