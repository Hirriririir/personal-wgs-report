"""扩展分析的人群参数：assets/population/<POP>.json（POP 来自 config）+ 1000 Genomes 样本表。"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402

ASSETS = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "assets")


def load():
    p = os.path.join(ASSETS, "population", f"{wgs.POP}.json")
    if not os.path.exists(p):
        sys.exit(f"没有 {p}：扩展分析的人群参数要先按 assets/population/EAS.json 的格式为 POP={wgs.POP} 写一份")
    return json.load(open(p))


def kg_samples():
    """1000 Genomes 3202 人：{IID: (超级人群, 人群)}，去掉二级以内亲属。"""
    rel = {x.split()[0] for x in open("ref/pca/deg2_hg38.king.cutoff.out.id") if not x.startswith("#")}
    out = {}
    for line in open("ref/pca/all_hg38.psam"):
        if line.startswith("#"):
            continue
        f = line.split()
        if f[0] not in rel:
            out[f[0]] = (f[4], f[5])
    return out
