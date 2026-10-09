#!/usr/bin/env python3
"""qpAdm：把样本（和对照人群）建模成几个古代来源人群的混合，给出比例和模型能否被拒绝（p 值）。
模型、来源人群池、外群（right，base / ext 两套）都在 assets/population/<POP>.json 的 aadr.qpadm 段。
先把用到的个体抽出来、来源人群换成池化标签，写成小 TGENO 加速；再逐个模型跑 ADMIXTOOLS qpAdm（allsnps: YES）。
  aadr_qpadm.py            样本 + 对照人群 × 全部模型 → work/ext/aadr/qpadm/results.tsv
  aadr_qpadm.py --indiv    个体级对照：从前 3 个对照人群各随机挑 6 人单独当目标（统计功效与单个样本相当），只跑主模型
                           → results_indiv.tsv（单人结果的误差天然比人群大，比较时要用这个）"""
import os
import random
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import popcfg  # noqa: E402
import wgs  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from aadr_lib import TGeno, write_tgeno  # noqa: E402

cfg = popcfg.load()["aadr"]["qpadm"]
D = "work/ext/aadr"
Q = f"{D}/qpadm"
os.makedirs(Q, exist_ok=True)
POOL = cfg["pools"]
RIGHT = {"base": cfg["right_base"], "ext": cfg["right_base"] + cfg.get("right_ext_add", [])}
SAMPLES = wgs.all_samples()
TARGETS = SAMPLES + cfg.get("compare_targets", [])
MODELS = [tuple(m) for m in cfg["models"]]
IMG = os.environ.get("IMG_ADMIXTOOLS")
INDIV = "--indiv" in sys.argv

if not os.path.exists(f"{Q}/sub.geno"):
    lab_of = {g: p for p, gs in POOL.items() for g in gs}
    need = set(lab_of) | set(RIGHT["ext"]) | set(TARGETS)
    t = TGeno(f"{D}/sub")
    rows, ind = [], []
    for i, (gid, sex, label) in enumerate(t.ind):
        if label in need:
            rows.append(np.asarray(t.mm[i]))
            ind.append((gid, sex, lab_of.get(label, label)))
    write_tgeno(f"{Q}/sub", rows, ind, f"{D}/sub.snp")
    print("qpAdm 子集个体", len(ind))

indfile = "sub.ind"
pick = {}
if INDIV:
    random.seed(42)
    ind = [x.split() for x in open(f"{Q}/sub.ind")]
    for p in cfg.get("compare_targets", [])[:3]:
        ids = sorted(i for i, _, lab in ind if lab == p)
        for i in random.sample(ids, min(6, len(ids))):
            pick[i] = p
    open(f"{Q}/sub.indiv.ind", "w").write("".join(f"{i}\t{sx}\t{('IND_' + i) if i in pick else lab}\n" for i, sx, lab in ind))
    indfile = "sub.indiv.ind"
    TARGETS = [f"IND_{i}" for i in pick] + SAMPLES
    MODELS = [tuple(cfg["main_model"])]


def run(target, model, rname, right):
    if target in model:
        return None
    right = [r for r in right if r not in model]
    tag = f"{'I__' if INDIV else ''}{target}__{'+'.join(model)}__{rname}"
    open(f"{Q}/{tag}.left", "w").write("\n".join((target,) + model) + "\n")
    open(f"{Q}/{tag}.right", "w").write("\n".join(right) + "\n")
    open(f"{Q}/{tag}.par", "w").write(
        f"genotypename: {Q}/sub.geno\nsnpname: {Q}/sub.snp\nindivname: {Q}/{indfile}\npopleft: {Q}/{tag}.left\n"
        f"popright: {Q}/{tag}.right\ndetails: NO\nallsnps: YES\ninbreed: NO\nhashcheck: NO\nnumchrom: 22\n")
    pr = subprocess.run(["scripts/dr.sh", IMG, "qpAdm", "-p", f"{Q}/{tag}.par"], capture_output=True, text=True)
    out = pr.stdout + ("\n#STDERR\n" + pr.stderr if pr.stderr else "")
    open(f"{Q}/{tag}.out", "w").write(out)
    k = len(model)
    coef = se = None
    p = float("nan")
    for line in out.splitlines():
        if line.startswith("best coefficients:"):
            coef = [float(x) for x in line.split(":")[1].split()]
        elif line.strip().startswith("std. errors:"):
            se = [float(x) for x in line.split(":")[1].split()]
        m = re.match(r"\s*f4rank:\s*(\d+)\s+dof:\s*(\d+)\s+chisq:\s*([\d.]+)\s+tail:\s*([\deE.+-]+)", line)
        if m and int(m.group(1)) == k - 1:
            p = float(m.group(4))
    if k == 1:
        coef, se = [1.0], [0.0]
    return dict(target=target, pop=pick.get(target.replace("IND_", ""), target), model="+".join(model), right=rname, p=p,
                **{f"w{i + 1}": (coef[i] if coef else float("nan")) for i in range(k)},
                **{f"se{i + 1}": (se[i] if se else float("nan")) for i in range(k)})


jobs = [(tg, m, rn, r) for tg in TARGETS for m in MODELS for rn, r in RIGHT.items()]
with ThreadPoolExecutor(int(os.environ.get("JOBS", "6"))) as ex:
    res = [r for r in ex.map(lambda a: run(*a), jobs) if r]
df = pd.DataFrame(res)
df.to_csv(f"{Q}/results{'_indiv' if INDIV else ''}.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(df.round(3).to_string(index=False))
