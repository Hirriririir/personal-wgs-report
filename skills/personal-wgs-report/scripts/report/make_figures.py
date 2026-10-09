#!/usr/bin/env python3
"""把已有的分析结果画成报告用的图 → results/<样本>/figs/*.png（同名 .csv 是图背后的数据表）。
哪个分析没跑就跳过哪张图；单张图出错不影响其他图。用法：make_figures.py [sample] [图名 ...]

图名：pca_1kg lai aadr_pca f3 qpadm archaic roh prs cnv virome kir
"""
import collections
import glob
import gzip
import json
import os
import re
import sys
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "lib"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "ext"))
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import wgs  # noqa: E402
import viz  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from labels import cn  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

S = sys.argv[1] if len(sys.argv) > 1 else wgs.SAMPLE
ONLY = set(sys.argv[2:])
viz.FIGS = f"results/{S}/figs"
os.makedirs(viz.FIGS, exist_ok=True)
POPCFG = {}
try:
    import popcfg
    POPCFG = popcfg.load()
except SystemExit:
    pass
POPNAME = POPCFG.get("pop_names", {})
ME_LABEL = os.environ.get("FIG_PERSON", "你")   # 图里怎么称呼样本本人：报告用第二人称，写成博客时可设 FIG_PERSON=我


def popname(p):
    return POPNAME.get(p) or cn(p)


def kg_pop():
    out = {}
    for line in open("ref/pca/all_hg38.psam"):
        if not line.startswith("#"):
            f = line.split()
            out[f[0]] = (f[4], f[5])
    return out


def chrom_shapes():
    """染色体长度与着丝粒（UCSC cytoBand）"""
    clen, cen = {}, {}
    with gzip.open("ref/annot/cytoBand.hg38.txt.gz", "rt") as f:
        for line in f:
            c, s, e, band, stain = line.split("\t")
            if not re.fullmatch(r"chr\d+", c):
                continue
            clen[c] = max(clen.get(c, 0), int(e))
            if stain.strip() == "acen":
                lo, hi = cen.get(c, (1e12, 0))
                cen[c] = (min(lo, int(s)), max(hi, int(e)))
    return clen, cen


def karyogram_axes(n=22, extra=0.0):
    fig, (ax, axl) = plt.subplots(2, 1, figsize=(7.6, 0.34 * n + 2.0 + extra),
                                  gridspec_kw={"height_ratios": [0.34 * n + 0.6, 0.8], "hspace": 0.18})
    return fig, ax, axl


def finish_karyogram(ax, clen, n=22):
    L = max(clen.values())
    ax.set_xlim(-12, L / 1e6 + 2)
    ax.set_ylim(-n + 0.4, 0.75)
    ax.set_yticks([])
    ticks = list(range(0, int(L / 1e6) + 1, 50))
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{x} Mb" if x == ticks[-1] else f"{x}" for x in ticks])
    ax.grid(False)
    ax.spines["left"].set_visible(False)


def label_groups(ax, items, fontsize=7.4, min_dist=0.075):
    """在各组中心放标签，和已放的标签（按坐标轴比例）太近就跳过。items: [(x, y, 文本)]，按优先级排序"""
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    placed = []
    for x, y, text in items:
        fx, fy = (x - x0) / (x1 - x0), (y - y0) / (y1 - y0)
        if not (0 <= fx <= 1 and 0 <= fy <= 1) or any(abs(fx - a) < min_dist * 2.2 and abs(fy - b) < min_dist for a, b in placed):
            continue
        placed.append((fx, fy))
        if not text:
            continue
        ax.annotate(text, (x, y), fontsize=fontsize, color=viz.INK, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.2", fc=viz.SURFACE, ec="none", alpha=0.85), zorder=4)


def log_ticks(ax, lo, hi, unit=""):
    ticks = [t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000) if lo <= t <= hi]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" + (f" {unit}" if unit and t == ticks[-1] else "") for t in ticks])
    ax.minorticks_off()


def strip_plot(ax, groups, values_of, me_value, me_text, log=False):
    """对照人群每人一个点 + 中位数竖线，样本用菱形单独一行。groups: [(代码, 标签)]"""
    rng = np.random.default_rng(3)
    for i, (p, lab) in enumerate(groups):
        v = np.asarray(values_of(p))
        if len(v) == 0:
            continue
        ax.scatter(v, i + rng.uniform(-0.15, 0.15, len(v)), s=18, color=viz.CAT3[0], alpha=0.8,
                   edgecolors=viz.SURFACE, linewidths=0.8, zorder=3)
        ax.plot([np.median(v)] * 2, [i - 0.3, i + 0.3], color=viz.INK2, lw=1.2, zorder=2)
    yi = len(groups)
    ax.scatter([me_value], [yi], s=70, marker="D", color=viz.ME, edgecolors=viz.SURFACE, linewidths=1.5, zorder=4)
    ax.annotate(me_text, (me_value, yi), xytext=(8, 0), textcoords="offset points", va="center", fontsize=8.5, color=viz.INK)
    ax.set_yticks(range(len(groups) + 1))
    ax.set_yticklabels([lab for _, lab in groups] + [f"{ME_LABEL}（{S}）"])
    ax.grid(axis="y", visible=False)
    if log:
        ax.set_xscale("log")


# ---------------------------------------------------------------- 1000 Genomes PCA
def fig_pca_1kg():
    summ = json.load(open("work/ancestry/pca_summary.json"))
    psam = kg_pop()
    pop = wgs.POP.lower()
    sets = [("kg_pruned", 0, "全球 1000 Genomes（5 个大洲人群）"), (f"kg_{pop}_pruned", 1, f"{wgs.POP} 内部（1000 Genomes 各人群）")]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), gridspec_kw={"wspace": 0.18})
    tab = []
    for ax, (setn, col, title) in zip(axes, sets):
        if setn not in summ:
            ax.axis("off")
            continue
        ev = pd.read_csv(f"ref/pca/{setn}.pca.eigenvec", sep="\t")
        ev["grp"] = [psam.get(i, ("?", "?"))[col] for i in ev["#IID"]]
        groups = sorted(ev["grp"].unique())
        cols = viz.ordinal(len(groups), 0.05, 0.85)
        for g, c in zip(groups, cols):
            sub = ev[ev["grp"] == g]
            ax.scatter(sub["PC1"], sub["PC2"], s=7, color=c, alpha=0.55, edgecolors="none", zorder=1)
        me = summ[setn]["samples"].get(S)
        items = []
        if me:
            mx, my = me["PC"]["PC1"], me["PC"]["PC2"]
            ax.scatter([mx], [my], s=80, marker="D", color=viz.ME, edgecolors=viz.SURFACE, linewidths=1.5, zorder=5)
            ax.annotate(ME_LABEL, (mx, my), xytext=(8, -12), textcoords="offset points", fontsize=9, fontweight="bold", color=viz.INK, zorder=6)
            items.append((mx, my, ""))   # 占位：别的标签不要压在“我”上
            tab.append(dict(set=setn, PC1=mx, PC2=my, nearest=",".join(me["nearest"])))
        ax.margins(0.08)
        label_groups(ax, [(x, y, t) for x, y, t in items] +
                     [(ev.loc[ev["grp"] == g, "PC1"].median(), ev.loc[ev["grp"] == g, "PC2"].median(), popname(g) if col else g)
                      for g in groups], fontsize=7.6)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title(title, fontsize=10, loc="left")
    fig.text(0.125, 1.0, f"{ME_LABEL}在人群遗传版图上的位置（主成分分析）", fontsize=12, fontweight="medium", ha="left")
    fig.text(0.125, 0.955, f"每个小点 = 1000 Genomes 的一个人；菱形 = {ME_LABEL}（按参考面板的权重投影）；标签放在各人群的中心",
             fontsize=8.5, color=viz.INK2, ha="left")
    return viz.save(fig, "pca_1kg", pd.DataFrame(tab))


# ---------------------------------------------------------------- 局部祖源
def fig_lai():
    seg = pd.read_csv(f"work/ext/lai/segments.{S}.tsv", sep="\t")
    summ = json.load(open(f"work/ext/lai/summary.{S}.json"))
    groups = summ["groups"]                      # 编号(str) → 组名
    labels = POPCFG.get("lai", {}).get("labels", {})
    clen, cen = chrom_shapes()
    chroms = sorted(seg["chm"].unique(), key=lambda c: int(c[3:]))
    cols = viz.CAT2 if len(groups) == 2 else viz.ordinal(len(groups), 0.1, 0.8)
    COL = {int(k): cols[i] for i, k in enumerate(sorted(groups, key=int))}
    fig, ax, axl = karyogram_axes(len(chroms))
    H, GAP = 0.34, 0.08
    for r, c in enumerate(chroms):
        y0 = -r
        sub = seg[seg["chm"] == c]
        for h, colname in enumerate(("hap1", "hap2")):
            yb = y0 + (GAP / 2 if h == 0 else -GAP / 2 - H)
            ax.add_patch(Rectangle((0, yb), clen[c] / 1e6, H, color=viz.NEUTRAL, lw=0))
            for st, en, a in zip(sub["spos"], sub["epos"], sub[colname]):
                ax.add_patch(Rectangle((st / 1e6, yb), (en - st) / 1e6, H, color=COL[int(a)], lw=0))
        if c in cen:
            mid = sum(cen[c]) / 2e6
            ax.add_patch(Rectangle((mid - 0.6, y0 - GAP / 2 - H - 0.02), 1.2, 2 * H + GAP + 0.04, color=viz.SURFACE, lw=0))
        ax.text(-4, y0, c[3:], ha="right", va="center", fontsize=8.5, color=viz.INK2)
    finish_karyogram(ax, clen, len(chroms))
    fr = summ["fractions"]
    viz.header(ax, f"{ME_LABEL}的 22 对常染色体：每一段像哪个参考人群",
               "；".join(f"{labels.get(g, g).split('（')[0]} {fr[g]:.0%}" for g in fr) +
               "。每条染色体上下两行是两份拷贝（分别来自父母，但跨染色体不对应同一位亲本）")
    viz.footer(axl, [Patch(color=COL[int(k)], label=labels.get(g, g)) for k, g in sorted(groups.items(), key=lambda x: int(x[0]))] +
               [Patch(color=viz.NEUTRAL, label="无可用位点（着丝粒 / 端粒 / 重复区）")],
               "方法：参考人群 MAF≥1% 的 SNP；Beagle 5.5 定相（1000 Genomes 3,202 人参考）；RFMix v2；参考组留出部分个体作对照")
    out = [viz.save(fig, "lai_karyogram", seg)]
    # 对照
    first = sorted(groups.values())[0]
    ctrl = summ["controls"]
    order = sorted([(p, popname(p)) for p in ctrl], key=lambda x: float(np.median(ctrl[x[0]][first])))
    fig, ax = plt.subplots(figsize=(7.2, 0.45 * len(order) + 1.6))
    strip_plot(ax, order, lambda p: ctrl[p][first], fr[first], f"{ME_LABEL} {fr[first]:.0%}")
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.0%}"))
    ax.set_xlabel(f"基因组中被判为“{labels.get(first, first).split('（')[0]}”的比例")
    viz.header(ax, "同一把尺子量别人：留出的 1000 Genomes 个体", "圆点 = 一个人，竖线 = 中位数；这些人没有参与训练参考，用来看这把尺子准不准")
    out.append(viz.save(fig, "lai_calibration", pd.DataFrame([dict(pop=p, **{g: v for g, v in zip([first], [x])})
                                                               for p in ctrl for x in ctrl[p][first]])))
    return out


# ---------------------------------------------------------------- 古 DNA：PCA / f3 / qpAdm
def fig_aadr_pca():
    D = "work/ext/aadr"
    ev = pd.read_csv(f"{D}/pca/pca.evec", sep=r"\s+", skiprows=1, header=None, names=["id", "pc1", "pc2", "pc3", "pc4", "pop"])
    meta = pd.read_csv(f"{D}/sub.meta.tsv", sep="\t", low_memory=False)
    ev = ev.merge(meta[["gid", "bp", "country"]], left_on="id", right_on="gid", how="left")
    modern = [x.strip() for x in open(f"{D}/pca/modern.pops") if x.strip()]
    mod = ev[ev["pop"].isin(modern)]
    me = ev[ev["id"] == S].iloc[0]
    lo1, hi1 = mod["pc1"].quantile([0.005, 0.995])
    lo2, hi2 = mod["pc2"].quantile([0.005, 0.995])
    pad1, pad2 = (hi1 - lo1) * 0.12, (hi2 - lo2) * 0.12
    xlim, ylim = (lo1 - pad1, hi1 + pad1), (lo2 - pad2, hi2 + pad2)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 5.0), sharex=True, sharey=True, gridspec_kw={"wspace": 0.06})
    ax = axes[0]
    ax.scatter(mod["pc1"], mod["pc2"], s=9, color=viz.NEUTRAL, edgecolors="none", zorder=1)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    big = mod.groupby("pop").size().sort_values(ascending=False).index
    label_groups(ax, [(me["pc1"], me["pc2"], "")] + [(mod.loc[mod["pop"] == p, "pc1"].median(), mod.loc[mod["pop"] == p, "pc2"].median(), popname(p))
                                                     for p in big], fontsize=7)
    ax.set_title("现代人群（定义坐标轴）", fontsize=10, loc="left")
    ax = axes[1]
    ax.scatter(mod["pc1"], mod["pc2"], s=7, color="#ecebe6", edgecolors="none", zorder=0)
    anc = ev[(ev["bp"] > 0) & ev["pc1"].between(*xlim) & ev["pc2"].between(*ylim)].copy()
    bins = [0, 3000, 5000, 10000, 1e9]
    labs = ["3,000 年前以来", "5,000–3,000 年前", "1 万–5,000 年前", "1 万年以前"]
    cols = viz.ordinal(4, 0.0, 0.8)
    anc["per"] = pd.cut(anc["bp"], bins=bins, labels=labs, right=False)
    for lab, c in zip(labs[::-1], cols[::-1]):
        sub = anc[anc["per"] == lab]
        ax.scatter(sub["pc1"], sub["pc2"], s=11, color=c, edgecolors=viz.SURFACE, linewidths=0.4, zorder=2, label=f"{lab}（{len(sub)}）")
    ax.legend(loc="best", fontsize=7.5, title="古代个体（投影）", title_fontsize=7.5)
    ax.set_title("古代个体投影到同一坐标系", fontsize=10, loc="left")
    for ax in axes:
        ax.scatter([me["pc1"]], [me["pc2"]], s=90, marker="D", color=viz.ME, edgecolors=viz.SURFACE, linewidths=1.6, zorder=5)
        ax.annotate(ME_LABEL, (me["pc1"], me["pc2"]), xytext=(7, -2), textcoords="offset points", fontsize=9, color=viz.INK, fontweight="bold")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    fig.text(0.125, 1.03, f"{ME_LABEL}和古今人群在同一张遗传版图上（古 DNA 主成分分析）", ha="left", fontsize=12, fontweight="medium", color=viz.INK)
    fig.text(0.125, 0.985, f"数据 AADR v66.p1（1240K）；坐标轴只用现代人群计算，古人与{ME_LABEL}按最小二乘投影（smartpca lsqproject）",
             fontsize=8.5, color=viz.INK2, ha="left")
    tab = ev[ev["id"].eq(S) | ev["pop"].isin(modern) | ev["bp"].gt(0)][["id", "pop", "bp", "pc1", "pc2"]]
    return viz.save(fig, "aadr_pca", tab)


def fig_f3():
    f3 = pd.read_csv(f"work/ext/aadr/f3.{S}.tsv", sep="\t")
    anc = f3[(f3["bp"] > 0) & (f3["nsnp"] >= 30000)].copy()
    out = []
    top = anc.sort_values("f3", ascending=False).head(20).iloc[::-1]
    bins = [0, 3000, 5000, 10000, 1e9]
    labels = ["3,000 年前以来", "5,000–3,000 年前", "1 万–5,000 年前", "1 万年以前"]
    cols = viz.ordinal(4, 0.0, 0.8)
    top["period"] = pd.cut(top["bp"], bins=bins, labels=labels, right=False)
    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    y = np.arange(len(top))
    for yi, (_, r) in zip(y, top.iterrows()):
        c = cols[labels.index(r["period"])]
        ax.plot([r["f3"] - 1.96 * r["se"], r["f3"] + 1.96 * r["se"]], [yi, yi], color=viz.AXIS, lw=1.2, zorder=1)
        ax.scatter([r["f3"]], [yi], s=46, color=c, edgecolors=viz.SURFACE, linewidths=1.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{cn(g)}（约 {bp:,.0f} 年前，{n:.0f} 人）" for g, bp, n in zip(top["group"], top["bp"], top["n"])], fontsize=8)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel(f"outgroup f3（越大 = 与{ME_LABEL}共享的遗传历史越多）")
    viz.header(ax, f"与{ME_LABEL}共享遗传漂移最多的 20 个古代人群", "点 = 估计值，灰线 = 95% 置信区间；区间重叠的人群之间分不出高下")
    present = set(top["period"].astype(str))
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=c, markersize=7, label=lab) for c, lab in zip(cols, labels) if lab in present],
              loc="lower right", fontsize=8, title="年代", title_fontsize=8)
    out.append(viz.save(fig, "f3_top20", top.iloc[::-1]))
    # 地图（有 Natural Earth 底图时）
    if os.path.exists("ref/geo/ne_50m_land.geojson"):
        from matplotlib.collections import PolyCollection

        def polys(path):
            res = []
            for ft in json.load(open(path))["features"]:
                geo = ft["geometry"]
                parts = geo["coordinates"] if geo["type"] == "MultiPolygon" else [geo["coordinates"]]
                res += [np.array(p[0]) for p in parts]
            return res
        land = polys("ref/geo/ne_50m_land.geojson")
        lon = POPCFG.get("aadr", {}).get("map", {}).get("lon", [anc["lon"].min() - 5, anc["lon"].max() + 5])
        lat = POPCFG.get("aadr", {}).get("map", {}).get("lat", [anc["lat"].min() - 5, anc["lat"].max() + 5])
        m = anc[anc["lat"].notna() & anc["lon"].between(*lon) & anc["lat"].between(*lat)]
        vmin, vmax = m["f3"].quantile(0.1), m["f3"].max()
        periods = [(10000, 1e9, "1 万年以前"), (5000, 10000, "1 万–5,000 年前"), (3000, 5000, "5,000–3,000 年前"), (0.5, 3000, "3,000 年前以来")]
        fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.4), gridspec_kw={"hspace": 0.18, "wspace": 0.04})
        sc = None
        for ax, (lo, hi, title) in zip(axes.flat, periods):
            ax.add_collection(PolyCollection(land, facecolors="#e9e8e2", edgecolors="none", zorder=0))
            sub = m[(m["bp"] >= lo) & (m["bp"] < hi)].sort_values("f3")
            sc = ax.scatter(sub["lon"], sub["lat"], c=sub["f3"], cmap=viz.VIRIDIS.reversed(), vmin=vmin, vmax=vmax,
                            s=18 + 6 * np.sqrt(sub["n"]), edgecolors=viz.SURFACE, linewidths=0.9, zorder=3)
            for k, (_, r) in enumerate(sub.nlargest(2, "f3").iterrows()):
                ax.annotate(cn(r["group"]), (r["lon"], r["lat"]), xytext=(10, 9 if k == 0 else -15), textcoords="offset points",
                            fontsize=7.2, color=viz.INK, arrowprops=dict(arrowstyle="-", color=viz.MUTED, lw=0.6, shrinkA=0, shrinkB=3))
            ax.set_xlim(*lon)
            ax.set_ylim(*lat)
            ax.set_aspect(1.15)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.set_title(f"{title}（{len(sub)} 个人群）", fontsize=9.5, pad=4)
        cax = fig.add_axes([0.30, 0.045, 0.40, 0.014])
        cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=7.5, colors=viz.INK2, length=0)
        cb.set_label(f"与{ME_LABEL}共享的遗传漂移 f3 —— 颜色越深越相似", fontsize=8, color=viz.INK2)
        fig.suptitle(f"和{ME_LABEL}最像的古人住在哪里：{len(m)} 个古代人群", x=0.125, y=0.995, ha="left", fontsize=12, fontweight="medium")
        out.append(viz.save(fig, "f3_map", m.sort_values("f3", ascending=False)))
    return out


def fig_qpadm():
    q = POPCFG["aadr"]["qpadm"]
    r = pd.read_csv("work/ext/aadr/qpadm/results.tsv", sep="\t")
    model = "+".join(q["main_model"])
    sub = r[(r["model"] == model) & (r["right"] == "base")].sort_values("w1").copy()
    names = q.get("pool_names", {})
    s1, s2 = q["main_model"][:2]
    fig, ax = plt.subplots(figsize=(7.4, 0.42 * len(sub) + 1.6))
    y = np.arange(len(sub))
    H = 0.56
    for yi, (_, row) in zip(y, sub.iterrows()):
        w1 = row["w1"]
        if not 0 <= w1 <= 1:   # 比例跑出 0–100%：这个两源模型对该人群不适用
            ax.barh(yi, 1, height=H, color=viz.NEUTRAL, lw=0)
            ax.text(0.5, yi, "模型不适用（估计比例超出 0–100%）", ha="center", va="center", fontsize=7.8, color=viz.INK2)
            ax.text(1.02, yi, f"p={row['p']:.2f} 被拒绝", va="center", fontsize=8, color=viz.INK2)
            continue
        ax.barh(yi, w1, height=H, color=viz.CAT2[0], lw=0)
        ax.barh(yi, 1 - w1, left=w1, height=H, color=viz.CAT2[1], lw=0)
        ax.plot([w1, w1], [yi - H / 2, yi + H / 2], color=viz.SURFACE, lw=2)
        ax.errorbar(w1, yi, xerr=1.96 * row["se1"], fmt="none", ecolor=viz.INK, elinewidth=1, capsize=2.5, zorder=4)
        ok = "可接受" if row["p"] >= 0.05 else "被拒绝"
        ax.text(1.02, yi, f"{w1:.0%} / {1 - w1:.0%}   p={row['p']:.2f} {ok}", va="center", fontsize=8, color=viz.INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{ME_LABEL}（{t}）" if t == S else popname(t) for t in sub["target"]])
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.0%}"))
    ax.grid(axis="y", visible=False)
    ax.spines["left"].set_visible(False)
    viz.header(ax, "把基因组拆成两个古代来源：qpAdm 混合模型",
               "黑线 = 第一来源比例的 95% 置信区间；p ≥ 0.05 = 两源模型可接受，p < 0.05 = 被拒绝（还缺别的成分）")
    fig.subplots_adjust(right=0.78)
    ax.legend(handles=[Patch(color=viz.CAT2[0], label=names.get(s1, s1)), Patch(color=viz.CAT2[1], label=names.get(s2, s2))],
              loc="upper center", bbox_to_anchor=(0.45, -0.08), ncol=1, fontsize=8)
    return viz.save(fig, "qpadm", sub)


# ---------------------------------------------------------------- 古人类片段
def fig_archaic():
    a = pd.read_csv(f"work/ext/archaic/segments.{S}.tsv", sep="\t")
    summ = json.load(open("work/ext/archaic/summary.json"))
    me = summ["samples"][S]
    clen, cen = chrom_shapes()
    COL = {"Neanderthal": viz.CAT3[0], "Denisovan": viz.CAT3[2], "ambiguous": viz.MUTED, "unclassified": "#b9b8b1"}
    LAB = {"Neanderthal": "尼安德特人型", "Denisovan": "丹尼索瓦人型", "ambiguous": "两者难分", "unclassified": "古人类来源、但无法归类"}
    chroms = [f"chr{i}" for i in range(1, 23)]
    fig, ax, axl = karyogram_axes(22)
    H, MINW = 0.42, 0.35
    for r, c in enumerate(chroms):
        y = -r
        ax.add_patch(Rectangle((0, y - H / 2), clen[c] / 1e6, H, color="#e9e8e2", lw=0))
        if c in cen:
            mid = sum(cen[c]) / 2e6
            ax.add_patch(Rectangle((mid - 0.6, y - H / 2 - 0.02), 1.2, H + 0.04, color=viz.SURFACE, lw=0))
        sub = a[a["chrom"] == c]
        for cls in ("unclassified", "ambiguous", "Neanderthal", "Denisovan"):
            for s_, e in zip(sub.loc[sub["cls"] == cls, "start"], sub.loc[sub["cls"] == cls, "end"]):
                w = max((e - s_) / 1e6, MINW)
                ax.add_patch(Rectangle(((s_ + e) / 2e6 - w / 2, y - H / 2), w, H, color=COL[cls], lw=0))
        ax.text(-4, y, c[3:], ha="right", va="center", fontsize=8.5, color=viz.INK2)
    # 自动标注：著名渗入基因上的片段（每条染色体最多一个）
    done = set()
    for h in me.get("known_gene_hits", []):
        if h["chrom"] in done:
            continue
        done.add(h["chrom"])
        r = chroms.index(h["chrom"])
        x0 = (h["start"] + h["end"]) / 2e6
        ax.plot([x0], [-r + H / 2 + 0.11], marker="v", ms=4.2, color=viz.INK, ls="", zorder=5)
        ax.text(x0 + 2.0, -r + H / 2 + 0.03, f"{'/'.join(h['genes'])} · {LAB[h['cls']]}", fontsize=7, color=viz.INK, va="bottom", zorder=5)
    finish_karyogram(ax, clen)
    viz.header(ax, f"{ME_LABEL}身上的古人类“遗产”：尼安德特人与丹尼索瓦人片段",
               f"共 {me['n_segments']} 段、{me['total_Mb']:.0f} Mb（任一份拷贝）；尼安德特型 {me['nea_Mb']:.0f} Mb、丹尼索瓦型 {me['den_Mb']:.1f} Mb；"
               f"片段很短（中位 {me['median_kb']:.0f} kb），图中最短按 0.35 Mb 画；▼ = 文献中著名的渗入基因")
    viz.footer(axl, [Patch(color=COL[k], label=LAB[k]) for k in ("Neanderthal", "Denisovan", "ambiguous", "unclassified")],
               "方法：hmmix（Skov 2018），后验 ≥0.8；外群 = 1000G+HGDP 非洲人；古人类基因组 Altai、Vindija33.19、Chagyrskaya、Denisova", ncol=2)
    out = [viz.save(fig, "archaic_karyogram", a)]
    ctrl = pd.DataFrame(summ["controls"].values())
    if len(ctrl):
        order = [(p, popname(p)) for p in POPCFG.get("controls", {}).get("pops", sorted(ctrl["pop"].unique()))]
        fig, axes = plt.subplots(1, 2, figsize=(8.8, 0.42 * len(order) + 1.8), sharey=True, gridspec_kw={"wspace": 0.08})
        for ax, col, title in ((axes[0], "nea_Mb", "尼安德特型片段总长（Mb）"), (axes[1], "den_Mb", "丹尼索瓦型片段总长（Mb）")):
            strip_plot(ax, order, lambda p: ctrl.loc[ctrl["pop"] == p, col], me[col], f"{ME_LABEL} {me[col]:.1f}")
            ax.set_title(title, fontsize=9.5, loc="left")
            ax.margins(x=0.1)
        fig.text(0.125, 1.02, f"{ME_LABEL}的古人类成分算多还是少？和 1000 Genomes 对照比一比", fontsize=12, fontweight="medium", ha="left")
        out.append(viz.save(fig, "archaic_compare", ctrl))
    return out


# ---------------------------------------------------------------- ROH
def fig_roh():
    t = pd.read_csv("work/ext/roh/per_individual.tsv", sep="\t")
    z = t[t["sample"] == S].iloc[0]
    ref = t[t["pop"] != "sample"]
    order = [(p, popname(p)) for p in sorted(ref["pop"].unique())]
    fig, ax = plt.subplots(figsize=(7.2, 0.45 * len(order) + 1.8))
    vals = np.clip(ref["sum1"], 1, None)
    strip_plot(ax, order, lambda p: np.clip(ref.loc[ref["pop"] == p, "sum1"], 1, None), max(z["sum1"], 1),
               f"{ME_LABEL} {z['sum1']:.0f} Mb（{z['n1']} 段 ≥1 Mb，最长 {z['longest']:.1f} Mb）", log=True)
    hi = max(450, float(vals.max()) * 1.3, z["sum1"] * 1.3)
    lo = max(1, min(float(vals.min()), z["sum1"]) * 0.7)
    ax.set_xlim(lo, hi)
    log_ticks(ax, lo, hi, "Mb")
    ax.axvline(180, color=viz.AXIS, lw=1, zorder=1)
    ax.text(190, len(order) - 0.4, "一级表亲婚配的后代\n平均约 180 Mb（1/16）", fontsize=7.5, color=viz.INK2, va="center")
    ax.set_xlabel("≥1 Mb 纯合片段总长（Mb，对数刻度）")
    viz.header(ax, "父母有没有血缘关系？看“纯合片段”的总长度",
               "圆点 = 1000 Genomes 无亲缘个体，竖线 = 中位数；同一套位点、同一参数（bcftools roh）")
    return viz.save(fig, "roh", t)


# ---------------------------------------------------------------- PRS
def fig_prs():
    rows = json.load(open("work/prs/summary.json"))
    sex = wgs.sex_of(S)
    t = pd.DataFrame([dict(trait=r["trait"], pct=r["pct"][S], z=r["z"].get(S), pgs=r["pgs_id"], note=r["note"], n_used=r["n_used"])
                      for r in rows if S in r["pct"] and r.get("sex_specific") in (None, sex)])
    t = t.sort_values("pct", ascending=False).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.4, 0.33 * len(t) + 1.9))
    ys = -np.arange(len(t))
    ax.axvspan(20, 80, color="#efeee8", lw=0, zorder=0)
    ax.text(50, 0.7, "人群中间 60%", ha="center", va="bottom", fontsize=7.5, color=viz.MUTED)
    for yy, (_, r) in zip(ys, t.iterrows()):
        ax.plot([50, r["pct"]], [yy, yy], color=viz.AXIS, lw=1, zorder=1)
        ax.scatter([r["pct"]], [yy], s=48, color=viz.ME, edgecolors=viz.SURFACE, linewidths=1.5, zorder=3)
        ax.text(r["pct"] + (2.2 if r["pct"] < 93 else -2.2), yy, f"{r['pct']:.0f}", va="center",
                ha="left" if r["pct"] < 93 else "right", fontsize=7.8, color=viz.INK)
    ax.set_yticks(ys)
    ax.set_yticklabels(t["trait"], fontsize=8.5)
    ax.set_xlim(0, 100)
    ax.set_ylim(min(ys) - 0.8, 1.4)
    ax.set_xticks([0, 20, 50, 80, 100])
    ax.set_xticklabels(["0", "20", "50", "80", "100 百分位"])
    ax.grid(axis="y", visible=False)
    viz.header(ax, f"多基因风险评分：{ME_LABEL}在 1000 Genomes {wgs.POP} 人群中排第几",
               f"百分位 = 参考人群里分数比{ME_LABEL}低的比例；线从 50（中位）画到{ME_LABEL}的位置")
    viz.note(fig, "PRS 只反映“基因背景”的相对高低，不是患病概率；生活方式、家族史、体检指标通常比它重要得多。", y=0.0)
    return viz.save(fig, "prs", t)


# ---------------------------------------------------------------- 拷贝数会变的基因
def fig_cnv():
    d = json.load(open(f"work/ext/cnv_genes.{S}.json"))
    P = d["profiles"]
    genes = {}
    for line in open("ref/annot/gencode/pc_genes.bed"):
        c, a, b, g = line.split()
        genes.setdefault(g, (c, int(a), int(b)))
    PANELS = [("UGT2B17", ["UGT2B17", "UGT2B15"], f"UGT2B17：约 {d['UGT2B17']:.1f} 份", 4),
              ("GSTM1", ["GSTM1"], f"GSTM1：约 {d['GSTM1']:.1f} 份", 4),
              ("CYP2A6", ["CYP2A6", "CYP2A7"], f"CYP2A6：约 {d['CYP2A6']:.1f} 份", 4),
              ("AMY", ["AMY2B", "AMY2A", "AMY1A", "AMY1B", "AMY1C"], f"唾液淀粉酶 AMY1：合计约 {d['AMY1_total']:.0f} 份", 12),
              ("LPA", ["LPA"], f"LPA：KIV-2 重复合计约 {d['LPA_KIV2_total_copies']:.0f} 份", 30)]
    fig, axes = plt.subplots(len(PANELS), 1, figsize=(7.4, 1.55 * len(PANELS) + 0.8), gridspec_kw={"hspace": 1.25})
    for ax, (key, gl, title, ymax) in zip(axes, PANELS):
        p = P[key]
        x = (p["start"] + np.arange(len(p["cn"])) * p["win"] + p["win"] / 2) / 1e6
        y = np.array(p["cn"])
        ymax = max(ymax, float(np.nanpercentile(y, 99)) * 1.15)
        ax.fill_between(x, 0, y, step="mid", color=viz.CAT3[0], alpha=0.12, lw=0)
        ax.step(x, y, where="mid", color=viz.CAT3[0], lw=1.4)
        ax.axhline(2, color=viz.AXIS, lw=0.9, zorder=0)
        ax.text(x[-1], 2, " 2 份", va="center", ha="left", fontsize=7.5, color=viz.MUTED)
        for g in gl:
            if g not in genes:
                continue
            c, a, b = genes[g]
            ax.add_patch(plt.Rectangle((a / 1e6, -0.16 * ymax), (b - a) / 1e6, 0.07 * ymax, color=viz.INK2, lw=0, clip_on=False))
            ax.text((a + b) / 2e6, -0.2 * ymax, g, ha="center", va="top", fontsize=7, color=viz.INK2)
        ax.set_ylim(0, ymax)
        ax.set_xlim(x[0], x[-1])
        ax.set_title(title, fontsize=9.5, loc="left", pad=4)
        ax.set_ylabel("拷贝数", fontsize=8)
        ax.tick_params(labelsize=7.5)
        ax.tick_params(axis="x", pad=27)
        ax.set_xlabel(f"{p['chrom']} 位置（Mb）", fontsize=7.5, labelpad=2)
    fig.text(0.125, 0.995, f"拷贝数会“变”的基因：数一数{ME_LABEL}有几份", fontsize=12, fontweight="medium", ha="left")
    fig.text(0.125, 0.972, "读段深度折算成拷贝数（1 kb 窗口；灰线 = 正常的 2 份）；多拷贝基因用全部读段，其余只用唯一比对读段",
             fontsize=8.3, color=viz.INK2, ha="left")
    return viz.save(fig, "cnv_profiles")


# ---------------------------------------------------------------- 未比对读段
def fig_virome():
    rep = pd.read_csv(f"work/ext/virome/{S}.k2.report", sep="\t", header=None, names=["pct", "clade", "direct", "rank", "taxid", "name"])
    rep["name"] = rep["name"].str.strip()

    def get(n):
        return int(rep.loc[rep["name"] == n, "clade"].iloc[0]) if (rep["name"] == n).any() else 0
    total = get("unclassified") + get("root")
    dom, doms = "", []
    for _, r in rep.iterrows():   # 报告按分类树顺序排列：记下每行所属的域（D 级）
        if r["rank"] == "D" or (r["rank"] == "R1" and r["name"] == "Viruses"):   # 新版 NCBI 分类里病毒是 R1（非细胞根）
            dom = r["name"]
        elif r["rank"] in ("U", "R", "R1"):
            dom = ""
        doms.append(dom)
    rep["domain"] = doms
    genera = rep[rep["rank"] == "G"].sort_values("clade", ascending=False)
    bact_g = genera[genera["domain"] == "Bacteria"].head(2)
    rows = [("未能归类", get("unclassified"))]
    rest_b = get("Bacteria")
    for _, r in bact_g.iterrows():
        rows.append((f"细菌 · {r['name']}", int(r["clade"])))
        rest_b -= int(r["clade"])
    rows += [("细菌 · 其他", rest_b), ("人类序列（参考基因组未收录的版本）", get("Homo sapiens")), ("古菌", get("Archaea"))]
    vir = rep[(rep["rank"] == "F") & (rep["domain"] == "Viruses")].sort_values("clade", ascending=False)
    vtot = get("Viruses")
    for _, r in vir.head(4).iterrows():
        rows.append((f"病毒 · {r['name']}", int(r["clade"])))
        vtot -= int(r["clade"])
    if vtot > 0:
        rows.append(("病毒 · 其他", vtot))
    t = pd.DataFrame([r for r in rows if r[1] > 0], columns=["类别", "读段数"])
    fig, ax = plt.subplots(figsize=(7.4, 0.36 * len(t) + 1.4))
    y = list(range(len(t)))[::-1]
    ax.barh(y, t["读段数"], height=0.55, color=viz.CAT3[0], lw=0)
    for yi, v in zip(y, t["读段数"]):
        ax.text(v * 1.15, yi, f"{v:,}", va="center", fontsize=8, color=viz.INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(t["类别"], fontsize=8.3)
    ax.set_xscale("log")
    hi = max(10, t["读段数"].max()) * 8
    ax.set_xlim(1, hi)
    ticks = [10 ** k for k in range(0, 9) if 10 ** k <= hi]
    ax.set_xticks(ticks)
    ax.set_xticklabels([{1: "1", 10: "10", 100: "100", 1000: "1千", 10 ** 4: "1万", 10 ** 5: "10万", 10 ** 6: "100万", 10 ** 7: "1000万",
                         10 ** 8: "1亿"}[x] for x in ticks])
    ax.minorticks_off()
    ax.grid(axis="y", visible=False)
    viz.header(ax, "血液 DNA 里“不属于人类参考基因组”的读段", f"共 {total:,} 条读段，用 Kraken2 逐条分类；横轴为对数刻度")
    return viz.save(fig, "virome", t)


# ---------------------------------------------------------------- KIR
def fig_kir():
    calls = {}
    for line in open(f"work/{S}/kir/{S}_genotype.tsv"):
        f = line.rstrip("\n").split("\t")
        q1 = float(f[4]) if f[4] not in ("-1", ".") else -1
        calls[f[0]] = (f[2] if q1 > 0 else ".", f[5] if len(f) > 7 and f[7] not in (".", "-1") and float(f[7]) > 0 else ".")
    ORDER = [("KIR3DL3", "框架"), ("KIR2DS2", "B"), ("KIR2DL2", "B"), ("KIR2DL3", "A"), ("KIR2DL5B", "B"), ("KIR2DS3", "B"),
             ("KIR2DP1", "假基因"), ("KIR2DL1", "A"), ("KIR3DP1", "框架"), ("KIR2DL4", "框架"), ("KIR3DL1", "A"), ("KIR3DS1", "B"),
             ("KIR2DL5A", "B"), ("KIR2DS5", "B"), ("KIR2DS1", "B"), ("KIR2DS4", "A"), ("KIR3DL2", "框架")]
    COL = {"框架": viz.CAT3[1], "A": viz.CAT3[0], "B": viz.CAT3[2], "假基因": viz.MUTED}
    fig, ax = plt.subplots(figsize=(9.4, 2.4))
    W, H, GAP = 0.9, 0.9, 0.1
    rows = []
    for i, (g, cat) in enumerate(ORDER):
        a1, a2 = calls.get(g, (".", "."))
        present = a1 != "."
        rows.append(dict(gene=g, category=cat, present=present, alleles=",".join(a for a in (a1, a2) if a != ".")))
        x = i * (W + GAP)
        ax.add_patch(Rectangle((x, 0), W, H, fc=COL[cat] if present else viz.SURFACE, ec=COL[cat], lw=1.2))
        ax.text(x + W / 2, H + 0.12, g.replace("KIR", ""), ha="center", va="bottom", fontsize=8, color=viz.INK)
        if present:
            al = sorted({a.split("*")[1] for a in (a1, a2) if a != "." and "*" in a})
            ax.text(x + W / 2, H / 2, "\n".join("*" + a for a in al), ha="center", va="center", fontsize=6.6,
                    color="white" if cat != "B" else viz.INK)
        else:
            ax.text(x + W / 2, H / 2, "无", ha="center", va="center", fontsize=8, color=viz.MUTED)
    ax.set_xlim(-0.2, len(ORDER) * (W + GAP))
    ax.set_ylim(-0.3, 1.45)
    ax.axis("off")
    fig.text(0.125, 1.03, f"{ME_LABEL}的 KIR 基因（NK 细胞读取 HLA 的“读卡器”）", fontsize=12, fontweight="medium", ha="left")
    fig.text(0.125, 0.965, "实心 = 有该基因（方块内为等位基因），空心 = 缺失；按染色体上的排列顺序。A 单倍型基因全在、B 单倍型基因全无 = AA 型",
             fontsize=8.5, color=viz.INK2, ha="left")
    ax.legend(handles=[Patch(fc=COL[k], label=lab) for k, lab in (("框架", "框架基因（人人都有）"), ("A", "A 单倍型基因"),
                                                                   ("B", "B 单倍型基因"), ("假基因", "假基因"))],
              loc="upper right", bbox_to_anchor=(1.0, -0.02), ncol=4, fontsize=7.8)
    return viz.save(fig, "kir", pd.DataFrame(rows))


FIGS = {"pca_1kg": fig_pca_1kg, "lai": fig_lai, "aadr_pca": fig_aadr_pca, "f3": fig_f3, "qpadm": fig_qpadm, "archaic": fig_archaic,
        "roh": fig_roh, "prs": fig_prs, "cnv": fig_cnv, "virome": fig_virome, "kir": fig_kir}
made = collections.OrderedDict()
for name, fn in FIGS.items():
    if ONLY and name not in ONLY:
        continue
    try:
        r = fn()
        made[name] = r if isinstance(r, list) else [r]
        print("ok  ", name, made[name])
    except FileNotFoundError as e:
        print("跳过", name, f"（缺输入：{e.filename}）")
    except Exception:  # noqa: BLE001
        print("失败", name)
        traceback.print_exc()
idx_path = f"results/{S}/figs/index.json"
idx = json.load(open(idx_path)) if ONLY and os.path.exists(idx_path) else {}
idx.update({k: [os.path.relpath(p, f"results/{S}") for p in v] for k, v in made.items()})
json.dump(idx, open(idx_path, "w"), ensure_ascii=False, indent=1)
