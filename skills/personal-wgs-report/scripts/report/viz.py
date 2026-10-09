"""报告统一图表风格：viridis 色板 + 中文黑体 + 细线条 / 发丝网格。
配色规则（已用 dataviz 校验脚本验证，浅色背景 #fcfcfb）：
  - 分类（身份）：viridis 中段取样，最多 3 类 CAT3（任意两两 CVD ΔE≥19.9）；2 类用 CAT2（ΔE 40）
    绿色 #6ece58 对背景对比度不足 3:1 → 必须配直接标注或数据表
  - 连续数值：完整 viridis 色带 + 色标
  - 有序类别（年代等）：viridis 0–0.85 均匀取样 + 图例
  - “我”（样本）在所有图里都用 viridis 最深端 ME=#440154，并直接标注
"""
import glob
import os

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402

FIGS = "results/figs"   # make_figures.py 按样本改成 results/<样本>/figs


def _cjk_font():
    """找一个能显示中文的字体：tools/fonts 里的 Noto Sans SC → 系统里常见的中文字体。"""
    for f in sorted(glob.glob("tools/fonts/*.otf") + glob.glob("tools/fonts/*.ttf")):
        fm.fontManager.addfont(f)
    have = {f.name for f in fm.fontManager.ttflist}
    for name in ("Noto Sans SC", "Noto Sans CJK SC", "Source Han Sans SC", "Source Han Sans CN", "WenQuanYi Micro Hei",
                 "WenQuanYi Zen Hei", "Microsoft YaHei", "PingFang SC", "Hiragino Sans GB", "SimHei", "Arial Unicode MS"):
        if name in have:
            return name
    print("提示：没找到中文字体，图里的中文会显示成方块；运行 scripts/setup/install_tools.sh 会下载 Noto Sans SC")
    return "DejaVu Sans"


FONT = _cjk_font()

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
VIRIDIS = mpl.colormaps["viridis"]
CAT3 = ["#3e4c8a", "#1f978b", "#6ece58"]
CAT2 = ["#3e4c8a", "#6ece58"]
ME = "#440154"
NEUTRAL = "#d9d8d2"   # 背景/其他（非数据身份）


def ordinal(n, lo=0.0, hi=0.85):
    """有序类别：viridis 上均匀取 n 个色"""
    if n == 1:
        return [mpl.colors.to_hex(VIRIDIS(lo))]
    return [mpl.colors.to_hex(VIRIDIS(lo + (hi - lo) * i / (n - 1))) for i in range(n)]


mpl.rcParams.update({
    "font.family": FONT,
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.titleweight": "medium",
    "axes.titlelocation": "left",
    "axes.titlepad": 10,
    "axes.labelsize": 9,
    "axes.labelcolor": INK2,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.75,
    "axes.facecolor": SURFACE,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "grid.linestyle": "-",
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2,
    "xtick.major.size": 0,
    "ytick.major.size": 0,
    "xtick.major.pad": 4,
    "ytick.major.pad": 4,
    "figure.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "legend.frameon": False,
    "legend.fontsize": 8.5,
    "lines.linewidth": 1.5,
    "lines.solid_capstyle": "round",
    "lines.solid_joinstyle": "round",
    "svg.fonttype": "none",
    "axes.unicode_minus": False,
})


def header(ax, title, sub=None):
    """左对齐标题 + 下方一行说明（次级墨色），两者不重叠"""
    ax.set_title(title, pad=26 if sub else 10, loc="left")
    if sub:
        ax.annotate(sub, xy=(0, 1), xycoords="axes fraction", xytext=(0, 9), textcoords="offset points",
                    ha="left", va="bottom", fontsize=8.5, color=INK2)


def footer(axl, handles=None, text=None, ncol=2):
    """底部独立一行：图例 + 方法注（axl 是 GridSpec 里专门留的空轴）"""
    axl.axis("off")
    if handles:
        axl.legend(handles=handles, loc="upper left", ncol=ncol, fontsize=8, handlelength=1.2,
                   columnspacing=1.6, borderaxespad=0)
    if text:
        axl.text(0, 0, text, transform=axl.transAxes, ha="left", va="bottom", fontsize=7.5, color=MUTED)


def note(fig, text, x=0.01, y=0.005):
    """图底部的来源 / 方法注"""
    fig.text(x, y, text, ha="left", va="bottom", fontsize=7.5, color=MUTED)


def save(fig, name, table=None):
    """导出 SVG + PNG（200 dpi）；table（DataFrame）另存 CSV 作为数据表视图"""
    os.makedirs(FIGS, exist_ok=True)
    fig.savefig(f"{FIGS}/{name}.svg", bbox_inches="tight")
    fig.savefig(f"{FIGS}/{name}.png", dpi=200, bbox_inches="tight")
    if table is not None:
        table.to_csv(f"{FIGS}/{name}.csv", index=False)
    plt.close(fig)
    return f"{FIGS}/{name}.png"
