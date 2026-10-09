#!/usr/bin/env python3
"""按 mosdepth 汇总的 X / Y 相对常染色体深度推断核型性别。
XY：X≈0.5、Y≈0.5；XX：X≈1、Y≈0（chrY 的 PAR 在参考里已屏蔽，女性 Y 上只剩少量错配读段）。"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "lib"))
import wgs  # noqa: E402,F401  进入工作目录、读取 config.env
import sys
s = sys.argv[1]
d = {}
for line in open(f"work/{s}/qc/{s}.mosdepth.summary.txt"):
    f = line.split("\t")
    if f[0] != "chrom":
        d[f[0]] = float(f[3])
auto = sum(d[f"chr{i}"] for i in range(1, 23)) / 22
x, y = d["chrX"] / auto, d["chrY"] / auto
sex = "male" if (y > 0.2 and x < 0.75) else "female" if (y < 0.05 and x > 0.8) else "ambiguous"
out = f"{sex}\tautosomal_mean={auto:.2f}\tX_ratio={x:.3f}\tY_ratio={y:.3f}\n"
open(f"work/{s}/qc/{s}.sex.txt", "w").write(out)
print(s, out, end="")
