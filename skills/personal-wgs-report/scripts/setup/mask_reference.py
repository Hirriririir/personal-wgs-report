#!/usr/bin/env python3
"""把 GRC 官方认定的假重复 / 污染区（GRC_exclusions.bed）硬屏蔽为 N。
与 GIAB v3 参考的 maskedGRC 做法一致：假重复的第二份拷贝被屏蔽后，
CBS、U2AF1、CRYAA 等基因的读段不再被分到两处、MAPQ 不再被压成 0。
坐标不变，.fai 与原文件相同。"""
import mmap, shutil, sys
src, bed, dst = sys.argv[1:4]
fai = {}
for line in open(src + ".fai"):
    name, length, offset, lb, lbytes = line.split("\t")[:5]
    fai[name] = (int(length), int(offset), int(lb), int(lbytes))
shutil.copyfile(src, dst)
with open(dst, "r+b") as fh:
    mm = mmap.mmap(fh.fileno(), 0)
    masked = 0
    for line in open(bed):
        if line.startswith(("#", "browser", "track")):
            continue
        c, s, e, why = line.rstrip("\n").split("\t")[:4]
        if c not in fai:
            print(f"skip {c} (not in reference)")
            continue
        length, off, lb, lbytes = fai[c]
        s, e = int(s), min(int(e), length)
        pos = s
        while pos < e:
            row, col = divmod(pos, lb)
            n = min(lb - col, e - pos)
            b = off + row * lbytes + col
            seg = mm[b:b + n]
            masked += n - seg.count(b"N") - seg.count(b"n")
            mm[b:b + n] = b"N" * n
            pos += n
        print(f"masked {c}:{s}-{e} {why}")
    mm.flush(); mm.close()
print(f"newly masked bases: {masked}")
