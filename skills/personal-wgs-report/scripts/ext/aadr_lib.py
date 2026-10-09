"""AADR TGENO（transpose_packed）读写小工具：每个个体一行，每位点 2 bit，值 = 第 5 列等位基因拷贝数，3 = 缺失。
文件头 48 字节（"TGENO nind nsnp ihash shash"），之后每行 ceil(nsnp/4) 字节；同一字节里第一个位点在高 2 位。"""
import numpy as np

HDR = 48


class TGeno:
    def __init__(self, prefix):
        self.prefix = prefix
        with open(prefix + ".geno", "rb") as f:
            head = f.read(HDR).split(b"\0")[0].split()
        assert head[0] == b"TGENO", head
        self.nind, self.nsnp = int(head[1]), int(head[2])
        self.rlen = (self.nsnp + 3) // 4
        self.mm = np.memmap(prefix + ".geno", dtype=np.uint8, mode="r", offset=HDR, shape=(self.nind, self.rlen))
        self.ind = [l.split() for l in open(prefix + ".ind")]
        assert len(self.ind) == self.nind

    def row(self, i):
        """解码第 i 个个体 → 长度 nsnp 的 uint8 数组"""
        b = np.asarray(self.mm[i])
        out = np.empty(self.rlen * 4, dtype=np.uint8)
        out[0::4] = (b >> 6) & 3
        out[1::4] = (b >> 4) & 3
        out[2::4] = (b >> 2) & 3
        out[3::4] = b & 3
        return out[: self.nsnp]


def _hashit(s):
    h = 0
    for ch in s.encode():
        h = (h * 23 + ch) & 0xFFFFFFFF
    return h


def calc_hash(ids):
    """ADMIXTOOLS 的 calcishash / calcsnphash（文件头里那两个十六进制数）"""
    h = 0
    for i in ids:
        h = ((h * 17) & 0xFFFFFFFF) ^ _hashit(i)
    return h


def pack(codes):
    """长度 nsnp 的 0/1/2/3 数组 → 一行 TGENO 字节"""
    n = len(codes)
    pad = (-n) % 4
    c = np.concatenate([codes.astype(np.uint8), np.full(pad, 3, np.uint8)])
    return ((c[0::4] << 6) | (c[1::4] << 4) | (c[2::4] << 2) | c[3::4]).astype(np.uint8)


def write_tgeno(prefix, rows, ind_lines, snp_src):
    """rows: 每个个体的 TGENO 字节行（np.uint8）；ind_lines: (id, sex, pop)；snp_src: 复制的 .snp 路径"""
    import shutil
    snp_ids = [l.split()[0] for l in open(snp_src)]
    nsnp = len(snp_ids)
    ih, sh = calc_hash([i for i, _, _ in ind_lines]), calc_hash(snp_ids)
    with open(prefix + ".geno", "wb") as f:
        head = f"TGENO {len(rows)} {nsnp} {ih:x} {sh:x}".encode()
        f.write(head + b"\0" * (HDR - len(head)))
        for r in rows:
            f.write(r.tobytes())
    with open(prefix + ".ind", "w") as f:
        for i, s, p in ind_lines:
            f.write(f"{i}\t{s}\t{p}\n")
    shutil.copyfile(snp_src, prefix + ".snp")
