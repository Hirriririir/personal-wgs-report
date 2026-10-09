#!/usr/bin/env python3
"""下载小工具：断点续传、重试、能拿到校验值时核对 MD5。

用法：
  fetch.py url <URL> <目标文件>
  fetch.py zenodo <record_id> <目标目录> [文件名通配 ...]      # 例：fetch.py zenodo 13368126 ref/hmmix/archaic 'individuals_highcov.*'
  fetch.py dataverse <DOI> <目标目录> [文件名通配 ...]         # 例：fetch.py dataverse doi:10.7910/DVN/FFIDCW ref/aadr 'v66.p1_1240K*'
需要 curl；装了 aria2c 时大文件自动用多连接下载。代理沿用环境变量 https_proxy。
"""
import fnmatch
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request


def md5sum(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def get(url, dest, size=None, md5=None):
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    done = dest + ".ok"
    if os.path.exists(dest) and os.path.exists(done):
        print(f"已有 {dest}")
        return
    big = size and size > 500 * 1024 * 1024
    for attempt in range(1, 6):
        if big and shutil.which("aria2c"):
            cmd = ["aria2c", "-c", "-x", "8", "-s", "8", "--file-allocation=none", "--summary-interval=60",
                   "-d", os.path.dirname(dest) or ".", "-o", os.path.basename(dest), url]
        else:
            cmd = ["curl", "-fSL", "-C", "-", "--retry", "5", "--retry-delay", "10", "-o", dest, url]
        if subprocess.call(cmd) == 0:
            break
        print(f"下载失败，第 {attempt} 次重试：{url}", file=sys.stderr)
    else:
        sys.exit(f"放弃：{url}")
    if size and os.path.getsize(dest) != size:
        sys.exit(f"大小不对：{dest} {os.path.getsize(dest)} != {size}")
    if md5:
        got = md5sum(dest)
        if got != md5:
            sys.exit(f"MD5 不对：{dest} {got} != {md5}")
    open(done, "w").write(url + "\n")
    print(f"ok {dest}")


def api(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


def match(name, pats):
    return not pats or any(fnmatch.fnmatch(name, p) for p in pats)


def main():
    kind = sys.argv[1]
    if kind == "url":
        get(sys.argv[2], sys.argv[3])
    elif kind == "zenodo":
        rec, out, pats = sys.argv[2], sys.argv[3], sys.argv[4:]
        d = api(f"https://zenodo.org/api/records/{rec}")
        for f in d["files"]:
            if match(f["key"], pats):
                ck = f.get("checksum", "")
                get(f["links"]["self"], os.path.join(out, f["key"]), f.get("size"),
                    ck.split(":", 1)[1] if ck.startswith("md5:") else None)
    elif kind == "dataverse":
        doi, out, pats = sys.argv[2], sys.argv[3], sys.argv[4:]
        d = api(f"https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId={doi}")["data"]["latestVersion"]
        for f in d["files"]:
            df = f["dataFile"]
            if match(df["filename"], pats):
                ck = df.get("checksum") or {}
                md5 = ck.get("value") if ck.get("type") == "MD5" else df.get("md5")
                get(f"https://dataverse.harvard.edu/api/access/datafile/{df['id']}", os.path.join(out, df["filename"]),
                    df.get("filesize"), md5)
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
