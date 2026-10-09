#!/usr/bin/env python3
"""把 AI 写好的 Markdown 报告转成单个 HTML 文件（图片以 base64 内嵌，可离线打开、直接发给别人）。
用法：render_html.py results/<样本>/report.md [输出.html]（默认同名 .html）"""
import base64
import html
import os
import re
import sys

import markdown

src = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(src)[0] + ".html"
base = os.path.dirname(os.path.abspath(src))
text = open(src, encoding="utf-8").read()
text = re.sub(r"<!--.*?-->", "", text, flags=re.S)   # 去掉模板里的写作提示
body = markdown.markdown(text, extensions=["tables", "fenced_code", "toc", "sane_lists", "attr_list"])


def inline(m):
    path = m.group(2)
    p = path if os.path.isabs(path) else os.path.join(base, path)
    if path.startswith(("http://", "https://", "data:")) or not os.path.exists(p):
        return m.group(0)
    ext = os.path.splitext(p)[1].lower().lstrip(".")
    mime = {"svg": "image/svg+xml", "jpg": "image/jpeg", "jpeg": "image/jpeg"}.get(ext, f"image/{ext}")
    return f'{m.group(1)}data:{mime};base64,{base64.b64encode(open(p, "rb").read()).decode()}"'


body = re.sub(r'(<img[^>]*?src=")([^"]+)"', inline, body)
title = re.search(r"<h1[^>]*>(.*?)</h1>", body)
CSS = """
body{max-width:860px;margin:2.5rem auto;padding:0 1.2rem;font:16px/1.75 "Noto Sans SC","PingFang SC","Microsoft YaHei",sans-serif;color:#1b1b1a;background:#fcfcfb}
h1{font-size:1.9rem;line-height:1.3;margin:0 0 1rem}h2{font-size:1.35rem;margin:2.4rem 0 .8rem;padding-top:.6rem;border-top:1px solid #e1e0d9}
h3{font-size:1.1rem;margin:1.6rem 0 .5rem}img{max-width:100%;height:auto;display:block;margin:1rem auto}
table{border-collapse:collapse;width:100%;margin:1rem 0;font-size:.9rem;display:block;overflow-x:auto}
th,td{border-bottom:1px solid #e1e0d9;padding:.35rem .6rem;text-align:left;vertical-align:top}th{background:#f3f2ee}
blockquote{margin:1rem 0;padding:.6rem 1rem;background:#f3f2ee;border-left:4px solid #3e4c8a}
code{font-size:.88em;background:#f3f2ee;padding:.1em .3em;border-radius:3px}pre code{display:block;padding:.8rem;overflow-x:auto}
.muted,small{color:#6b6a65}
"""
doc = (f'<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
       f'<title>{html.escape(re.sub("<[^>]+>", "", title.group(1)) if title else "全基因组报告")}</title><style>{CSS}</style></head>'
       f'<body>{body}</body></html>')
open(out, "w", encoding="utf-8").write(doc)
print("写出", out, f"{len(doc) / 1e6:.1f} MB")
