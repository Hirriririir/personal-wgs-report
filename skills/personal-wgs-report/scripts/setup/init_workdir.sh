#!/usr/bin/env bash
# 新建工作目录：子目录、scripts 软链接（指向本仓库的 scripts）、config.env 模板。
# 用法：bash <仓库>/skills/personal-wgs-report/scripts/setup/init_workdir.sh <工作目录>
set -euo pipefail
W=${1:?用法：init_workdir.sh <工作目录>}
SCR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd -P)
mkdir -p "$W"/{raw,ref,work,qc,logs,results,tools,tmp}
ln -sfn "$SCR" "$W/scripts"
[ -f "$W/config.env" ] || cp "$SCR/../config.example.env" "$W/config.env"
cat <<MSG
工作目录已就绪：$W
下一步：
  1. 编辑 $W/config.env（样本名、FASTQ 路径、线程数、是否用 GPU）
  2. cd $W && bash scripts/setup/check_env.sh
MSG
