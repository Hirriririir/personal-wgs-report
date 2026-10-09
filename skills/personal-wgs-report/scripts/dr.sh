#!/usr/bin/env bash
# CPU 容器统一入口：按当前用户身份运行，路径与宿主机一致。用法：scripts/dr.sh <image> <cmd...>
# 需要从管道读标准输入时加 DR_STDIN=1（默认不接 stdin，免得在 while read 循环里把输入吃掉）
set -euo pipefail
. "$(dirname "$(readlink -f "$0")")/lib/mounts.sh"
stdin=(); [ "${DR_STDIN:-0}" = 1 ] && stdin=(-i)
exec docker run --rm "${stdin[@]}" -u "$(id -u):$(id -g)" -e HOME=/tmp "${MOUNTS[@]}" "$@"
