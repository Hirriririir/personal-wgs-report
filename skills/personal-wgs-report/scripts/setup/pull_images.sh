#!/usr/bin/env bash
# 预先拉取全部容器镜像（约 25 GB），避免跑到一半才下载。不用 GPU 时跳过 Parabricks。
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
for v in $(compgen -v IMG_); do
  img=${!v}
  [ "$v" = IMG_PARABRICKS ] && ! use_gpu && continue
  [ "$v" = IMG_DEEPVARIANT ] && use_gpu && continue
  docker image inspect "$img" > /dev/null 2>&1 && { echo "已有 $img"; continue; }
  docker pull "$img" || echo "拉取失败：$img（国内可在 config.env 里把 BIOC 换成镜像站，或手动 docker pull 后 docker tag）"
done
echo "IMAGES DONE"
