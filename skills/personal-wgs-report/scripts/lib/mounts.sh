# 容器挂载参数（被 dr.sh / gpu_docker.sh 引用）：
#   工作目录按原路径挂载；工作目录顶层的软链接（scripts→仓库、ref→共享参考目录等）把真实目录也按原路径挂上；
#   配置文件（WGS_CONFIG，默认 config.env）里 FASTQ 所在的目录只读挂载。同一目录只挂一次。
W=$(pwd -P)
declare -A _seen=()
MOUNTS=(-v "$W:$W" -w "$W")
_add() {  # $1 目录 $2 ro/rw
  local d; d=$(readlink -f "$1") || return 0
  [ -d "$d" ] || return 0
  case "$d/" in "$W"/*) return 0;; esac
  [ -n "${_seen[$d]:-}" ] && return 0
  _seen[$d]=1
  MOUNTS+=(-v "$d:$d${2:+:$2}")
}
for e in "$W"/*; do [ -L "$e" ] && _add "$e" ""; done
_cfg=${WGS_CONFIG:-config.env}
if [ -f "$_cfg" ]; then
  for f in $( . "./$_cfg"; echo "${FASTQ_R1:-},${FASTQ_R2:-}" | tr ',' ' '); do
    _add "$(dirname "$f")" ro
  done
fi
