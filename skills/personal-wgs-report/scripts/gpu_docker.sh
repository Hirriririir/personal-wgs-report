#!/usr/bin/env bash
# GPU 容器入口（Parabricks）。用法：scripts/gpu_docker.sh <image> <cmd...>
#   GPU_MODE=toolkit：机器装了 nvidia-container-toolkit，用 --gpus all（推荐）
#   GPU_MODE=manual ：没装 toolkit 时手动挂 /dev/nvidia* 设备和宿主机驱动库
#   GPU_MODE=auto   ：docker 里有 nvidia 运行时就用 toolkit，否则 manual
set -euo pipefail
. "$(dirname "$(readlink -f "$0")")/lib/mounts.sh"
mode=${GPU_MODE:-auto}
_cfg=${WGS_CONFIG:-config.env}
[ -f "$_cfg" ] && mode=$( . "./$_cfg"; echo "${GPU_MODE:-$mode}")
if [ "$mode" = auto ]; then
  docker info 2>/dev/null | grep -qi "runtimes:.*nvidia" && mode=toolkit || mode=manual
fi
args=(--rm -u "$(id -u):$(id -g)" -e HOME=/tmp)
if [ "$mode" = toolkit ]; then
  args+=(--gpus all)
else
  L=/usr/lib/x86_64-linux-gnu
  V=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -1)
  for dev in /dev/nvidia0 /dev/nvidiactl /dev/nvidia-uvm /dev/nvidia-uvm-tools /dev/nvidia-modeset; do
    [ -e "$dev" ] && args+=(--device "$dev")
  done
  for lib in libcuda.so libnvidia-ml.so libnvidia-ptxjitcompiler.so libnvidia-nvvm.so libcudadebugger.so; do
    [ -e "$L/$lib.$V" ] || continue
    case $lib in
      libcuda.so) args+=(-v "$L/$lib.$V:/usr/local/nvidia/lib64/libcuda.so.1:ro" -v "$L/$lib.$V:/usr/local/nvidia/lib64/libcuda.so:ro");;
      libnvidia-nvvm.so) args+=(-v "$L/$lib.$V:/usr/local/nvidia/lib64/$lib.4:ro");;
      *) args+=(-v "$L/$lib.$V:/usr/local/nvidia/lib64/$lib.1:ro");;
    esac
  done
  args+=(-v /usr/bin/nvidia-smi:/usr/bin/nvidia-smi:ro -e LD_LIBRARY_PATH=/usr/local/nvidia/lib64:/usr/local/cuda/lib64)
fi
exec docker run "${args[@]}" "${MOUNTS[@]}" "$@"
