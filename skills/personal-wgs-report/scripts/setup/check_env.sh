#!/usr/bin/env bash
# 开跑前的体检：配置、Docker、磁盘、内存、GPU、FASTQ。只检查，不改任何东西。
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
ok() { echo "  [ok]   $*"; }; warn() { echo "  [注意] $*"; }; bad() { echo "  [不行] $*"; FAIL=1; }
FAIL=0
echo "工作目录 $PWD，样本 $SAMPLE"
[[ "$SAMPLE" =~ ^[A-Za-z0-9_]+$ ]] && ok "样本名合法" || bad "样本名只能有字母、数字、下划线"
command -v docker > /dev/null && docker info > /dev/null 2>&1 && ok "Docker 可用（$(docker --version | cut -d, -f1)）" \
  || bad "Docker 不可用：装 Docker 并把当前用户加入 docker 组"
command -v python3 > /dev/null && ok "python3 $(python3 -V 2>&1 | cut -d' ' -f2)" || bad "没有 python3"
for t in curl pigz; do command -v $t > /dev/null && ok "$t" || warn "缺 $t（apt install $t）"; done
free_gb=$(df -BG --output=avail . | tail -1 | tr -dc 0-9)
[ "$free_gb" -ge 600 ] && ok "可用磁盘 ${free_gb} GB" || warn "可用磁盘只有 ${free_gb} GB；全流程峰值约 500–600 GB（参考与数据库约 250 GB + 每个样本 150–250 GB）"
mem_gb=$(awk '/MemTotal/{print int($2/1048576)}' /proc/meminfo)
[ "$mem_gb" -ge 60 ] && ok "内存 ${mem_gb} GB" || warn "内存 ${mem_gb} GB；建议 ≥64 GB（GLnexus、VEP、Parabricks 都吃内存）"
ok "线程 THREADS=$THREADS（本机 $(nproc) 个）"
if use_gpu; then
  ok "将用 GPU：$(nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader | head -1)"
  docker info 2>/dev/null | grep -qi "runtimes:.*nvidia" && ok "检测到 nvidia-container-toolkit" \
    || warn "没检测到 nvidia-container-toolkit，将用手动直通（GPU_MODE=manual）；装 toolkit 更省事"
else
  warn "不用 GPU：比对 + 变异检测走 CPU（bwa-mem2 + DeepVariant），30× 数据约 8–15 小时"
fi
fastq_lists
for f in "${FQ1[@]}" "${FQ2[@]}"; do
  [ -r "$f" ] && ok "FASTQ $(basename "$f")（$(du -h "$f" | cut -f1)）" || bad "读不到 FASTQ：$f"
done
[ "$FAIL" = 0 ] && echo "体检通过。下一步：bash scripts/setup/install_tools.sh" || { echo "有必须先解决的问题（见 [不行]）"; exit 1; }
