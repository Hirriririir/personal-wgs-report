#!/usr/bin/env bash
# 比对 + 标重复 + 变异检测。有 GPU 用 NVIDIA Parabricks（约 1 小时），没有就走 CPU（约 8–15 小时）：
#   GPU：pbrun fq2bam（BWA-MEM 比对、排序、标重复、QC 指标）→ pbrun deepvariant（WGS 模型，输出 gVCF）
#   CPU：bwa-mem2 → samtools fixmate / sort / markdup（光学重复距离同上）→ DeepVariant 官方 CPU 镜像
# 用法：scripts/02_align_call.sh [fq2bam|deepvariant ...]（默认两步都做；已有结果会跳过）
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
fastq_lists
s=$SAMPLE; steps=${*:-fq2bam deepvariant}
O=work/$s; T=tmp/$s; mkdir -p $O $T logs
RG="@RG\tID:$s.$LIBRARY\tLB:$LIBRARY\tPL:ILLUMINA\tPM:$PLATFORM_MODEL\tSM:$s\tPU:$s.$LIBRARY"
W=$(pwd -P)
for step in $steps; do case $step in
fq2bam)
  [ -s $O/$s.bam.bai ] && { log "$s.bam 已存在，跳过比对"; continue; }
  if use_gpu; then
    log "GPU 比对 $s"
    infq=(); for i in "${!FQ1[@]}"; do infq+=(--in-fq "$(readlink -f "${FQ1[$i]}")" "$(readlink -f "${FQ2[$i]}")" "${RG/ID:$s./ID:$s.$i.}"); done
    scripts/gpu_docker.sh $IMG_PARABRICKS pbrun fq2bam --ref $W/ref/GRCh38/bwa/GRCh38.fa "${infq[@]}" \
      --bwa-options="-K 100000000" --out-bam $W/$O/$s.bam --out-duplicate-metrics $W/$O/$s.dup_metrics.txt \
      --optical-duplicate-pixel-distance $OPTICAL_DUP_PIXEL_DISTANCE --out-qc-metrics-dir $W/$O/qc_metrics \
      --tmp-dir $W/$T --num-gpus 1 --memory-limit $(( MEM_GB * 3 / 4 )) --logfile $W/logs/$s.fq2bam.log
  else
    log "CPU 比对 $s（bwa-mem2，${THREADS} 线程）"
    cat_cmd() { for f in "$@"; do echo -n "$(readlink -f "$f") "; done; }
    r1="<(cat $(cat_cmd "${FQ1[@]}"))"; r2="<(cat $(cat_cmd "${FQ2[@]}"))"
    bash -c "$BIN/bwa-mem2 mem -t $THREADS -K 100000000 -R '$RG' $REF $r1 $r2" 2> logs/$s.bwamem2.log \
      | DR_STDIN=1 $D $IMG_SAMTOOLS sh -c "samtools fixmate -u -m - - | samtools sort -u -@ 8 -m 2G -T $T/sort - \
          | samtools markdup -@ 8 -d $OPTICAL_DUP_PIXEL_DISTANCE -f $O/$s.markdup_stats.txt --write-index - $O/$s.bam##idx##$O/$s.bam.bai"
  fi ;;
deepvariant)
  [ -s $O/$s.dv.g.vcf.gz ] && { log "$s.dv.g.vcf.gz 已存在，跳过变异检测"; continue; }
  if use_gpu; then
    log "GPU DeepVariant $s"
    scripts/gpu_docker.sh $IMG_PARABRICKS pbrun deepvariant --ref $W/ref/GRCh38/bwa/GRCh38.fa --in-bam $W/$O/$s.bam \
      --out-variants $W/$O/$s.dv.g.vcf.gz --gvcf --num-gpus 1 --tmp-dir $W/$T --logfile $W/logs/$s.deepvariant.log
    # Parabricks 加 --gvcf 时会同时写出 $s.dv.vcf.gz（只含变异位点）和 $s.dv.g.vcf.gz
    [ -s $O/$s.dv.vcf.gz ] || { echo "没找到 $O/$s.dv.vcf.gz（Parabricks 版本不同？）" >&2; exit 1; }
  else
    log "CPU DeepVariant $s（${THREADS} 个分片）"
    $D $IMG_DEEPVARIANT /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$REF --reads=$O/$s.bam \
      --output_vcf=$O/$s.dv.vcf.gz --output_gvcf=$O/$s.dv.g.vcf.gz --num_shards=$THREADS \
      --intermediate_results_dir=$T/dv > logs/$s.deepvariant.log 2>&1
  fi ;;
*) echo "未知步骤 $step"; exit 1;;
esac; done
log "ALIGN/CALL DONE $s"
