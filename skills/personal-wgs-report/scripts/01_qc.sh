#!/usr/bin/env bash
# FASTQ 质控：seqkit 读数统计 + fastp 质控报告（只出报告，不写过滤后的读段；比对软件自己会处理接头和低质量末端）。
# 测序公司给了 MD5 文件的话，先自己核对一遍（放在 FASTQ 同目录、名为 MD5.txt / md5.txt）。
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
fastq_lists
mkdir -p qc/fastp logs
for d in $(for f in "${FQ1[@]}" "${FQ2[@]}"; do dirname "$f"; done | sort -u); do
  for m in "$d"/MD5.txt "$d"/md5.txt; do
    [ -s "$m" ] && (cd "$d" && md5sum -c "$(basename "$m")") 2>&1 | tee -a logs/fastq_md5_check.txt
  done
done
$BIN/seqkit stats -a -T -j 4 "${FQ1[@]}" "${FQ2[@]}" > qc/$SAMPLE.seqkit_stats.tsv &
for i in "${!FQ1[@]}"; do
  $BIN/fastp -i "${FQ1[$i]}" -I "${FQ2[$i]}" -w 8 --detect_adapter_for_pe --overrepresentation_analysis \
    -j qc/fastp/$SAMPLE.$i.json -h qc/fastp/$SAMPLE.$i.html -R "$SAMPLE pair $i" 2> qc/fastp/$SAMPLE.$i.log &
done
# 读段名里的 仪器:运行号:芯片:泳道（Illumina）或 芯片+泳道（华大，形如 V350…L1）—— 数据来自几台仪器、几条泳道
( for f in "${FQ1[@]}"; do pigz -dc "$f" 2>/dev/null || gzip -dc "$f"; done \
  | awk 'NR%4==1{ n=split($1,a,":"); if (n>=4) k=a[1]":"a[2]":"a[3]":"a[4]; else if (match($1,/^@[A-Z]*[0-9]+L[0-9]+/)) k=substr($1,1,RLENGTH); else k="other"; c[k]++ }
         END{for (k in c) print c[k], k}' | sort -rn > qc/$SAMPLE.lanes.txt ) &
wait
echo "QC DONE $SAMPLE"
