#!/usr/bin/env bash
# 三个小分析（各自独立，一个失败不影响其他）：
#   bloodgroups.py   红细胞血型扩展分型（Duffy、Kidd、Kell、MNS、Diego、Rh CcEe、Lewis、分泌型、RhD-DEL 等；需联网查 Ensembl）
#   cnv_genes.py     拷贝数会变的基因：AMY1（唾液淀粉酶）、C4、LPA KIV-2、UGT2B17、GSTM1、CYP2A6、RHD 等
#   genome_stats.py  变异总数、杂合 / 纯合、Ti/Tv、1000 Genomes 里没有的“新”变异、杂合度沿染色体分布
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
set +e
s=$SAMPLE; mkdir -p work/ext/bloodgroup
$PYB scripts/ext/bloodgroups.py $s > work/ext/bloodgroup/$s.log 2>&1 || log "血型分型失败（见 work/ext/bloodgroup/$s.log）"
$PYB scripts/ext/cnv_genes.py $s > work/ext/cnv_genes.$s.log 2>&1 || log "基因拷贝数失败（见 work/ext/cnv_genes.$s.log）"
$PYB scripts/ext/genome_stats.py $s > work/ext/genome_stats.$s.log 2>&1 || log "基因组概览失败（见 work/ext/genome_stats.$s.log）"
log "MISC DONE $s"
