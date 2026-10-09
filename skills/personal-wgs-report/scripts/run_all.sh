#!/usr/bin/env bash
# 一键跑完整流程（在工作目录里执行；准备工作 setup/*.sh 要先做完，见 SKILL.md 第 2 步）。
# 每一步成功后写 logs/done/<步骤>.<样本>，重跑时自动跳过已完成的步骤；想重做某一步就删掉对应的 done 文件。
# 用法：
#   scripts/run_all.sh                 核心流程（质控 → 比对 / 变异检测 → 各项分析 → 解读 → 图 → 结果底稿）
#   scripts/run_all.sh --with-sv       加结构变异（Manta / Delly / CNVpytor，+3–5 小时，需 setup/get_resources_sv.sh）
#   scripts/run_all.sh --with-ext      加扩展分析（局部祖源、古 DNA、古人类片段、ROH、衰老、血型、病毒等，+半天，需 setup/get_resources_ext.sh）
#   scripts/run_all.sh --only 10_annotate   只跑某一步（名字见下面 STEPS）
#   第二个人：WGS_CONFIG=config.<名字>.env scripts/run_all.sh（多人步骤会把工作目录里所有人一起重算）
. "$(dirname "$(readlink -f "$0")")/lib/env.sh"
SV=0; EXT=0; ONLY=""
while [ $# -gt 0 ]; do case $1 in
  --with-sv) SV=1;; --with-ext) EXT=1;; --only) ONLY=$2; shift;;
  *) echo "未知参数 $1"; exit 1;;
esac; shift; done
mkdir -p logs/done
s=$SAMPLE
FAILED=()

# step <名字> <是否必须成功 1/0> <命令...>
step() {
  local name=$1 must=$2; shift 2
  [ -n "$ONLY" ] && [ "$ONLY" != "$name" ] && return 0
  local done=logs/done/$name.$s
  if [ -z "$ONLY" ] && [ -e "$done" ]; then log "跳过 $name（已完成）"; return 0; fi
  log "开始 $name"
  if "$@" > logs/$name.$s.log 2>&1; then
    touch "$done"; log "完成 $name"
  else
    log "失败 $name（见 logs/$name.$s.log）"
    FAILED+=("$name")
    [ "$must" = 1 ] && { echo "关键步骤失败，停止。修好后重跑 scripts/run_all.sh 会从这里继续。" >&2; exit 1; }
  fi
}
# 多人步骤：工作目录里所有人一起算；有人新加入后要重算
cohort_step() {
  local name=$1; shift
  local key; key=$(all_samples | tr '\n' '_')
  if [ -e logs/done/$name.cohort ] && [ "$(cat logs/done/$name.cohort)" != "$key" ]; then rm -f logs/done/$name.cohort; fi
  [ -n "$ONLY" ] && [ "$ONLY" != "$name" ] && return 0
  if [ -z "$ONLY" ] && [ -e logs/done/$name.cohort ]; then log "跳过 $name（已完成）"; return 0; fi
  log "开始 $name（$(all_samples | wc -l) 个样本）"
  if "$@" > logs/$name.cohort.log 2>&1; then echo "$key" > logs/done/$name.cohort; log "完成 $name"
  else log "失败 $name（见 logs/$name.cohort.log）"; FAILED+=("$name"); return 1; fi
}

# ---- 单样本：质控、比对、变异检测
step 01_qc            1 scripts/01_qc.sh
step 02_align_call    1 scripts/02_align_call.sh
step 03_bam_qc        1 scripts/03_bam_qc.sh
# ---- 单样本：各项专门分析（互相独立，失败不影响其他）
step 05_str           0 scripts/05_str.sh
step 06_special_loci  0 scripts/06_special_loci.sh
step 07_hla_kir       0 scripts/07_hla_kir.sh
step 08_mito          0 scripts/08_mito.sh
step 09_ychr          0 scripts/09_ychr.sh
# ---- 多人：联合分型 → 注释 → 解读
cohort_step 04_cohort_genotype scripts/04_cohort_genotype.sh || exit 1
cohort_step 10_annotate        scripts/10_annotate.sh || exit 1
cohort_step 10b_spliceai       scripts/10b_spliceai.sh || log "SpliceAI 失败：继续，只是少了剪接预测"
cohort_step 11_findings        scripts/11_findings.sh || exit 1
rm -f logs/done/11b_evidence.$s   # 依赖多人解读结果，很快，每次重做
step 11b_evidence     0 scripts/11b_evidence.sh
step 12_pgx           0 scripts/12_pgx.sh
cohort_step 13_ancestry        scripts/13_ancestry.sh || true
cohort_step 14_prs             scripts/14_prs.sh || true
# ---- 可选
if [ "$SV" = 1 ] || [ "$ONLY" = 15_sv ]; then step 15_sv 0 scripts/15_sv.sh; fi
if [ "$EXT" = 1 ] || [[ "$ONLY" == x* ]]; then
  cohort_step x00_panel scripts/ext/x00_panel.sh || true
  step x01_lai                 0 scripts/ext/x01_lai.sh
  cohort_step x02_aadr  scripts/ext/x02_aadr.sh || true
  cohort_step x03_roh   scripts/ext/x03_roh.sh || true
  cohort_step x04_archaic scripts/ext/x04_archaic.sh || true
  step x05_aging               0 scripts/ext/x05_aging.sh
  step x06_bloodgroup_cnv_stats 0 scripts/ext/x06_bloodgroup_cnv_stats.sh
  step x07_virome              0 scripts/ext/x07_virome.sh
fi
# ---- 报告素材（每次都重做，很快）
if [ -z "$ONLY" ] || [ "$ONLY" = report ]; then
  rm -f logs/done/report.$s
  step report 0 bash -c "$PYB scripts/report/make_figures.py $s && $PYB scripts/report/collect_results.py $s"
fi
if [ ${#FAILED[@]} -gt 0 ]; then
  log "完成，但这些步骤失败了：${FAILED[*]}（不影响其他结果；见 logs/<步骤>.*.log）"
else
  log "全部完成。结果底稿：results/$s/summary.md；图：results/$s/figs/；接下来按 SKILL.md 第 5 步写报告"
fi
