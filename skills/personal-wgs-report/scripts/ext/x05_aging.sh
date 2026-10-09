#!/usr/bin/env bash
# 衰老相关指标：TelSeq 端粒长度（全 BAM 扫一遍，约 1 小时）+ aging.py（线粒体拷贝数、Y 丢失、克隆性造血筛查）
# 注意：TelSeq 的绝对值受建库、读长、测序平台影响很大，只适合同一批数据里比较；不要拿来和体检报告的“端粒年龄”对照
. "$(dirname "$(readlink -f "$0")")/../lib/env.sh"
s=$SAMPLE; O=work/ext/aging; mkdir -p $O
RL=$($PYB -c "import json,sys; q=json.load(open('work/findings/qc.json')).get('$s',{}); print(q.get('read_length') or 150)" 2>/dev/null || echo 150)
[ -s $O/telseq.$s.txt ] || $D $IMG_TELSEQ telseq -r $RL -m -o $O/telseq.$s.txt work/$s/$s.bam > $O/telseq.$s.log 2>&1 || log "TelSeq 失败"
$PYB scripts/ext/aging.py $s > $O/aging.$s.log 2>&1 && head -3 $O/aging.$s.log
log "AGING DONE $s"
