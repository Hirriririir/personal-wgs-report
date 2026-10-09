# 给改这个仓库的 AI / 开发者

- 用户使用说明在 `skills/personal-wgs-report/SKILL.md`；这里是**维护仓库本身**的约定。
- 仓库里**不能有任何个人数据**：不要提交真实样本名、变异、结果文件、测序订单号、主机名 / IP / 路径、凭证。提交前跑一遍下面的检查。
- Shell 脚本：开头 `. "$(dirname "$(readlink -f "$0")")/lib/env.sh"`（子目录里是 `../lib/env.sh`），在工作目录里运行；容器一律经 `scripts/dr.sh`（CPU）或 `scripts/gpu_docker.sh`（GPU）；
  镜像版本只在 `scripts/lib/env.sh` 里定义。Python 脚本开头插入 `lib/wgs.py` 的引导（进入工作目录、读配置），路径一律相对工作目录。
- 注释和输出信息用中文，写清“为什么”而不是“做了什么”。
- 新增或改动分析步骤时同步更新：`scripts/run_all.sh`、`references/pipeline.md`、`scripts/report/collect_results.py`（底稿）、`assets/report_template.md`、必要时 `references/interpretation.md`。
- 人群相关参数放 `assets/population/<POP>.json`，不要写死在脚本里。

检查：

```bash
cd skills/personal-wgs-report
for f in $(find scripts -name '*.sh'); do bash -n "$f" || echo "语法错误 $f"; done
python3 -m py_compile $(find scripts -name '*.py') && find scripts -name __pycache__ -prune -exec rm -rf {} +
```
