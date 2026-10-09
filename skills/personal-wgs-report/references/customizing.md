# 定制：其他人群、家人、评分与面板

## 非东亚人群

核心流程（质控、比对、变异检测、注释、ACMG / 携带 / 药物基因组、STR、HLA、线粒体、Y）与人群无关，改 `POP=` 即可。和人群有关的地方：

1. **config.env**：`POP=EUR / AFR / AMR / SAS`（1000 Genomes 超级人群）。PCA 的第二套主成分、PRS 百分位、ROH / 古人类对照都会换成该人群。
2. **PRS 评分清单**：仓库只带了东亚的 `assets/pgs_selected_EAS.tsv`。其他人群：
   `tools/venv-annot/bin/python scripts/interpret/pgs_search.py "Coronary artery disease" "Type 2 diabetes" ...` 列出 GWAS 来源 / 验证人群含该人群的评分，
   挑 GRCh38 harmonized、位点数合理、在该人群验证过的，按同样三列格式（性状 TAB PGS 编号 TAB 说明）写成 `assets/pgs_selected_<POP>.tsv` 或直接写 `ref/pgs/selected.tsv`。
3. **扩展分析参数**：照 `assets/population/EAS.json` 写 `assets/population/<POP>.json`：
   - `lai.ref`：局部祖源的参考组（组名 → 1000 Genomes 人群代码列表），选该人群内部遗传上分得开的两三组；`controls` 选同超级人群的各人群；
   - `controls.pops`：古人类片段对照用的人群；
   - `aadr.countries` / `outgroups` / `modern_pops` / `map`：古 DNA 的地区、外群、定 PCA 轴的现代人群、地图范围；
   - `aadr.qpadm`：来源人群池、外群（right）、模型。这一部分需要群体遗传学背景：来源要是该人群历史上公认的祖先成分，外群要和来源有不同的亲缘关系。
     拿不准就删掉 `qpadm` 段（x02 会跳过 qpAdm），只做 PCA 和 f3。
   - `pop_names`：图里人群代码的中文名。
4. **性状位点**（`assets/genelists/trait_snps.tsv`）的说明文字有些是针对东亚写的频率，换人群时酌情改说明。

## 加入家人 / 伴侣

```bash
cp config.env config.partner.env    # 改 SAMPLE、FASTQ_R1/R2、SEX
WGS_CONFIG=config.partner.env bash scripts/run_all.sh   # 第二个人：单样本步骤 + 多人步骤（自动按两个人重算）
bash scripts/run_all.sh                                 # 再给第一个人跑一次：已完成的单样本步骤跳过，证据页和报告素材按两人结果更新
```

- 两个人时会多出 `work/findings/pair.*.tsv`（同一隐性基因都携带 / 一方致病另一方 VUS）。
- 亲子 / 同胞：somalier relate 会给出亲缘系数（亲子 / 同胞 ≈0.5），可用来核对样本没拿错。本流程没有做家系定相和新发变异检测。
- 每个人的报告单独写；涉及另一个人的结果（如伴侣的携带情况）只在两人都同意时写进对方的报告。

## 换 / 加基因面板

- ACMG 次要发现：`assets/genelists/acmg_sf_v3.3.tsv`（ACMG 出新版时更新基因和 report 规则）。
- 神经肌肉面板：PanelApp Australia 的面板编号在 `get_resources.sh`（3101 肌病、4092 神经肌肉）；想换成心血管、肿瘤等面板，改编号并相应修改 `build_gene_tables.py` 与 `findings.py` 里的 `nmd_*` 列用法。
- 性状位点：在 `assets/genelists/trait_snps.tsv` 加行（rsID、基因、主题、效应等位、说明）；坐标列留空时首次运行会查 Ensembl 补全。

## 更新数据库后重跑解读

ClinVar 每周更新、VEP / gnomAD 每年更新。只重跑解读部分：重新下载 ClinVar（`rm -f ref/annot/clinvar/clinvar*` 后跑 `get_resources.sh` 和 `prep_annotation.sh`；
旧文件一定要删掉，否则断点续传会把新旧两个版本拼在一起），
再删 `logs/done/10_annotate.cohort logs/done/10b_spliceai.cohort logs/done/11_findings.cohort`，运行 `run_all.sh`。
