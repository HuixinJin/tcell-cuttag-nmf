# T 细胞 Cut&Tag 数据的 NMF 聚类分析

对 4 种组蛋白修饰（H3K27ac、H3K4me1、H3K4me3、H3K9me3）的 T 细胞 Cut&Tag 数据
进行 NMF（非负矩阵分解）聚类，并对得到的 peak 簇进行 motif 富集、GO 富集和
信号可视化分析。

分析实际运行在集群上，工作目录为 `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/`。
本仓库按分析步骤整理代码，下方每个文件均标注其在服务器上的绝对路径，方便查找。

## 目录结构

```
├── 01_nmf/                  # NMF 聚类主分析
│   ├── t_nimfa.py           #   ★ 核心：rank=3 NMF + consensus 聚类
│   ├── t_nimfa_estimate.py  #   rank 估计（rank 2–6，五指标评估）
│   ├── cluster_peaks.py     #   簇分配 → 每簇 BED 文件（下游分析输入）
│   ├── debug_H3K9me3_nimfa.py  # H3K9me3 consensus 热图（大矩阵低内存版）
│   ├── exploration/         #   nimfa 探索性分析 notebook
│   └── run/                 #   SLURM 集群提交脚本
├── 02_deeptools/            # 各簇 peak 在 16 种 T 细胞状态的信号热图
│   ├── t_deeptools.sh       #   computeMatrix + plotHeatmap（通用版）
│   ├── t_deeptools-H3K9me3-rank3.sh  # H3K9me3 rank3 版
│   ├── t_plotheatmap.sh     #   单独重画热图
│   └── run/                 #   SLURM 提交脚本
├── 03_motif_enrichment/     # motif 富集分析（MEME-ChIP / HOMER / cisTarget）
└── 04_go_and_annotation/    # GO 富集与基因注释
```

## 各步骤说明

### 1. NMF 聚类（01_nmf）

- **输入矩阵**：`/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/matrix/<mark>_variable_peak_matrix.tsv`
  （行为 peak：`chr:start-end`；列为样本，列名末尾 5 字符为重复编号，脚本自动去掉后按样本取均值）
- **rank 估计**：服务器 `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/t_nimfa_estimate.py`
  （rank 2–6，评估 RSS / cophenetic / dispersion / sparsity / evar 五项指标）
- **NMF 主分析（核心脚本）**：服务器 `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/t_nimfa.py`
  - 预处理：`sklearn.preprocessing.Normalizer` 对每个 peak 的信号向量做 L2 归一化
  - 参数：`nimfa.Nmf(data.T, rank=3, seed='nndsvd', max_iter=100,
    update='divergence', objective='div', n_run=50, track_factor=True)`
  - 聚类：consensus 矩阵取 `1 - consensus()` 距离 → average linkage 层次聚类 →
    `fcluster(t=0.5, criterion='distance')` 切簇
  - 输出目录：`/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/t_rank3/result/<mark>/`
    （`cluster_assignments.txt`、`C.pkl`、`Y_linkage.pkl`）
- **簇分配转 BED**：服务器 `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/t_rank3/cluster_peaks.py`
  （输出每簇 `cluster_bed/cluster_<N>_peaks.bed`，是所有下游分析的输入）
- H3K9me3 consensus 热图（大矩阵低内存渲染）：服务器
  `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/debug_H3K9me3_nimfa.py`
- 探索性分析 notebook：服务器 `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/nimfa.ipynb`
- 集群提交脚本（服务器 `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/`）：
  `sbatch_t_consensus.sh`、`do_python.sh`、`sbatch_t_debug_H3K9me3.sh`、`debug_python.sh`

### 2. deeptools 可视化（02_deeptools）

- 脚本（服务器 `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/`）：
  `t_deeptools.sh`、`t_deeptools-H3K9me3-rank3.sh`、`t_plotheatmap.sh`，
  对应 SLURM 提交脚本 `sbatch_t_deeptools.sh`、`sbatch_t_deeptools-H3K9me3-rank3.sh`、
  `sbatch_t_computematrix.sh`
- 对每个簇的 BED 文件：`computeMatrix reference-point --referencePoint center -b 2000 -a 2000`
- 信号 bigWig 目录：`/dshare/xielab/analysisdata/ThreeDGene/zhangjk/TCellProject/CUT_TAG/pool/`
- 16 种 T 细胞状态：Naive、Arm-MP、Arm-Tcm、Arm-Trm、Arm-Tem、Arm-Teff、
  TdLN-Ttsm、TdLN-Tpex、TME-Tpex、TME-Ttran、TME-Tex、
  Cl-Tpex(d7)、Cl-Tpex(d21)、Cl-Ttran、Cl-Tex、Cl-Teff

### 3. motif 富集（03_motif_enrichment）

- **cisTarget**：服务器
  `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/cistarget/t_nimfa_cluster_cistarget.ipynb`
  （pycistarget + 预计算数据库
  `/dshare/home/jinhx/database/pycistarget/mm10/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather`）
- **HOMER 结果整理**：服务器
  `/dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis/homer/Untitled.ipynb`
  （读取各簇 `knownResults.txt`，按 q-value < 0.05 过滤）
- **MEME-ChIP / HOMER 结果汇总** `summary_motif.ipynb`：该 notebook 只在本地存档中，
  服务器上未找到，见本仓库 `03_motif_enrichment/summary_motif.ipynb`

### 4. GO 富集与注释（04_go_and_annotation）

- peak 注释：ChIPseeker / HOMER annotatePeaks / GREAT
- `agg_into_xlsx.ipynb`（每簇 GO 结果合并为 Excel）、`unique_gene.ipynb`
  （每簇注释基因 Venn / UpSet、簇特异基因提取）：这两个 notebook 同样只在本地存档中，
  服务器上未找到，见本仓库 `04_go_and_annotation/`

## 环境依赖

**Python 包**：nimfa、fastcluster、pandas、numpy、scipy、scikit-learn、matplotlib、
seaborn、openpyxl、upsetplot、matplotlib-venn、pyranges、pycistarget、jupyter

**外部工具**：deepTools、HOMER、MEME suite（MEME-ChIP）、SLURM（集群提交）

**参考基因组**：mm10

## 使用方式

集群上（SLURM）：

```bash
cd /dshare/home/jinhx/Tcell_project/cut_tag/nimfa/analysis
bash sbatch_t_consensus.sh                          # 提交 NMF 主分析（10 核 / 400G / 7 天）
bash sbatch_t_deeptools.sh <mark> <cluster>         # 如: bash sbatch_t_deeptools.sh H3K27ac 2
```

本地/notebook 步骤按 `01 → 03 → 04` 顺序逐个运行 notebook 即可。

## 注意

- 脚本中的 bigWig 路径、feather 数据库路径为集群本地路径，使用前请改为自己的路径
- 大型中间结果（`*.pkl`、`*.rds`、`*.bw`、矩阵文件）未纳入版本管理，见 `.gitignore`
