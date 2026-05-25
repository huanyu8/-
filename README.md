# 单细胞转录组数据分析流程复现（GSE267764）

## 项目简介

本项目复现论文：

Kimura K, Subramanian A, Yin Z, Khalinezhad A et al.  
*Immune checkpoint TIM-3 regulates microglia and Alzheimer's disease.*  
Nature, 2025, 641(8063):718–731.

数据来源：  
GEO 数据集 GSE267764

- GEO 数据链接：  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267764


本项目主要复现论文中的单细胞转录组（scRNA-seq）分析流程，包括：

- Seurat 数据预处理
- Harmony 批次矫正
- UMAP 聚类分析
- 细胞类型注释
- DoubletFinder 双细胞检测
- 差异表达分析（DGE）
- 火山图可视化

由于原始数据规模较大，本项目主要针对论文中的核心单细胞分析流程进行复现与验证。

---

# 方法概述

本项目基于 Seurat v5 单细胞分析框架，对 GSE267764 数据集进行了完整的 scRNA-seq 标准分析流程复现。

相较于传统 Seurat v4 教程，本复现特别针对 Seurat v5 的 Layer 机制进行了兼容性修复，包括：

- 使用 `LayerData()` 代替旧版 counts 提取方式
- 使用 `JoinLayers()` 解决差异分析中的多图层问题
- 修复 DoubletFinder 与 Seurat v5 的兼容问题
- Harmony 批次校正整合多样本
- 构建 microglia 亚群注释体系
- 复现论文中的 DEG 火山图

---

# 复现结果

| 模块 | 复现结果 |
|---|---|
| 数据类型 | 单细胞转录组（scRNA-seq） |
| 分析框架 | Seurat v5 |
| 批次矫正 | Harmony |
| 双细胞检测 | DoubletFinder |
| 聚类参数 | npc = 50, resolution = 0.4 |
| 差异分析方法 | Wilcoxon rank-sum test |
| 可视化结果 | UMAP、火山图 |
| 输出文件 | DEG_of_HMG.pdf |

成功复现内容包括：

- UMAP 聚类结果
- 小胶质细胞（Microglia）亚群划分
- DoubletFinder 双细胞识别
- 差异表达基因筛选
- HMG vs CC HMG1 火山图

---

# 小组信息
| 成员 | 学号 | 分工 |github|
|---|---|---|---|
| 邹林 | 2025303110139 | 合适论文寻找、初始框架建立、README 编写 | [@huanyu8](https://github.com/huanyu8) |
| 林观清 | 2025303120114 | 论文图复现 | [@LGQ559492](https://github.com/LGQ559492) |
| 郑永琪 | 2025303120147 | 论文图复现 | [@zcguksjv65](https://github.com/zcguksjv65)  |
| 李迁 | 2025303120110 | 论文图复现 |[li20010229-eng](https://github.com/li20010229-eng)
| 黄子凯 | 2025303110141 | 项目整理、Quarto 文档转换、最终提交 |[@hirasauayui](https://github.com/hirasauayui) |

---

# 目录结构

```text
scRNAseq-GSE267764/
├── README.md
├── requirements.txt
├── environment.yml
├── .gitignore
├── data/
│   └── README.md                  # 数据获取说明
├── notebooks/
│   ├── 01_preprocessing.qmd       # 数据预处理
│   ├── 02_harmony_cluster.qmd     # Harmony 聚类分析
│   ├── 03_cell_annotation.qmd     # 细胞类型注释
│   ├── 04_doubletfinder.qmd       # 双细胞检测
│   ├── 05_DGE_analysis.qmd        # 差异表达分析
│   └── 06_volcano_plot.qmd        # 火山图绘制
├── outputs/
│   ├── figures/
│   └── DEG_of_HMG.pdf
└── scripts/
    └── utils.R
```

---

# 完整复现步骤

## 第一步：克隆仓库

```bash
git clone https://github.com/your-repo/scRNAseq-GSE267764.git
cd scRNAseq-GSE267764
```

---

## 第二步：安装 R 环境

本项目使用：

```text
R version 4.6.0
```

推荐使用：

- RStudio
- Quarto
- Rtools（Windows）

---

## 第三步：安装依赖包

在 R 中运行：

```r
install.packages(c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "patchwork",
  "harmony",
  "DoubletFinder",
  "tibble",
  "ggrepel"
))
```

---

## 第四步：下载 GEO 数据

从 GEO 下载原始数据：

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267764

解压后目录结构如下：

```text
GSE267764_RAW/
├── GSMxxxxxxx_barcodes.tsv.gz
├── GSMxxxxxxx_features.tsv.gz
├── GSMxxxxxxx_matrix.mtx.gz
└── ...
```

---

## 第五步：整理 10X 数据目录

运行如下代码：

```r
data_dir <- "你的数据路径/GSE267764_RAW"

fs <- list.files(data_dir, pattern = "GSM", full.names = FALSE)

samples <- substr(fs, 1, 10)

lapply(unique(samples), function(x) {

  y <- fs[grepl(x, fs)]

  folder <- file.path(
    data_dir,
    paste(strsplit(y[1], split = "_")[[1]][1:2], collapse = "_")
  )

  dir.create(folder, recursive = TRUE, showWarnings = FALSE)

  file.rename(file.path(data_dir, y[1]), file.path(folder, "barcodes.tsv.gz"))
  file.rename(file.path(data_dir, y[2]), file.path(folder, "features.tsv.gz"))
  file.rename(file.path(data_dir, y[3]), file.path(folder, "matrix.mtx.gz"))
})
```

⚠️ 注意：

该步骤第二次运行时会报错。  
原因是原始文件已经被移动。

解决办法：

- 重新解压 GEO 原始文件
- 或备份原始数据目录

---

## 第六步：运行 Seurat 预处理

主要步骤包括：

- 创建 Seurat 对象
- 合并样本
- 线粒体比例计算
- 质量控制过滤

过滤条件：

```r
nCount_RNA >= 1000
nFeature_RNA >= 200
percent.mt < 80
```

---

## 第七步：Harmony 批次矫正与聚类

核心参数：

```r
npc = 50
resolution = 0.4
```

主要流程：

```r
NormalizeData()
FindVariableFeatures()
ScaleData()
RunPCA()
RunHarmony()
RunUMAP()
FindNeighbors()
FindClusters()
```

随后去除低质量 Cluster6：

```r
scRNAseq_cortex <- subset(
  scRNAseq_cortex,
  seurat_clusters %in% c(0:5, 7:10)
)
```

---

## 第八步：细胞类型注释

根据原论文定义 microglia 亚群：

| Cluster | Cell Type |
|---|---|
| 0 | HMG |
| 1 | HMG |
| 2 | MGnD 0 |
| 3 | MGnD1 |
| 4 | CC HMG1 |
| 5 | CC HMG2 |
| 7 | IFN HMG |
| 8 | IFN HMGD |
| 9 | CC MGnD |
| 10 | MGnD 2 |

---

## 第九步：DoubletFinder 双细胞检测

本项目特别修复了 Seurat v5 与 DoubletFinder 的兼容问题。

关键修改：

### 修复 1：使用 LayerData 提取 counts

```r
counts_sub <- LayerData(
  sc_sub_cortex,
  assay = "RNA",
  layer = "counts"
)
```

### 修复 2：重新创建单层 Seurat 对象

```r
sc_sub_simple <- CreateSeuratObject(
  counts = counts_sub,
  meta.data = metadata_sub
)
```

### 修复 3：运行 DoubletFinder

```r
sc_sub_simple <- doubletFinder(
  sc_sub_simple,
  PCs = 1:50,
  pK = pK,
  nExp = nExp_poi.adj
)
```

---

## 第十步：差异表达分析（DGE）

使用：

```r
FindMarkers()
```

统计方法：

```text
Wilcoxon rank-sum test
```

关键修复：

```r
subset_obj <- JoinLayers(subset_obj)
```

用于解决 Seurat v5 多 Layer 差异分析问题。

显著性筛选条件：

```r
p_val_adj < 0.05
abs(avg_log2FC) > 0.25
pct.1 > 0.3
```

---

## 第十一步：绘制火山图

复现内容：

- HMG vs CC HMG1
- 上调基因
- 下调基因
- 基因标签标注

输出文件：

```text
DEG_of_HMG.pdf
```

---

# 常见问题

## Q1：Read10X 无法读取数据？

A：请确认目录结构为：

```text
sample/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

---

## Q2：Harmony 报错？

A：请确认已安装：

```r
install.packages("harmony")
```

---

## Q3：DoubletFinder 与 Seurat v5 不兼容？

A：必须使用：

```r
LayerData()
JoinLayers()
```

避免 multi-layer 报错。

---

## Q4：第二次运行文件整理代码失败？

A：因为原始文件已被移动。  
请重新解压 GEO 原始数据。

---

## Q5：FindMarkers 报 Layer 错误？

A：需要先执行：

```r
subset_obj <- JoinLayers(subset_obj)
```

---

# 参考文献

```bibtex
@article{kimura2025tim3,
  title={Immune checkpoint TIM-3 regulates microglia and Alzheimer's disease},
  author={Kimura, K and Subramanian, A and Yin, Z and Khalinezhad, A},
  journal={Nature},
  volume={641},
  number={8063},
  pages={718--731},
  year={2025}
}
```
