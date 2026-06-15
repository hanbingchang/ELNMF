# ELNMF

ELNMF is an ensemble learning-driven framework for spatial domain identification in spatial transcriptomics data. By jointly leveraging spatial context and transcriptional information, ELNMF identifies biologically meaningful spatial domains with enhanced robustness and reliability.

## Overview

ELNMF employs a graph-regularized non-negative matrix factorization framework to learn informative latent representations from spatial transcriptomics data by jointly considering spatial proximity and gene expression similarity. This integration enables the characterization of both tissue architecture and transcriptional heterogeneity underlying spatial organization.

To further enhance the reliability of spatial domain identification, ELNMF introduces an ensemble learning-driven consensus strategy. Multiple latent representations generated under different parameter settings are first used to derive diverse base clustering solutions. High-quality clustering results are subsequently selected and integrated through consensus learning to aggregate complementary clustering information. By leveraging the diversity of latent representations, ELNMF effectively reduces sensitivity to parameter selection and stochastic variation, thereby improving the stability, robustness, and reliability of the final spatial domain assignments.

The framework was implemented in **R (version 4.5.2)**. 

---
## Framework 
<p align="center">
  <img src="ELNMF_overview.png" width="900">
</p>
 
 **Figure 1.** Schematic overview of the ELNMF framework. 
(A) Hybrid graph construction: spatial coordinates and the preprocessed gene expression matrix are integrated to construct a hybrid graph that characterizes local spatial structure and expression associations. (B) Non-negative matrix factorization: the expression matrix is decomposed to learn multiple low-dimensional latent representations, with the constructed hybrid graph incorporated as a graph regularization term. (C) Consensus clustering: multiple low-dimensional representations are integrated, and a consensus matrix is constructed through consensus clustering to obtain stable clustering results. (D) Downstream applications: the optimized low-dimensional representations can be further used for downstream tasks, including spatial domain identification, functional enrichment analysis, and cell–cell communication analysis.
 
 ---
## Requirements

### Software

- R ≥ 4.5.2

### R packages (direct dependencies)

The following R packages are required to run ELNMF:

```r
required_packages <- c(
  "Seurat",   
  "hdf5r",    
  "ggplot2", 
  "patchwork",
  "Matrix",   
  "mclust",   
  "aricode",  
  "FNN",      
  "dbscan",   
  "igraph",  
  "kernlab",  
  "cluster"   
)

install.packages(
  setdiff(required_packages,
          rownames(installed.packages()))
)
```
---
## Description of Scripts

| Script | Description |
|----------|-------------|
| `ELNMF.R` | Main script for running the complete ELNMF pipeline |
| `hybrid_graph.R` | Construction of spatial–expression similarity structure |
| `NMF.R` | Graph-regularized non-negative matrix factorization implementation |
| `ensemble_clustering_filtered.R` | Generation of clustering results from multiple latent representations |
| `aggregative_indicator.R` | Computation of clustering quality score based on Silhouette Score and CHAOS (Silhouette/CHAOS ratio) |
| `filter_top.R` | Selection of high-quality clustering results based on the composite score and removal of low-quality solutions for consensus construction |

---

## Tutorials

Here is a quick example of running ELNMF on the Human Breast Cancer dataset (10X Visium).

```r
library(Seurat)
library(Matrix)
library(hdf5r)  

# 1. Load example dataset
data_path <- "data/Human_Breast_Cancer/"

expr <- Read10X_h5(paste0(data_path, "filtered_feature_bc_matrix.h5"))
coords <- read.csv(paste0(data_path, "spatial/tissue_positions_list.csv"))

# 2. Run ELNMF
source("main/ELNMF.R")

result <- ELNMF(
  expr_matrix = expr,
  spatial_coord = coords,
  k_range = 10:28,
  lambda_range = 40:60,  
  n_cluster = 20          
)

# 3. Save clustering results to CSV
cluster_df <- data.frame(
  barcode = colnames(expr),        
  cluster = result$cluster         
)
write.csv(cluster_df, file = "ELNMF_clustering_results.csv", row.names = FALSE)

# Optional: print first few rows to check
head(cluster_df)
```
---
## License

This project is released under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).  
See the [LICENSE](LICENSE) file for the full text.

