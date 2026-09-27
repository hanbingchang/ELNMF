setwd("D:/Papercode/my/GR")

# ---------------- 加载包 ----------------
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(Matrix)  
library(mclust)
library(aricode)
library(dplyr)

# ---------------- 多随机种子 ----------------
seeds <- 1:20
k_cutoff_values <- c(27)

# ---------------- 结果存储 ----------------
all_results <- data.frame(
  seed = integer(),
  k_cutoff = integer(),
  ARI = numeric(),
  NMI = numeric(),
  n_kept = integer(),
  stringsAsFactors = FALSE
)

# ================= ⭐ 新增：标签存储 =================
pred_list <- list()
run_info <- data.frame(
  run_id = integer(),
  seed = integer(),
  k_cutoff = integer(),
  ARI = numeric(),
  NMI = numeric(),
  stringsAsFactors = FALSE
)

run_id <- 1

# ================= 数据加载 =================
cat("========== 数据加载 ==========\n")
data_path <- "D:/Papercode/my/GR/data/STAR/filtered_feature_bc_matrix_fixed_csc.h5"
st_data <- Read10X_h5(data_path)

seurat_obj <- CreateSeuratObject(
  counts = st_data,
  project = "MVC",
  min.cells = 0,
  min.features = 0
)

# ---------------- 基因过滤 ----------------
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
zero_ratio <- Matrix::rowSums(counts_mat == 0) / ncol(counts_mat)
genes_keep <- names(zero_ratio)[zero_ratio < 0.86]
seurat_obj <- subset(seurat_obj, features = genes_keep)

# ---------------- 归一化 + log ----------------
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
spot_totals <- Matrix::colSums(counts_mat)
spot_totals[spot_totals == 0] <- NA

norm_mat <- t(t(counts_mat) / spot_totals) * 1e4
log_norm_mat <- log1p(norm_mat)

seurat_obj <- SetAssayData(
  object = seurat_obj,
  assay = "RNA",
  layer = "data",
  new.data = log_norm_mat
)

# ---------------- HVG ----------------
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
X_full <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[VariableFeatures(seurat_obj), ]

# ---------------- Spot过滤 ----------------
coord_df <- read.csv("D:/Papercode/my/GR/data/STAR/tissue_positions_list.csv", header = FALSE)
colnames(coord_df) <- c("barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol")

tissue_spots <- coord_df$barcode[coord_df$in_tissue == 1]

metadata <- read.csv("D:/Papercode/my/GR/data/STAR/metadata.csv")
label_map <- setNames(metadata$label, metadata$ID)
label_map[label_map %in% c("", "NA", "na")] <- NA

valid_labels <- label_map[!is.na(label_map)]
labeled_spots <- names(valid_labels)

keep_spots <- colnames(X_full) %in% tissue_spots & colnames(X_full) %in% labeled_spots

X_full <- X_full[, keep_spots]
true_labels <- factor(valid_labels[colnames(X_full)])

# ================= 提前加载函数（避免循环重复source） =================
source("D:/Papercode/my/GR/混合图.R") 
source("D:/Papercode/my/GR/随机NMF.R")
source("D:/Papercode/my/GR/集成聚类.R")
source("D:/Papercode/my/GR/综合指标.R")
source("D:/Papercode/my/GR/ensemble_clustering_filtered.R")

# ================= 外层随机种子循环 =================
for (seed in seeds) {
  
  set.seed(seed)
  cat("\n=============================\n")
  cat(sprintf(">>> Seed = %d\n", seed))
  cat("=============================\n")
  
  X <- X_full
  
  for (k_cutoff in k_cutoff_values) {
    
    cat(sprintf("\n>>> k_cutoff = %d\n", k_cutoff))
    
    # ---------------- 混合图 ----------------
    hybrid_df <- hybrid_graph(
      X = X,
      coord_file = "D:/Papercode/my/GR/data/STAR/tissue_positions_list.csv",
      model = "KNN",
      k_cutoff = k_cutoff,
      alpha = 0.8,
      sigma = 1
    )
    
    # ---------------- 构建W ----------------
    current_barcodes <- colnames(X)
    idx_1 <- match(hybrid_df$Cell1_barcode, current_barcodes)
    idx_2 <- match(hybrid_df$Cell2_barcode, current_barcodes)
    
    W <- sparseMatrix(
      i = idx_1,
      j = idx_2,
      x = hybrid_df$Similarity,
      dims = c(ncol(X), ncol(X))
    )
    W <- (W + t(W)) / 2
    
    # ---------------- 集成聚类 ----------------
    n_clusters <- length(unique(true_labels))
    
    ensemble_res <- ensemble_clustering_filtered(
      X = X,
      W = W,
      k_values = seq(10, 28, 2),
      lambda_values = seq(180, 200, 5),
      n_clusters = n_clusters,
      q = 0.3
    )
    
    kept_S <- ensemble_res$all_S_matrices
    S_consensus <- Reduce("+", kept_S) / length(kept_S)
    
    # ---------------- 层次聚类 ----------------
    hc <- hclust(as.dist(1 - S_consensus), method = "average")
    pred_labels <- cutree(hc, k = n_clusters)
    
    # ⭐ 强制加名字（非常关键）
    names(pred_labels) <- colnames(X)
    
    # ⭐ 严格对齐
    pred_labels <- pred_labels[colnames(X)]
    if (any(is.na(pred_labels))) {
      warning("出现NA标签，本次run跳过")
      next
    }
    # ---------------- 评估 ----------------
    ari <- ARI(true_labels, pred_labels)
    nmi <- NMI(true_labels, pred_labels)
    
    cat(sprintf("ARI=%.4f | NMI=%.4f\n", ari, nmi))
    
    # ---------------- 保存指标 ----------------
    all_results <- rbind(
      all_results,
      data.frame(
        seed = seed,
        k_cutoff = k_cutoff,
        ARI = ari,
        NMI = nmi,
        n_kept = ensemble_res$n_kept
      )
    )
    
    # ================= ⭐ 保存标签 =================
    pred_list[[run_id]] <- pred_labels
    
    run_info <- rbind(
      run_info,
      data.frame(
        run_id = run_id,
        seed = seed,
        k_cutoff = k_cutoff,
        ARI = ari,
        NMI = nmi
      )
    )
    
    run_id <- run_id + 1
  }
}

# ================= 统计分析 =================
summary_stats <- all_results %>%
  group_by(k_cutoff) %>%
  summarise(
    ARI_mean = mean(ARI),
    ARI_sd   = sd(ARI),
    ARI_var  = var(ARI),
    NMI_mean = mean(NMI),
    NMI_sd   = sd(NMI),
    NMI_var  = var(NMI),
    n = n()
  ) %>%
  mutate(
    t_value = qt(0.975, df = n - 1),
    ARI_ci_lower = ARI_mean - t_value * ARI_sd / sqrt(n),
    ARI_ci_upper = ARI_mean + t_value * ARI_sd / sqrt(n),
    NMI_ci_lower = NMI_mean - t_value * NMI_sd / sqrt(n),
    NMI_ci_upper = NMI_mean + t_value * NMI_sd / sqrt(n)
  )

print(summary_stats)

# ================= 保存结果 =================
write.csv(summary_stats, "D:/Papercode/my/GR/results/MVC/robustness_results.csv", row.names = FALSE)

write.csv(all_results, "D:/Papercode/my/GR/results/MVC/all_seed_results.csv", row.names = FALSE)

# ⭐ 保存标签（最重要）
saveRDS(pred_list, "D:/Papercode/my/GR/results/MVC/pred_labels_list.rds")

write.csv(run_info, "D:/Papercode/my/GR/results/MVC/run_info_with_labels.csv", row.names = FALSE)