# Run with Rscript or source(); locate the repository independently of the working directory.
repo_root <- local({
  candidates <- c(Sys.getenv("ELNMF_ROOT", unset = ""), getwd())
  for (frame in sys.frames()) {
    if (!is.null(frame$ofile)) candidates <- c(dirname(frame$ofile), candidates)
  }
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg)) candidates <- c(dirname(sub("^--file=", "", script_arg[[1]])), candidates)
  found <- NULL
  for (candidate in candidates[nzchar(candidates)]) {
    current <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(current, "main", "shared")) && dir.exists(file.path(current, "data"))) {
        found <- current; break
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    if (!is.null(found)) break
  }
  if (is.null(found)) stop("Cannot locate ELNMF repository. Run from its root or set ELNMF_ROOT.")
  found
})
setwd(repo_root)
# Outputs are created only when this script is run, outside main/.


library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(Matrix)  
library(mclust)
library(aricode)
library(dplyr)


seeds <- 1:20
output_dir <- "results/MVC/stability"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
k_cutoff_values <- c(27)


all_results <- data.frame(
  seed = integer(),
  k_cutoff = integer(),
  ARI = numeric(),
  NMI = numeric(),
  n_kept = integer(),
  stringsAsFactors = FALSE
)


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

data_path <- "data/MVC/filtered_feature_bc_matrix.h5"
st_data <- Read10X_h5(data_path)

seurat_obj <- CreateSeuratObject(
  counts = st_data,
  project = "MVC",
  min.cells = 0,
  min.features = 0
)


counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
zero_ratio <- Matrix::rowSums(counts_mat == 0) / ncol(counts_mat)
genes_keep <- names(zero_ratio)[zero_ratio < 0.86]
seurat_obj <- subset(seurat_obj, features = genes_keep)


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

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
X_full <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[VariableFeatures(seurat_obj), ]


coord_df <- read.csv("data/MVC/tissue_positions_list.csv", header = FALSE)
colnames(coord_df) <- c("barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol")

tissue_spots <- coord_df$barcode[coord_df$in_tissue == 1]

metadata <- read.csv("data/MVC/metadata.csv")
label_map <- setNames(metadata$label, metadata$ID)
label_map[label_map %in% c("", "NA", "na")] <- NA

valid_labels <- label_map[!is.na(label_map)]
labeled_spots <- names(valid_labels)

keep_spots <- colnames(X_full) %in% tissue_spots & colnames(X_full) %in% labeled_spots

X_full <- X_full[, keep_spots]
true_labels <- factor(valid_labels[colnames(X_full)])


source("main/shared/hybrid graph.R") 
source("main/shared/NMF.R")
source("main/shared/ensemble_clustering_filtered.R")
source("main/shared/aggregative_indicator.R")
source("main/shared/filter_top.R")


for (seed in seeds) {
  
  set.seed(seed)
  cat("\n=============================\n")
  cat(sprintf(">>> Seed = %d\n", seed))
  cat("=============================\n")
  
  X <- X_full
  
  for (k_cutoff in k_cutoff_values) {
    
    cat(sprintf("\n>>> k_cutoff = %d\n", k_cutoff))
    

    hybrid_df <- hybrid_graph(
  coord_header = "legacy", # Preserve the reader used for this historical experiment.
      X = X,
      coord_file = "data/MVC/tissue_positions_list.csv",
      model = "KNN",
      k_cutoff = k_cutoff,
      alpha = 0.8,
      sigma = 1
    )
    

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
    

    n_clusters <- length(unique(true_labels))
    
    ensemble_res <- ensemble_clustering_online(
      X = X,
      W = W,
      k_values = seq(10, 28, 2),
      lambda_values = seq(180, 200, 5),
      n_clusters = n_clusters,
      seed = seed,
      q = 0.3
    )
    
    kept_S <- ensemble_res$all_S_matrices
    S_consensus <- Reduce("+", kept_S) / length(kept_S)
    

    hc <- hclust(as.dist(1 - S_consensus), method = "average")
    pred_labels <- cutree(hc, k = n_clusters)
    

    names(pred_labels) <- colnames(X)
    

    pred_labels <- pred_labels[colnames(X)]
    if (any(is.na(pred_labels))) {
      warning("出现NA标签，本次run跳过")
      next
    }

    ari <- ARI(true_labels, pred_labels)
    nmi <- NMI(true_labels, pred_labels)
    
    cat(sprintf("ARI=%.4f | NMI=%.4f\n", ari, nmi))
    

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
write.csv(summary_stats, "results/MVC/stability/robustness_results.csv", row.names = FALSE)

write.csv(all_results, "results/MVC/stability/all_seed_results.csv", row.names = FALSE)

# ⭐ 保存标签（最重要）
saveRDS(pred_list, "results/MVC/stability/pred_labels_list.rds")

write.csv(run_info, "results/MVC/stability/run_info_with_labels.csv", row.names = FALSE)