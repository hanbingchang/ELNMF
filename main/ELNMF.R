setwd("D:/Papercode/my/GR")
getwd()

# Human_Breast_Cancer data analysis

# Load required packages --------------------------------------------------
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(Matrix)
library(mclust)
library(aricode)

set.seed(1234)

# ------------------------------ Load data ---------------------------------
cat("========== Data loading ==========\n")
data_path <- "D:/Papercode/my/GR/data/Human_Breast_Cancer/filtered_feature_bc_matrix.h5"
st_data <- Read10X_h5(data_path)
seurat_obj <- CreateSeuratObject(
  counts = st_data,
  project = "HBC",
  min.cells = 0,
  min.features = 0
)
cat("Original data dimensions (genes x spots):", dim(seurat_obj), "\n")

# ------------------------------ Low-quality filtering --------------------
cat("========== Low-quality filtering ==========\n")
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
n_spots <- ncol(counts_mat)
zero_ratio <- Matrix::rowSums(counts_mat == 0) / n_spots
tau <- 0.91
genes_keep <- names(zero_ratio)[zero_ratio < tau]
seurat_obj <- subset(seurat_obj, features = genes_keep)
cat("Genes after zero-ratio filtering:", nrow(seurat_obj), "\n")

# ------------------------------ Normalization ----------------------------
cat("========== Normalization ==========\n")
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
spot_totals <- Matrix::colSums(counts_mat)
spot_totals[spot_totals == 0] <- NA
scale_factor <- 1e4
norm_mat <- t(t(counts_mat) / spot_totals) * scale_factor

# ------------------------------ Log transformation -----------------------
log_norm_mat <- log1p(norm_mat)
seurat_obj <- SetAssayData(
  object = seurat_obj,
  assay = "RNA",
  layer = "data",
  new.data = log_norm_mat
)

# ------------------------------ Select highly variable genes -------------
cat("========== Selecting highly variable genes ==========\n")
seurat_obj <- FindVariableFeatures(seurat_obj,
                                   selection.method = "vst",
                                   nfeatures = 2000)
X <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[VariableFeatures(seurat_obj), ]
cat("HVG expression matrix dimensions:", dim(X), "\n")

# ------------------------------ Align spots ------------------------------
# Keep only spots that are in tissue and have a ground-truth label
cat("========== Spot alignment ==========\n")

coord_file <- "D:/Papercode/my/GR/data/Human_Breast_Cancer/spatial/tissue_positions_list.csv"
coord_df <- read.csv(coord_file, header = FALSE, stringsAsFactors = FALSE)
colnames(coord_df) <- c("barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol")
tissue_spots <- coord_df$barcode[coord_df$in_tissue == 1]
cat("Spots in tissue:", length(tissue_spots), "\n")

meta_path <- "D:/Papercode/my/GR/data/Human_Breast_Cancer/metadata.tsv"
metadata <- read.delim(meta_path, stringsAsFactors = FALSE)
label_map <- setNames(metadata$ground_truth, metadata$ID)
na_patterns <- c("", "NA", "na", "Na", "nA")
for (pattern in na_patterns) {
  label_map[label_map == pattern] <- NA
}
valid_labels <- label_map[!is.na(label_map)]
labeled_spots <- names(valid_labels)
cat("Spots with ground-truth labels:", length(labeled_spots), "\n")

all_spots <- colnames(X)
final_spots <- intersect(intersect(all_spots, tissue_spots), labeled_spots)
idx_keep <- match(final_spots, all_spots)

cat("\n=== Spot filtering summary ===\n")
cat("Spots in expression matrix:", length(all_spots), "\n")
cat("Spots in tissue:", length(tissue_spots), "\n")
cat("Spots with labels:", length(labeled_spots), "\n")
cat("Final spots kept (tissue + labeled):", length(final_spots), "\n")
cat("Filtered out:", length(all_spots) - length(final_spots), 
    sprintf("(%.1f%%)", (length(all_spots) - length(final_spots)) / length(all_spots) * 100), "\n")

X <- X[, idx_keep, drop = FALSE]
cat("Filtered X dimensions:", dim(X), "\n")

true_labels <- valid_labels[final_spots]
true_labels <- factor(true_labels)
cat("True classes:", paste(levels(true_labels), collapse=", "), "\n")
cat("Number of true clusters:", nlevels(true_labels), "\n")

# ------------------------------ Build hybrid graph -----------------------
cat("========== Building hybrid graph ==========\n")
source("D:/Papercode/my/CL-NMF/hybrid graph.R")   # function name: hybrid_graph

hybrid_df <- hybrid_graph(
  X = X,
  coord_file = "D:/Papercode/my/GR/data/Human_Breast_Cancer/spatial/tissue_positions_list.csv",
  model = "KNN",
  k_cutoff = 12,
  alpha = 0.8,
  sigma = 1,
  verbose = TRUE
)

n_spots_total <- max(hybrid_df$Cell1_index, hybrid_df$Cell2_index)
W <- sparseMatrix(
  i = hybrid_df$Cell1_index,
  j = hybrid_df$Cell2_index,
  x = hybrid_df$Similarity,
  dims = c(n_spots_total, n_spots_total)
)
W <- (W + t(W)) / 2
stopifnot(nrow(W) == ncol(X), ncol(W) == ncol(X))
cat("Hybrid adjacency matrix built, dimensions:", dim(W), "\n")

# ------------------------------ Set parameters ---------------------------
cat("========== Parameter setting ==========\n")

k_values <- seq(10, 28, by = 2)
lambda_values <- c(40, 45, 50, 55, 60)

max_iter <- 500
tol <- 1e-4
n_clusters <- length(unique(na.omit(true_labels)))

cat("k_values:", k_values, "\n")
cat("lambda_values:", lambda_values, "\n")

source("D:/Papercode/my/GR/NMF.R")   

# ------------------------------ Generate and save V matrices -------------
cat("\n========== Generating and saving all V matrices ==========\n")

save_dir <- "D:/Papercode/my/GR/results/HBC"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Created directory:", save_dir, "\n")
}

V_list <- list()
V_metadata <- data.frame()
counter <- 1
start_time <- Sys.time()

for (lambda1 in lambda_values) {
  for (k in k_values) {
    cat(sprintf("Generating V [%d/%d]: lambda=%d, k=%d\n", 
                counter, length(k_values) * length(lambda_values),
                lambda1, k))
    
    nmf_res <- regularized_nmf(
      X = X,
      W = W,
      k = k,
      lambda1 = lambda1,
      max_iter = max_iter,
      tol = tol,
      verbose = FALSE
    )
    
    V <- nmf_res$V
    v_name <- sprintf("V_lambda%d_k%d", lambda1, k)
    V_list[[v_name]] <- V
    
    V_metadata <- rbind(V_metadata, data.frame(
      id = v_name,
      lambda = lambda1,
      k = k,
      n_rows = nrow(V),
      n_cols = ncol(V),
      file_name = paste0(v_name, ".rds"),
      stringsAsFactors = FALSE
    ))
    
    counter <- counter + 1
  }
}

saveRDS(V_list, file.path(save_dir, "all_V_matrices.rds"))
write.csv(V_metadata, 
          file.path(save_dir, "V_metadata.csv"),
          row.names = FALSE)

end_time <- Sys.time()
time_diff <- difftime(end_time, start_time, units = "mins")

cat(sprintf("\n========== V matrices saved ==========\n"))
cat(sprintf("Total V matrices generated: %d\n", length(V_list)))
cat(sprintf("Save directory: %s\n", save_dir))
cat(sprintf("Main files:\n"))
cat(sprintf("  - all_V_matrices.rds\n"))
cat(sprintf("  - V_metadata.csv\n"))
cat(sprintf("Time taken: %.2f minutes\n", time_diff))

# ------------------------------ Load external functions ------------------
source("D:/Papercode/my/GR/NMF.R")
source("D:/Papercode/my/GR/ensemble_clustering_filtered.R")  
source("D:/Papercode/my/GR/aggregative_indicator.R")  
source("D:/Papercode/my/GR/filter_top.R")    # may be used inside ensemble_clustering_filtered

# ------------------------------ Ensemble clustering with filtering --------
cat("========== Ensemble clustering with quality filtering ==========\n")

ensemble_res <- ensemble_clustering_filtered(
  X = X,
  V_dir = save_dir,
  k_values = k_values,
  lambda_values = lambda_values,
  n_clusters = n_clusters,
  clustering_methods = c("kmeans", "pam"),
  q = 0.45,          # quantile threshold for keeping
  verbose = TRUE
)

kept_S_matrices <- ensemble_res$all_S_matrices

if (length(kept_S_matrices) > 0) {
  S_consensus <- Reduce("+", kept_S_matrices) / length(kept_S_matrices)
  S_consensus <- as.matrix(S_consensus)
  cat("Number of kept high-quality clusterings:", ensemble_res$n_kept, "\n")
} else {
  stop("No clustering passed quality filtering! Consider lowering q.")
}

# ------------------------------ Hierarchical clustering on consensus -----
cat("\n========== Hierarchical clustering ==========\n")

hierarchical_clustering <- function(S, n_clusters, method = "average") {
  dist_matrix <- as.dist(1 - S)
  hc <- hclust(dist_matrix, method = method)
  pred_labels <- cutree(hc, k = n_clusters)
  return(list(
    labels = pred_labels,
    hclust_obj = hc,
    method = method
  ))
}

hc_result <- hierarchical_clustering(
  S = S_consensus,
  n_clusters = n_clusters,
  method = "average"
)

pred_labels <- hc_result$labels
names(pred_labels) <- colnames(X)

# ------------------------------ Evaluation --------------------------------
cat("\n========== Clustering evaluation ==========\n")

stopifnot(length(true_labels) == length(pred_labels))
cat(sprintf("Number of evaluated spots: %d\n", length(true_labels)))

ari <- ARI(true_labels, pred_labels)
nmi <- NMI(true_labels, pred_labels)

cat(sprintf("Clustering performance: ARI = %.4f, NMI = %.4f\n", ari, nmi))