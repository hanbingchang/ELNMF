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

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript HBC_ELNMF_scaling_benchmark.R <target_n> <output_root>")
}

target_n <- as.integer(args[[1]])
output_root <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
if (is.na(target_n) || target_n < 20) stop("target_n must be an integer >= 20")

suppressPackageStartupMessages({
  library(Seurat)
  library(hdf5r)
  library(Matrix)
  library(cluster)
})

# Optional third argument "consensus" reuses the 50 saved V matrices.
if (length(args) >= 3 && identical(args[[3]], "consensus")) {
  run_saved_consensus <- function() {
set.seed(1234)
helper_root <- "main/shared"
run_dir <- file.path(output_root, "work", sprintf("N%05d", target_n))
v_dir <- file.path(run_dir, "V_matrices")
coord_path <- file.path(run_dir, "benchmark_positions.csv")
start_marker <- file.path(run_dir, "CONSENSUS_RESUME_STARTED.marker")
end_marker <- file.path(run_dir, "CONSENSUS_RESUME_COMPLETED.marker")
runtime_path <- file.path(run_dir, "consensus_resume_runtime_seconds.txt")
error_path <- file.path(run_dir, "CONSENSUS_RESUME_ERROR.txt")
unlink(c(start_marker, end_marker, runtime_path, error_path), force = TRUE)

tryCatch({
  v_files <- list.files(v_dir, pattern = "^V_lambda[0-9]+_k[0-9]+\\.rds$", full.names = TRUE)
  if (length(v_files) != 50L) stop(sprintf("Expected 50 V matrices, found %d", length(v_files)))
  coords <- read.csv(coord_path, header = FALSE, stringsAsFactors = FALSE)
  if (nrow(coords) != target_n) stop("Coordinate row count does not match target_n")
  spot_ids <- coords[[1]]
  X_stub <- matrix(0, nrow = 1L, ncol = target_n, dimnames = list("stub", spot_ids))

  writeLines(sprintf("N=%d;started=%s", target_n, format(Sys.time(), "%Y-%m-%d %H:%M:%S")), start_marker)
  t0 <- proc.time()[["elapsed"]]

  source(file.path(helper_root, "ensemble_clustering_filtered.R"))
  source(file.path(helper_root, "aggregative_indicator.R"))
  source(file.path(helper_root, "filter_top.R"))
  ensemble_res <- ensemble_clustering_filtered(
    X = X_stub,
    V_dir = v_dir,
    k_values = seq(10, 28, by = 2),
    lambda_values = c(40, 45, 50, 55, 60),
    n_clusters = 20L,
    clustering_methods = c("kmeans", "pam"),
    q = 0.45,
    verbose = TRUE
  )
  S_consensus <- Reduce("+", ensemble_res$all_S_matrices) / ensemble_res$n_kept
  hc <- hclust(as.dist(1 - as.matrix(S_consensus)), method = "average")
  pred_labels <- cutree(hc, k = 20L)
  names(pred_labels) <- spot_ids

  elapsed <- proc.time()[["elapsed"]] - t0
  writeLines(sprintf("%.6f", elapsed), runtime_path)
  writeLines(sprintf("N=%d;completed=%s", target_n, format(Sys.time(), "%Y-%m-%d %H:%M:%S")), end_marker)
  write.csv(data.frame(spot_id = spot_ids, cluster = as.integer(pred_labels)),
            file.path(run_dir, "final_cluster_labels.csv"), row.names = FALSE)
  write.csv(data.frame(N_spots = target_n, Base_clusterings = 100L,
                       Retained_clusterings = ensemble_res$n_kept,
                       Consensus_runtime_sec = elapsed),
            file.path(run_dir, "consensus_resume_summary.csv"), row.names = FALSE)
}, error = function(e) {
  writeLines(conditionMessage(e), error_path)
  message(conditionMessage(e))
  stop(e)
})


  }
  run_saved_consensus()
} else {
seed_value <- 1234L
set.seed(seed_value)

data_root <- "data/Human_Breast_Cancer"
helper_root <- "main/shared"
nmf_script <- "main/shared/NMF.R"

run_name <- sprintf("N%05d", target_n)
run_dir <- file.path(output_root, "work", run_name)
v_dir <- file.path(run_dir, "V_matrices")
dir.create(v_dir, recursive = TRUE, showWarnings = FALSE)

start_marker <- file.path(run_dir, "TIMED_REGION_STARTED.marker")
end_marker <- file.path(run_dir, "TIMED_REGION_COMPLETED.marker")
error_marker <- file.path(run_dir, "ERROR.txt")
runtime_path <- file.path(run_dir, "runtime_seconds.txt")
summary_path <- file.path(run_dir, "run_summary.csv")
labels_path <- file.path(run_dir, "final_cluster_labels.csv")
coord_path <- file.path(run_dir, "benchmark_positions.csv")
unlink(c(start_marker, end_marker, error_marker, runtime_path), force = TRUE)

tryCatch({
  # ---------------------------------------------------------------------------
  # Input construction is intentionally outside the timed region.
  # ---------------------------------------------------------------------------
  raw_counts <- Read10X_h5(file.path(data_root, "filtered_feature_bc_matrix.h5"))

  coord_df <- read.csv(
    file.path(data_root, "spatial", "tissue_positions_list.csv"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(coord_df) <- c(
    "barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol"
  )

  metadata <- read.delim(file.path(data_root, "metadata.tsv"), stringsAsFactors = FALSE)
  label_map <- setNames(metadata$ground_truth, metadata$ID)
  label_map[trimws(as.character(label_map)) %in% c("", "NA", "na", "Na", "nA")] <- NA
  valid_labels <- label_map[!is.na(label_map)]

  tissue_ids <- coord_df$barcode[coord_df$in_tissue == 1]
  base_ids <- colnames(raw_counts)[
    colnames(raw_counts) %in% tissue_ids & colnames(raw_counts) %in% names(valid_labels)
  ]
  if (length(base_ids) != 3798L) {
    stop(sprintf("Expected 3798 aligned HBC spots, found %d", length(base_ids)))
  }

  base_labels <- factor(valid_labels[base_ids])
  base_coord <- coord_df[match(base_ids, coord_df$barcode), , drop = FALSE]
  if (anyNA(base_coord$barcode)) stop("Coordinate alignment failed")

  if (target_n < length(base_ids)) {
    class_counts <- table(base_labels)
    exact_alloc <- target_n * as.numeric(class_counts) / sum(class_counts)
    alloc <- floor(exact_alloc)
    remainder <- target_n - sum(alloc)
    if (remainder > 0) {
      add_order <- order(exact_alloc - alloc, decreasing = TRUE)
      alloc[add_order[seq_len(remainder)]] <- alloc[add_order[seq_len(remainder)]] + 1L
    }

    selected_ids <- unlist(lapply(seq_along(levels(base_labels)), function(i) {
      ids_i <- base_ids[base_labels == levels(base_labels)[i]]
      sample(ids_i, size = alloc[i], replace = FALSE)
    }), use.names = FALSE)
    selected_ids <- sample(selected_ids, length(selected_ids), replace = FALSE)
    source_ids <- selected_ids
    output_ids <- selected_ids
    target_labels <- valid_labels[source_ids]
    target_coord <- coord_df[match(source_ids, coord_df$barcode), , drop = FALSE]
  } else if (target_n == length(base_ids)) {
    source_ids <- base_ids
    output_ids <- base_ids
    target_labels <- valid_labels[source_ids]
    target_coord <- base_coord
  } else {
    extra_n <- target_n - length(base_ids)
    # Uniform sampling over original spots preserves the observed domain proportions.
    extra_source <- sample(base_ids, size = extra_n, replace = TRUE)
    extra_ids <- sprintf("%s__rep%05d", extra_source, seq_len(extra_n))

    source_ids <- c(base_ids, extra_source)
    output_ids <- c(base_ids, extra_ids)
    target_labels <- valid_labels[source_ids]

    extra_coord <- coord_df[match(extra_source, coord_df$barcode), , drop = FALSE]
    spatial_span <- max(base_coord$imagerow) - min(base_coord$imagerow)
    array_span <- max(base_coord$array_row) - min(base_coord$array_row)
    extra_coord$imagerow <- extra_coord$imagerow + spatial_span + 100
    extra_coord$array_row <- extra_coord$array_row + array_span + 10
    # Deterministic tiny offsets prevent exact coordinate ties among replicates.
    extra_coord$imagecol <- extra_coord$imagecol + seq_len(extra_n) * 1e-6
    extra_coord$array_col <- extra_coord$array_col + seq_len(extra_n) * 1e-6
    extra_coord$barcode <- extra_ids
    target_coord <- rbind(base_coord, extra_coord)
  }

  counts_target <- raw_counts[, match(source_ids, colnames(raw_counts)), drop = FALSE]
  colnames(counts_target) <- output_ids
  target_coord$barcode <- output_ids
  target_coord$in_tissue <- 1L
  write.table(
    target_coord[, c("barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol")],
    file = coord_path,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  if (ncol(counts_target) != target_n || nrow(target_coord) != target_n) {
    stop("Constructed benchmark input has an unexpected size")
  }

  rm(raw_counts, coord_df, metadata, base_coord)
  invisible(gc())

  # ---------------------------------------------------------------------------
  # Timed method region: preprocessing through final label generation.
  # ---------------------------------------------------------------------------
  writeLines(
    sprintf("N=%d;seed=%d;started=%s", target_n, seed_value, format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    start_marker
  )
  timed_start <- proc.time()[["elapsed"]]

  seurat_obj <- CreateSeuratObject(
    counts = counts_target,
    project = sprintf("HBC_scaling_%d", target_n),
    min.cells = 0,
    min.features = 0
  )

  counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  zero_ratio <- Matrix::rowSums(counts_mat == 0) / ncol(counts_mat)
  genes_keep <- names(zero_ratio)[zero_ratio < 0.91]
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
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  X <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")[VariableFeatures(seurat_obj), ]
  stopifnot(ncol(X) == target_n)

  source(file.path(helper_root, "hybrid graph.R"))
  hybrid_df <- hybrid_graph(
    X = X,
    coord_file = coord_path,
    model = "KNN",
    k_cutoff = 12,
    alpha = 0.8,
    sigma = 1,
    verbose = FALSE
  )
  W <- sparseMatrix(
    i = hybrid_df$Cell1_index,
    j = hybrid_df$Cell2_index,
    x = hybrid_df$Similarity,
    dims = c(target_n, target_n)
  )
  W <- (W + t(W)) / 2
  stopifnot(nrow(W) == ncol(X), ncol(W) == ncol(X))

  k_values <- seq(10, 28, by = 2)
  lambda_values <- c(40, 45, 50, 55, 60)
  n_clusters <- 20L
  max_iter <- 500L
  tol <- 1e-4

  source(nmf_script)
  V_list <- list()
  V_metadata <- data.frame()
  counter <- 1L
  for (lambda1 in lambda_values) {
    for (k in k_values) {
      cat(sprintf("N=%d, V[%d/50], lambda=%d, K=%d\n", target_n, counter, lambda1, k))
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
      saveRDS(V, file.path(v_dir, paste0(v_name, ".rds")))
      V_metadata <- rbind(V_metadata, data.frame(
        id = v_name,
        lambda = lambda1,
        k = k,
        n_rows = nrow(V),
        n_cols = ncol(V),
        file_name = paste0(v_name, ".rds")
      ))
      counter <- counter + 1L
    }
  }

  source(file.path(helper_root, "ensemble_clustering_filtered.R"))
  source(file.path(helper_root, "aggregative_indicator.R"))
  source(file.path(helper_root, "filter_top.R"))
  ensemble_res <- ensemble_clustering_filtered(
    X = X,
    V_dir = v_dir,
    k_values = k_values,
    lambda_values = lambda_values,
    n_clusters = n_clusters,
    clustering_methods = c("kmeans", "pam"),
    q = 0.45,
    verbose = TRUE
  )

  S_consensus <- Reduce("+", ensemble_res$all_S_matrices) / ensemble_res$n_kept
  S_consensus <- as.matrix(S_consensus)
  dist_matrix <- as.dist(1 - S_consensus)
  hc <- hclust(dist_matrix, method = "average")
  pred_labels <- cutree(hc, k = n_clusters)
  names(pred_labels) <- colnames(X)

  runtime_seconds <- proc.time()[["elapsed"]] - timed_start
  writeLines(sprintf("%.6f", runtime_seconds), runtime_path)
  writeLines(
    sprintf("N=%d;completed=%s", target_n, format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    end_marker
  )

  # Output writing and evaluation are outside the timed region.
  write.csv(V_metadata, file.path(run_dir, "V_metadata.csv"), row.names = FALSE)
  write.csv(
    data.frame(spot_id = names(pred_labels), cluster = as.integer(pred_labels)),
    labels_path,
    row.names = FALSE
  )
  run_summary <- data.frame(
    Dataset = "HBC",
    N_spots = target_n,
    Source_spots = length(base_ids),
    Input_mode = if (target_n < length(base_ids)) "stratified_subsample" else if (target_n == length(base_ids)) "original" else "original_plus_resampled_spots",
    Seed = seed_value,
    Clusters = n_clusters,
    K_values = "10,12,14,16,18,20,22,24,26,28",
    Lambda_values = "40,45,50,55,60",
    Base_clusterings = 100L,
    Retained_clusterings = ensemble_res$n_kept,
    Runtime_sec = runtime_seconds,
    stringsAsFactors = FALSE
  )
  write.csv(run_summary, summary_path, row.names = FALSE)
  cat(sprintf("SCALING_COMPLETE|N=%d|runtime_sec=%.6f|retained=%d\n", target_n, runtime_seconds, ensemble_res$n_kept))
}, error = function(e) {
  writeLines(conditionMessage(e), error_marker)
  message("SCALING_ERROR|N=", target_n, "|", conditionMessage(e))
  quit(status = 1L, save = "no")
})

}
