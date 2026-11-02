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

suppressPackageStartupMessages({
  library(Seurat)
  library(rliger)
  library(RANN)
  library(Matrix)
  library(proxy)
})

data_dir <- "data/Human_Breast_Cancer"
output_dir <- Sys.getenv("HBC_BENCH_OUTPUT", unset = "results/HBC_baseline_computational_benchmark/labels")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Input preparation: read and align exactly the 3,798 HBC spots used by ELNMF.
expr <- as.matrix(Read10X_h5(file.path(data_dir, "filtered_feature_bc_matrix.h5")))
meta <- read.delim(file.path(data_dir, "metadata.tsv"), stringsAsFactors = FALSE)
coord <- read.csv(file.path(data_dir, "spatial/tissue_positions_list.csv"),
                  header = FALSE, stringsAsFactors = FALSE)
colnames(coord) <- c("barcode", "in_tissue", "array_row", "array_col", "x", "y")
valid_meta <- meta[!is.na(meta$ground_truth) & nzchar(meta$ground_truth), , drop = FALSE]
tissue_ids <- coord$barcode[coord$in_tissue == 1]
final_spots <- colnames(expr)[colnames(expr) %in% tissue_ids & colnames(expr) %in% valid_meta$ID]
if (length(final_spots) != 3798L) stop("Expected 3798 HBC spots, found ", length(final_spots))
expr <- expr[, final_spots, drop = FALSE]
coord <- coord[match(final_spots, coord$barcode), , drop = FALSE]

cat("BENCHMARK_START|CGNMF|spots=3798|target_clusters=20\n")
flush.console()
timer_start <- proc.time()[["elapsed"]]

# Existing CGNMF preprocessing and parameters from the supplied script.
coord_mat <- as.matrix(coord[, c("x", "y")])
mode(coord_mat) <- "numeric"
coord_mat[, 1] <- coord_mat[, 1] - min(coord_mat[, 1])
coord_mat[, 2] <- coord_mat[, 2] - min(coord_mat[, 2])
coord_mat <- coord_mat / max(coord_mat)

gene_means <- rowMeans(expr)
expr <- expr[gene_means > quantile(gene_means, 0.69), , drop = FALSE]

liger_obj <- createLiger(list(data = expr))
liger_obj <- rliger::normalize(liger_obj)
liger_obj <- rliger::selectGenes(liger_obj)
liger_obj <- rliger::scaleNotCenter(liger_obj)
Y <- as.matrix(liger_obj@datasets$data@scaleData)

# Parameters currently selected by the supplied entry (sample 151672 row).
K <- 10L
kA <- 33L
k1 <- 30L
k2 <- 30L
lambda1 <- 3
lambda2 <- 6
tau <- 0.3
max_iter <- 200L

createA <- function(coords, neighbors) {
  neighbors <- min(neighbors, nrow(coords) - 1L)
  nn <- nn2(coords, k = neighbors + 1L)$nn.idx[, -1, drop = FALSE]
  n <- nrow(coords)
  A <- Matrix(0, n, n, sparse = TRUE)
  for (i in seq_len(ncol(nn))) {
    A <- A + sparseMatrix(i = seq_len(n), j = nn[, i], x = 1, dims = c(n, n))
  }
  as.matrix(A %*% t(A))
}

A <- createA(coord_mat, kA)
D <- diag(rowSums(A))

set.seed(123)
B <- matrix(runif(nrow(Y) * K), nrow = nrow(Y))
P <- matrix(runif(ncol(Y) * K), nrow = ncol(Y))

compute_QW <- function(P) {
  N <- nrow(P)
  Q <- matrix(0, N, K)
  W <- matrix(0, N, K)
  sim_cos <- 1 - as.matrix(proxy::dist(P, method = "cosine"))
  sim_euc <- 1 / (0.001 + as.matrix(dist(coord_mat)))
  for (i in seq_len(N)) {
    idx_p <- order(sim_cos[i, ], decreasing = TRUE)[-1][seq_len(k1)]
    idx_a <- order(sim_euc[i, ], decreasing = TRUE)[-1][seq_len(k1)]
    pos_idx <- intersect(idx_p, idx_a)
    if (!length(pos_idx)) pos_idx <- idx_a
    pos <- P[pos_idx, , drop = FALSE]
    d <- sqrt(rowSums((pos - matrix(P[i, ], nrow(pos), K, byrow = TRUE))^2))
    w <- exp(1 / (d + 1e-6))
    w <- w / sum(w)
    P_pos <- drop(t(w) %*% pos)
    idx_pn <- order(sim_cos[i, ], decreasing = FALSE)[seq_len(k2)]
    idx_an <- order(sim_euc[i, ], decreasing = FALSE)[seq_len(k2)]
    neg_idx <- intersect(idx_pn, idx_an)
    if (!length(neg_idx)) neg_idx <- idx_pn
    # Preserve the supplied stable-gradient implementation.
    Q[i, ] <- (P_pos - P[i, ]) / (tau * N)
    W[i, ] <- (P_pos - P[i, ]) / (tau * N)
  }
  list(Q = Q, W = W)
}

for (iter in seq_len(max_iter)) {
  B <- (Y %*% P) * B / (B %*% t(P) %*% P)
  B[!is.finite(B) | B < 0] <- 0
  tmp <- compute_QW(P)
  P <- (t(Y) %*% B + lambda1 * (A %*% P) + lambda2 * tmp$W / 2) * P /
    (P %*% t(B) %*% B + lambda1 * (D %*% P) + lambda2 * tmp$Q / 2)
  P[!is.finite(P) | P < 0] <- 0
  if (iter %% 50L == 0L) cat("CGNMF iteration", iter, "of", max_iter, "\n")
}

# CGNMF accepts a fixed cluster count directly.
set.seed(123)
labels <- kmeans(P, centers = 20L, nstart = 50)$cluster
runtime_sec <- proc.time()[["elapsed"]] - timer_start
cat(sprintf("BENCHMARK_END|CGNMF|runtime_sec=%.6f|clusters=%d\n",
            runtime_sec, length(unique(labels))))
flush.console()

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(Barcode = final_spots, Cluster = labels),
          file.path(output_dir, "HBC_CGNMF_labels.csv"), row.names = FALSE)

