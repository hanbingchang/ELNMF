library(FNN)
library(dbscan)
library(Matrix)
library(igraph)
library(cluster)

gaussian_kernel <- function(distances, sigma = 1.0) {
  exp(-(distances^2) / (2 * sigma^2))
}

cosine_similarity_01 <- function(vec1, vec2) {
  if (all(vec1 == 0) || all(vec2 == 0)) return(0)
  cos_sim <- sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
  (cos_sim + 1) / 2
}

hybrid_graph <- function(
    X, 
    coord_file,
    model = "Radius", 
    rad_cutoff = NULL, 
    k_cutoff = NULL, 
    alpha = 0.8, 
    sigma = 1.0, 
    verbose = TRUE,
    return_intermediate = FALSE
) {
  stopifnot(model %in% c("Radius", "KNN"))
  stopifnot(alpha >= 0 && alpha <= 1)
  beta <- 1 - alpha
  if (verbose) cat("------ Calculating hybrid graph...\n")
  
  # Load spatial coordinates
  coord_df <- read.csv(coord_file, header = FALSE, stringsAsFactors = FALSE)
  if (!grepl("^[0-9.+-]", coord_df[1, 1])) {
    coord_df <- coord_df[-1, ]
  }
  
  if (ncol(coord_df) == 6) {
    colnames(coord_df) <- c("barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol")
  } else if (ncol(coord_df) == 4) {
    colnames(coord_df) <- c("barcode", "imagerow", "imagecol", "label")
    coord_df$in_tissue <- 1
    coord_df$array_row <- 0
    coord_df$array_col <- 0
    coord_df <- coord_df[, c("barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol")]
  } else {
    stop("Unexpected number of columns in coordinate file. Expected 4 or 6 columns.")
  }
  
  coord_df$imagerow <- as.numeric(coord_df$imagerow)
  coord_df$imagecol <- as.numeric(coord_df$imagecol)
  coor <- coord_df[, c("imagerow", "imagecol")]
  rownames(coor) <- coord_df$barcode
  
  # Align expression matrix with coordinates
  spot_names <- colnames(X)
  if (is.null(spot_names)) stop("Expression matrix X must have column names (spot barcodes).")
  common_spots <- intersect(spot_names, rownames(coor))
  if (length(common_spots) == 0) stop("No common spots between X and coordinate file!")
  
  X <- X[, common_spots, drop = FALSE]
  coor_matrix <- coor[colnames(X), , drop = FALSE]
  
  # Build neighbor relationships
  if (model == "Radius") {
    if (is.null(rad_cutoff)) stop("rad_cutoff must be provided when model='Radius'")
    if (verbose) cat("Calculating distance matrix...\n")
    dist_mat <- as.matrix(dist(coor_matrix))
    indices_euclidean <- vector("list", nrow(coor_matrix))
    distances_euclidean <- vector("list", nrow(coor_matrix))
    for (i in seq_len(nrow(coor_matrix))) {
      neighbors <- which(dist_mat[i, ] <= rad_cutoff & seq_len(nrow(dist_mat)) != i)
      indices_euclidean[[i]] <- neighbors
      distances_euclidean[[i]] <- dist_mat[i, neighbors]
    }
  } else if (model == "KNN") {
    if (is.null(k_cutoff)) stop("k_cutoff must be provided when model='KNN'")
    if (verbose) cat("Finding KNN neighbors...\n")
    knn_result <- FNN::get.knnx(coor_matrix, coor_matrix, k = k_cutoff + 1, algorithm = "kd_tree")
    indices_euclidean <- lapply(seq_len(nrow(coor_matrix)), function(i) knn_result$nn.index[i, -1])
    distances_euclidean <- lapply(seq_len(nrow(coor_matrix)), function(i) knn_result$nn.dist[i, -1])
  }
  
  # Spatial similarity matrix (EW)
  if (verbose) cat("Building spatial similarity matrix...\n")
  spatial_neighbors_list <- list()
  for (it in seq_along(indices_euclidean)) {
    idxs <- indices_euclidean[[it]]
    dists <- distances_euclidean[[it]]
    if (length(idxs) == 0) next
    sims <- gaussian_kernel(dists, sigma)
    df <- data.frame(Cell1 = it, Cell2 = idxs, Similarity_euclidean = sims, stringsAsFactors = FALSE)
    spatial_neighbors_list[[it]] <- df
  }
  EW_ij <- do.call(rbind, spatial_neighbors_list)
  
  # Expression similarity matrix (CW)
  if (verbose) cat("Building expression similarity matrix...\n")
  expr_mat <- t(X)
  expr_sim_list <- list()
  for (it in seq_along(indices_euclidean)) {
    idxs <- indices_euclidean[[it]]
    if (length(idxs) == 0) next
    expr_sims <- sapply(idxs, function(j) cosine_similarity_01(expr_mat[it, ], expr_mat[j, ]))
    df <- data.frame(Cell1 = it, Cell2 = idxs, Similarity_expr = expr_sims, stringsAsFactors = FALSE)
    expr_sim_list[[it]] <- df
  }
  CW_ij <- do.call(rbind, expr_sim_list)
  
  # Merge and combine similarities
  if (verbose) cat("Merging similarity matrices...\n")
  W_ij <- merge(EW_ij, CW_ij, by = c("Cell1", "Cell2"))
  W_ij$Similarity <- alpha * W_ij$Similarity_euclidean + beta * W_ij$Similarity_expr
  
  # Remove duplicates and average
  if (verbose) cat("Processing edges...\n")
  hybrid_graph_df <- W_ij[W_ij$Similarity > 0, c("Cell1", "Cell2", "Similarity"), drop = FALSE]
  hybrid_graph_df$pair <- apply(hybrid_graph_df[, c("Cell1", "Cell2")], 1, function(x) paste(sort(x), collapse = "_"))
  hybrid_graph_df <- aggregate(Similarity ~ pair, data = hybrid_graph_df, FUN = mean)
  tmp <- do.call(rbind, strsplit(hybrid_graph_df$pair, "_"))
  hybrid_graph_df$Cell1 <- as.integer(tmp[, 1])
  hybrid_graph_df$Cell2 <- as.integer(tmp[, 2])
  hybrid_graph_df$pair <- NULL
  
  # Map back to barcodes
  id_to_spot <- colnames(X)
  hybrid_graph_df$Cell1_barcode <- id_to_spot[hybrid_graph_df$Cell1]
  hybrid_graph_df$Cell2_barcode <- id_to_spot[hybrid_graph_df$Cell2]
  hybrid_graph_df$Cell1_index <- hybrid_graph_df$Cell1
  hybrid_graph_df$Cell2_index <- hybrid_graph_df$Cell2
  hybrid_graph_df$Cell1 <- NULL
  hybrid_graph_df$Cell2 <- NULL
  
  if (verbose) {
    n_edges <- nrow(hybrid_graph_df)
    n_cells <- ncol(X)
    cat(sprintf("The graph contains %d edges and %d cells.\n", n_edges, n_cells))
    cat(sprintf("%.4f neighbors per cell on average.\n", n_edges / n_cells))
  }
  
  if (return_intermediate) {
    return(list(
      hybrid_graph = hybrid_graph_df,
      EW_ij = EW_ij,
      CW_ij = CW_ij,
      W_ij = W_ij,
      indices = indices_euclidean,
      distances = distances_euclidean
    ))
  } else {
    return(hybrid_graph_df)
  }
}