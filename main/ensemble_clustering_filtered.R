library(igraph)
library(kernlab)
library(cluster)

# Compute association (connectivity) matrix from labels
compute_association_matrix <- function(labels) {
  labels <- as.numeric(as.factor(labels))
  n <- length(labels)
  mat <- matrix(labels, n, n, byrow = FALSE)
  mat_t <- matrix(labels, n, n, byrow = TRUE)
  S <- (mat == mat_t) + 0
  return(S)
}

# K-means clustering on V matrix
perform_kmeans <- function(V, n_clusters, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  V <- as.matrix(V)
  km_result <- kmeans(V, centers = n_clusters, nstart = 25, iter.max = 100)
  return(km_result$cluster)
}

# PAM clustering (fixed number of clusters)
perform_pam <- function(V, n_clusters, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  V <- as.matrix(V)
  pam_result <- pam(x = V, k = n_clusters, diss = FALSE, metric = "euclidean")
  return(pam_result$clustering)
}

# Main ensemble clustering function
ensemble_clustering <- function(
    X, 
    W, 
    k_values, 
    n_clusters,
    lambda1 = 1.0,
    clustering_methods = c("kmeans", "pam"),
    max_iter = 100,
    tol = 1e-6,
    verbose = TRUE
) {
  cat("========== Starting ensemble clustering ==========\n")
  cat(sprintf("Number of reduction ranks: %d\n", length(k_values)))
  cat(sprintf("Clustering methods: %s\n", paste(clustering_methods, collapse = ", ")))
  cat(sprintf("Target number of clusters: %d\n", n_clusters))
  
  n <- ncol(X)
  all_S_matrices <- list()
  all_labels <- list()
  matrix_counter <- 1
  
  for (k in k_values) {
    if (verbose) {
      cat(sprintf("\nReduction rank k=%d: ", k))
    }
    
    # Run NMF to obtain V
    nmf_result <- regularized_nmf(
      X = X,
      W = W,
      k = k,
      lambda1 = lambda1,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose,
      check_interval = 100
    )
    
    V <- nmf_result$V
    
    for (method in clustering_methods) {
      seed_val <- 1234 + k * 10 + match(method, clustering_methods)
      cat(sprintf("%s ", method))
      
      if (method == "kmeans") {
        labels <- perform_kmeans(V, n_clusters, seed = seed_val)
      } else if (method == "pam") {
        labels <- perform_pam(V, n_clusters, seed = seed_val)
      }
      
      labels <- as.numeric(as.factor(labels))
      S <- compute_association_matrix(labels)
      
      all_S_matrices[[matrix_counter]] <- S
      all_labels[[matrix_counter]] <- labels
      matrix_counter <- matrix_counter + 1
    }
  }
  
  cat(sprintf("\nTotal number of clustering results: %d\n", length(all_S_matrices)))
  cat("Computing consensus matrix...\n")
  
  S_consensus <- matrix(0, n, n)
  for (S in all_S_matrices) {
    S_consensus <- S_consensus + S
  }
  S_consensus <- S_consensus / length(all_S_matrices)
  
  diag(S_consensus) <- 1
  S_consensus <- (S_consensus + t(S_consensus)) / 2
  
  return(list(
    consensus_matrix = S_consensus,
    all_labels = all_labels,
    all_S_matrices = all_S_matrices,
    n_clusters = n_clusters
  ))
}