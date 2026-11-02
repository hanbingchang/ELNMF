compute_CHAOS_V <- function(V, labels) {
  clusters <- unique(labels)
  chaos_vec <- c()
  
  for (cl in clusters) {
    idx <- which(labels == cl)
    if (length(idx) < 2) next
    
    V_cl <- V[idx, , drop = FALSE]
    center <- colMeans(V_cl)
    diff_mat <- sweep(V_cl, 2, center, FUN = "-")
    dist_to_center <- sqrt(rowSums(diff_mat^2))
    chaos_vec <- c(chaos_vec, mean(dist_to_center))
  }
  
  mean(chaos_vec)
}

compute_cluster_score <- function(V, labels) {
  sil <- NA_real_
  sil_obj <- try(silhouette(labels, dist(V)), silent = TRUE)
  if (!inherits(sil_obj, "try-error")) {
    sil <- mean(sil_obj[, 3])
  }
  
  chaos <- compute_CHAOS_V(V, labels)
  
  if (is.na(sil) || is.na(chaos) || chaos == 0) {
    score <- NA_real_
  } else {
    score <- sil / chaos
  }
  
  return(list(
    silhouette = as.numeric(sil),
    chaos = as.numeric(chaos),
    score = as.numeric(score)
  ))
}