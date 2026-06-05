ensemble_clustering_filtered <- function(
    X,
    V_dir,
    k_values,
    lambda_values,
    n_clusters,
    clustering_methods = c("kmeans", "pam"),
    q = 0.75,                 # quantile threshold
    verbose = TRUE
) {
  set.seed(1234)
  n <- ncol(X)
  
  all_results <- list()      # store S + score
  evaluation_df <- data.frame()
  
  counter <- 1
  
  # ================= Phase 1: compute all =================
  for (lambda1 in lambda_values) {
    for (k in k_values) {
      
      if (verbose) {
        cat(sprintf("lambda=%d | k=%d\n", lambda1, k))
      }
      
      v_name <- sprintf("V_lambda%d_k%d", lambda1, k)
      v_path <- file.path(V_dir, paste0(v_name, ".rds"))
      
      if (!file.exists(v_path)) {
        warning(sprintf("Missing %s, skipping", v_name))
        next
      }
      
      V <- readRDS(v_path)
      
      # ---- Unify V orientation: spots x k ----
      if (nrow(V) == n) {
        # OK
      } else if (ncol(V) == n) {
        V <- t(V)
      } else {
        stop(sprintf(
          "Invalid V dimensions: %s, dim(V) = %d x %d, n = %d",
          v_name, nrow(V), ncol(V), n
        ))
      }
      
      for (method in clustering_methods) {
        
        labels <- if (method == "kmeans") {
          perform_kmeans(V, n_clusters)
        } else {
          perform_pam(V, n_clusters)
        }
        
        labels <- as.numeric(as.factor(labels))
        
        # ---- Compute cluster quality score ----
        score_res <- compute_cluster_score(V, labels)
        
        evaluation_df <- rbind(
          evaluation_df,
          data.frame(
            lambda = lambda1,
            k = k,
            method = method,
            silhouette = score_res$silhouette,
            CHAOS = score_res$chaos,
            score = score_res$score
          )
        )
        
        all_results[[counter]] <- list(
          S = compute_association_matrix(labels),
          score = score_res$score
        )
        
        counter <- counter + 1
      }
    }
  }
  
  # ================= Phase 2: dynamic quantile filtering =================
  scores <- evaluation_df$score
  scores <- scores[!is.na(scores)]
  
  if (length(scores) == 0) {
    stop("All clustering scores are NA")
  }
  
  thr <- as.numeric(quantile(scores, q))
  
  if (verbose) {
    cat(sprintf("\nUsing %.0f%% quantile as threshold: %.4f\n", q * 100, thr))
  }
  
  kept_S <- list()
  kept_scores <- c()
  
  for (res in all_results) {
    if (!is.na(res$score) && res$score >= thr) {
      kept_S[[length(kept_S) + 1]] <- res$S
      kept_scores <- c(kept_scores, res$score)
    }
  }
  
  # ================= Consensus matrix =================
  if (length(kept_S) == 0) {
    stop("No clustering results passed the quantile filter")
  }
  
  S_consensus <- Reduce("+", kept_S) / length(kept_S)
  diag(S_consensus) <- 1
  S_consensus <- (S_consensus + t(S_consensus)) / 2
  
  # ================= Statistics before and after filtering =================
  summary_before <- data.frame(
    silhouette_mean = mean(evaluation_df$silhouette, na.rm = TRUE),
    silhouette_median = median(evaluation_df$silhouette, na.rm = TRUE),
    CHAOS_mean = mean(evaluation_df$CHAOS, na.rm = TRUE),
    CHAOS_median = median(evaluation_df$CHAOS, na.rm = TRUE),
    score_mean = mean(evaluation_df$score, na.rm = TRUE),
    score_median = median(evaluation_df$score, na.rm = TRUE)
  )
  
  eval_after <- evaluation_df[evaluation_df$score >= thr, ]
  
  summary_after <- data.frame(
    silhouette_mean = mean(eval_after$silhouette, na.rm = TRUE),
    silhouette_median = median(eval_after$silhouette, na.rm = TRUE),
    CHAOS_mean = mean(eval_after$CHAOS, na.rm = TRUE),
    CHAOS_median = median(eval_after$CHAOS, na.rm = TRUE),
    score_mean = mean(eval_after$score, na.rm = TRUE),
    score_median = median(eval_after$score, na.rm = TRUE)
  )
  
  if (verbose) {
    cat(sprintf("Before filtering: %d clusterings\n", nrow(evaluation_df)))
    cat(sprintf("After filtering: %d clusterings\n", length(kept_S)))
  }
  
  return(list(
    consensus_matrix = S_consensus,
    all_S_matrices = kept_S,
    evaluation = evaluation_df,
    threshold = thr,
    q = q,
    summary_before = summary_before,
    summary_after = summary_after,
    n_kept = length(kept_S)
  ))
}