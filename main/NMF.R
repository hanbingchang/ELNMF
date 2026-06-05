compute_objective <- function(
    X, U, V,
    L_w, 
    lambda1 = 1
) {
  if (nrow(U) != nrow(X)) stop("U rows must equal X rows")
  if (nrow(V) != ncol(X)) stop("V rows must equal X columns")
  if (ncol(U) != ncol(V)) stop("U and V must have same number of columns")
  
  # Reconstruction error ||X - U V^T||_F^2
  UVt <- U %*% t(V)
  recon_error <- sum((X - UVt)^2)
  
  # Graph regularization term: lambda1 * Tr(V^T L_w V)
  graph_reg1 <- 0
  if (lambda1 > 0) {
    LwV <- L_w %*% V
    graph_reg1 <- lambda1 * sum(V * LwV)
  }
  
  total_obj <- recon_error + graph_reg1 
  
  return(list(
    total = total_obj,
    reconstruction = recon_error,
    graph_w = graph_reg1
  ))
}

regularized_nmf <- function(X, W, 
                            k = 20, 
                            lambda1 = 53, 
                            max_iter = 100, 
                            tol = 1e-6,
                            verbose = TRUE,
                            check_interval = 100) {
  
  cat("Initializing regularized NMF parameters...\n")
  
  m <- nrow(X)
  n <- ncol(X)
  
  cat(sprintf("Data dimensions: X = %d x %d\n", m, n))
  cat(sprintf("Reduction rank: k = %d\n", k))
  
  if (nrow(W) != n || ncol(W) != n) {
    stop(sprintf("Adjacency matrix W dimensions (%d x %d) do not match X columns (%d)", 
                 nrow(W), ncol(W), n))
  }
  
  W_norm <- W / mean(W[W > 0])
  degree_vec <- rowSums(W_norm)
  D_w <- Diagonal(x = degree_vec)
  L_w <- D_w - W_norm
  W <- W_norm
  
  cat("Initializing U and V matrices...\n")
  set.seed(1234)
  U <- matrix(runif(m * k, min = 0, max = 1), nrow = m, ncol = k)
  V <- matrix(runif(n * k, min = 0, max = 1), nrow = n, ncol = k)
  
  cat(sprintf("Starting iterative optimization (max iter: %d)...\n", max_iter))
  
  obj_history <- numeric(max_iter)
  loss_components <- matrix(0, nrow = max_iter, ncol = 3)
  colnames(loss_components) <- c("total", "reconstruction", "graph_w")
  
  convergence_info <- list(
    converged = FALSE,
    iter = max_iter,
    max_change = NA,
    reason = "Reached max iterations"
  )
  
  for (iter in 1:max_iter) {
    U_old <- U
    V_old <- V
    
    # Update U
    XV <- as.matrix(X %*% V)
    VTV <- crossprod(V)
    UVTV <- U %*% VTV
    U <- U * (XV / (UVTV + 1e-10))
    U[U < 0] <- 0
    
    # Update V
    XTU <- crossprod(X, U)
    if (lambda1 > 0) {
      WV <- as.matrix(W %*% V)
    } else {
      WV <- 0
    }
    numerator_V <- XTU + lambda1 * WV 
    
    UTU <- crossprod(U)
    VUTU <- V %*% UTU
    if (lambda1 > 0) {
      D_w_V <- as.matrix(D_w %*% V)
    } else {
      D_w_V <- 0
    }
    denominator_V <- VUTU + lambda1 * D_w_V
    V <- V * (numerator_V / (denominator_V + 1e-10))
    V[V < 0] <- 0
    
    # Compute objective
    if (iter %% check_interval == 0 || iter == 1) {
      obj_result <- compute_objective(X, U, V, L_w, lambda1)
    }
    obj_history[iter] <- obj_result$total
    loss_components[iter, ] <- c(obj_result$total, obj_result$reconstruction, obj_result$graph_w)
    
    if (verbose && (iter %% check_interval == 0 || iter == 1)) {
      cat(sprintf("Iter %3d: objective = %.4e", iter, obj_result$total))
      if (iter %% 10 == 0) {
        cat(sprintf(" (recon:%.2e, graph:%.2e)", obj_result$reconstruction, obj_result$graph_w))
      }
      cat("\n")
    }
    
    # Convergence check
    if (iter %% check_interval == 0 && iter > 10) {
      U_change <- norm(U - U_old, "F") / (norm(U_old, "F") + 1e-10)
      V_change <- norm(V - V_old, "F") / (norm(V_old, "F") + 1e-10)
      max_change <- max(U_change, V_change)
      if (max_change < tol) {
        message(paste("Algorithm converged at iteration:", iter, "max_change =", max_change))
        convergence_info$converged <- TRUE
        convergence_info$iter <- iter
        convergence_info$max_change <- max_change
        convergence_info$reason <- "Convergence tolerance reached"
        break
      }
      if (verbose) {
        message(paste("Iter", iter, "| max_change =", max_change, "| tol =", tol))
      }
    }
  }
  
  if (!convergence_info$converged) {
    message(paste("Reached max iterations", max_iter, "without convergence | final max_change =", max_change))
  }
  
  obj_history <- obj_history[1:convergence_info$iter]
  loss_components <- loss_components[1:convergence_info$iter, ]
  
  final_obj <- compute_objective(X, U, V, L_w, lambda1)
  
  result <- list(
    U = U,
    V = V,
    D_w = D_w,
    L_w = L_w,
    objective_history = obj_history,
    loss_components = loss_components,
    final_objective = final_obj,
    convergence = convergence_info,
    params = list(
      k = k,
      lambda1 = lambda1,
      max_iter = max_iter,
      tol = tol
    ),
    dimensions = list(
      m = m,
      n = n,
      k = k
    )
  )
  
  class(result) <- "regularized_nmf"
  
  if (verbose) {
    cat("\n========== NMF optimization complete ==========\n")
    cat(sprintf("Iterations: %d\n", convergence_info$iter))
    cat(sprintf("Converged: %s\n", ifelse(convergence_info$converged, "Yes", "No")))
    if (convergence_info$converged) {
      cat(sprintf("Reason: %s\n", convergence_info$reason))
    }
    cat(sprintf("Final objective value: %.4e\n", final_obj$total))
    cat("Loss components:\n")
    cat(sprintf("  Reconstruction: %.4e\n", final_obj$reconstruction))
    cat(sprintf("  Graph regularization: %.4e\n", final_obj$graph_w))
    cat("==============================================\n")
  }
  
  return(result)
}