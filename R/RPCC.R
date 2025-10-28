#' Random Projection Consensus Clustering for Multivariate Functional Data
#'
#' @param X_list List of length d0. Each element is an n x grid_length matrix
#'        of functional observations for one dimension
#' @param B1 Integer. Number of consensus iterations (default: 50)
#' @param B2 Integer. Number of random projections per iteration (default: 100)
#' @param r Integer. Number of basis functions for random projections
#' @param basis Character. Basis type: "bspline" (default) or "fourier"
#' @param grid Vector. Grid points where functions are evaluated
#' @param bspline_degree Integer. Degree of B-splines if basis = "bspline" (default 3)
#' @param bspline_knots Vector. Optional knot positions for B-splines
#' @param kmeans_nstart Integer. Number of random starts for k-means (default: 25)
#' @param kmeans_iter_max Integer. Max iterations for k-means (default: 100)
#' @param verbose Logical. Print progress messages (default: FALSE)
#'
#' @return List containing:
#'   \item{partition}{Vector of cluster labels (1 or 2) for each observation}
#'   \item{I1}{Indices of observations in cluster 1}
#'   \item{I2}{Indices of observations in cluster 2}
#'   \item{all_partitions}{Matrix B1 x n of all aligned cluster labels}
#'   \item{consensus_matrix}{n x n matrix of pairwise co-clustering frequencies}
#'   \item{B1}{Number of consensus iterations used}
#'   \item{B2}{Number of random projections per iteration}
#'
#' @details
#' This function implements the Random Projection Consensus Clustering (RPCC)
#' algorithm for multivariate functional data. The algorithm:
#' 1. Generates B1 different clustering results using random projections
#' 2. Aligns cluster labels across iterations using the Hungarian algorithm
#' 3. Aggregates aligned labels via mean partition (majority voting)
#'
#' @examples
#' \dontrun{
#'   # Simulate bivariate functional data from two populations
#'   n1 <- 50; n2 <- 50; n <- n1 + n2
#'   grid <- seq(0, 1, length.out = 100)
#'
#'   # Population 1: smooth increasing functions
#'   X1_dim1 <- matrix(rnorm(n1 * 100, mean = grid, sd = 0.1), nrow = n1)
#'   X1_dim2 <- matrix(rnorm(n1 * 100, mean = grid^2, sd = 0.1), nrow = n1)
#'
#'   # Population 2: smooth decreasing functions
#'   X2_dim1 <- matrix(rnorm(n2 * 100, mean = 1 - grid, sd = 0.1), nrow = n2)
#'   X2_dim2 <- matrix(rnorm(n2 * 100, mean = (1-grid)^2, sd = 0.1), nrow = n2)
#'
#'   # Combine into list format
#'   X_list <- list(
#'     rbind(X1_dim1, X2_dim1),
#'     rbind(X1_dim2, X2_dim2)
#'   )
#'
#'   # Apply RPCC
#'   result <- RPCC(X_list, B1 = 50, B2 = 100, r = 8,
#'                  basis = "bspline", grid = grid)
#'
#'   # Check clustering accuracy
#'   true_labels <- c(rep(1, n1), rep(2, n2))
#'   table(true_labels, result$partition)
#' }
#'
#' @export
RPCC <- function(X_list,
                 B1 = 50,
                 B2 = 100,
                 r,
                 basis = "bspline",
                 grid = seq(0, 1, length.out = 100),
                 bspline_degree = 3,
                 bspline_knots = NULL,
                 kmeans_nstart = 25,
                 kmeans_iter_max = 100,
                 verbose = FALSE) {

  # Validate inputs
  if (!is.list(X_list)) {
    stop("X_list must be a list of matrices")
  }

  d0 <- length(X_list)
  n <- nrow(X_list[[1]])

  # Check all dimensions have same number of observations
  for (j in 1:d0) {
    if (!is.matrix(X_list[[j]])) {
      stop(paste("X_list[[", j, "]] must be a matrix", sep = ""))
    }
    if (nrow(X_list[[j]]) != n) {
      stop("All matrices in X_list must have the same number of rows")
    }
  }

  if (B1 <= 0 || B2 <= 0 || r <= 0) {
    stop("B1, B2, and r must be positive integers")
  }

  # Store all cluster labels from B1 iterations
  all_labels <- matrix(0, nrow = B1, ncol = n)

  # Store optimal projections for inspection (optional)
  optimal_projections <- vector("list", B1)

  if (verbose) {
    cat(sprintf("Running RPCC with B1=%d consensus iterations, B2=%d projections per iteration\n",
                B1, B2))
    cat(sprintf("Data: n=%d observations, d0=%d dimensions, grid length=%d\n",
                n, d0, length(grid)))
  }

  # ============================================================================
  # STEP 1: GENERATE B1 CLUSTERING RESULTS
  # ============================================================================

  for (b1 in 1:B1) {
    if (verbose && b1 %% 10 == 0) {
      cat(sprintf("  Iteration %d/%d...\n", b1, B1))
    }

    # Generate B2 random projections
    Q_B2 <- generateRP(B = B2,
                       d0 = d0,
                       r = r,
                       basis = basis,
                       grid = grid,
                       bspline_degree = bspline_degree,
                       bspline_knots = bspline_knots)

    # Find optimal projection and cluster via functional k-means
    kmeans_result <- funKMeans(X_list = X_list,
                               Q_B = Q_B2,
                               kmeans_nstart = kmeans_nstart,
                               kmeans_iter_max = kmeans_iter_max)

    # Store cluster labels (1 or 2)
    all_labels[b1, ] <- kmeans_result$partition

    # Optionally store the optimal projection
    optimal_projections[[b1]] <- kmeans_result$optimal_projection
  }

  if (verbose) {
    cat("All clustering iterations complete.\n")
  }

  # ============================================================================
  # STEP 2: ALIGN CLUSTER LABELS USING SIMPLE VOTING METHOD
  # ============================================================================

  if (verbose) {
    cat("Aligning cluster labels...\n")
  }

  # Use first result as reference
  aligned_labels <- matrix(0, nrow = B1, ncol = n)
  aligned_labels[1, ] <- all_labels[1, ]

  # Align all other results to the reference using Hungarian algorithm
  for (b1 in 2:B1) {
    # Compute contingency matrix: omega[k, k'] = number of observations
    # labeled k in reference and k' in current iteration
    omega <- matrix(0, nrow = 2, ncol = 2)

    for (k in 1:2) {
      for (k_prime in 1:2) {
        omega[k, k_prime] <- sum((aligned_labels[1, ] == k) &
                                   (all_labels[b1, ] == k_prime))
      }
    }

    # Find optimal label permutation using Hungarian algorithm
    # We want to maximize sum of omega[k, Theta(k)]
    # This is equivalent to minimizing the negative
    assignment <- solve_assignment(-omega)

    # Relabel current iteration according to assignment
    new_labels <- all_labels[b1, ]
    for (old_label in 1:2) {
      new_label <- assignment[old_label]
      aligned_labels[b1, new_labels == old_label] <- new_label
    }
  }

  if (verbose) {
    cat("Label alignment complete.\n")
  }

  # ============================================================================
  # STEP 3: AGGREGATE VIA MEAN PARTITION (MAJORITY VOTING)
  # ============================================================================

  if (verbose) {
    cat("Computing consensus partition...\n")
  }

  # For each observation, calculate average aligned label
  avg_labels <- colMeans(aligned_labels)

  # Assign to cluster 2 if average > 1.5, otherwise cluster 1
  final_partition <- ifelse(avg_labels > 1.5, 2, 1)

  # Create index sets
  I1 <- which(final_partition == 1)
  I2 <- which(final_partition == 2)

  # Compute consensus matrix (co-clustering frequency)
  consensus_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      # Proportion of times i and j were in same cluster
      consensus_matrix[i, j] <- mean(aligned_labels[, i] == aligned_labels[, j])
    }
  }

  if (verbose) {
    cat(sprintf("Clustering complete. Cluster sizes: %d and %d\n",
                length(I1), length(I2)))
  }

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  result <- list(
    partition = final_partition,
    I1 = I1,
    I2 = I2,
    all_partitions = aligned_labels,
    consensus_matrix = consensus_matrix,
    avg_labels = avg_labels,
    optimal_projections = optimal_projections,
    B1 = B1,
    B2 = B2,
    r = r,
    d0 = d0,
    n = n
  )

  class(result) <- c("RPCC_result", "list")

  return(result)
}


#' Solve Assignment Problem (Hungarian Algorithm)
#'
#' @param cost_matrix Matrix. Cost matrix for assignment problem
#'
#' @return Vector of assignments
#'
#' @details
#' Solves the linear assignment problem using a simple implementation.
#' For a 2x2 matrix (our case), we can use brute force.
#' For larger problems, consider using the clue::solve_LSAP function.
#'
solve_assignment <- function(cost_matrix) {

  # For 2x2 case, we can use brute force
  # Two possible assignments: [1,2] or [2,1]

  if (nrow(cost_matrix) != 2 || ncol(cost_matrix) != 2) {
    stop("This implementation only handles 2x2 assignment problems")
  }

  # Option 1: 1->1, 2->2
  cost1 <- cost_matrix[1, 1] + cost_matrix[2, 2]

  # Option 2: 1->2, 2->1
  cost2 <- cost_matrix[1, 2] + cost_matrix[2, 1]

  if (cost1 <= cost2) {
    return(c(1, 2))
  } else {
    return(c(2, 1))
  }
}


#' Print method for RPCC results
#'
#' @param x Object of class RPCC_result
#' @param ... Additional arguments (unused)
#'
#' @export
print.RPCC_result <- function(x, ...) {
  cat("Random Projection Consensus Clustering Results\n")
  cat("==============================================\n\n")

  cat(sprintf("Number of observations: %d\n", x$n))
  cat(sprintf("Data dimension (d0): %d\n", x$d0))
  cat(sprintf("Consensus iterations (B1): %d\n", x$B1))
  cat(sprintf("Random projections per iteration (B2): %d\n", x$B2))
  cat(sprintf("Basis functions used (r): %d\n\n", x$r))

  cat("Cluster sizes:\n")
  cat(sprintf("  Cluster 1: %d observations (%.1f%%)\n",
              length(x$I1), 100 * length(x$I1) / x$n))
  cat(sprintf("  Cluster 2: %d observations (%.1f%%)\n",
              length(x$I2), 100 * length(x$I2) / x$n))

  cat("\nConsensus strength:\n")
  cat(sprintf("  Mean co-clustering frequency: %.3f\n",
              mean(x$consensus_matrix)))
  cat(sprintf("  Within-cluster consensus: %.3f\n",
              mean(c(x$consensus_matrix[x$I1, x$I1][lower.tri(x$consensus_matrix[x$I1, x$I1])],
                     x$consensus_matrix[x$I2, x$I2][lower.tri(x$consensus_matrix[x$I2, x$I2])]))))
  cat(sprintf("  Between-cluster consensus: %.3f\n",
              mean(x$consensus_matrix[x$I1, x$I2])))

  invisible(x)
}


#' Summary method for RPCC results
#'
#' @param object Object of class RPCC_result
#' @param ... Additional arguments (unused)
#'
#' @export
summary.RPCC_result <- function(object, ...) {
  print(object)

  cat("\nStability across iterations:\n")

  # Calculate adjusted Rand index between consecutive iterations
  if (requireNamespace("mclust", quietly = TRUE)) {
    ari_values <- numeric(object$B1 - 1)
    for (b in 1:(object$B1 - 1)) {
      ari_values[b] <- mclust::adjustedRandIndex(object$all_partitions[b, ],
                                                 object$all_partitions[b + 1, ])
    }
    cat(sprintf("  Mean ARI between consecutive iterations: %.3f\n",
                mean(ari_values)))
    cat(sprintf("  SD of ARI: %.3f\n", stats::sd(ari_values)))
  } else {
    cat("  (Install 'mclust' package for ARI calculation)\n")
  }

  # Identify observations with low consensus
  obs_consensus <- numeric(object$n)
  for (i in 1:object$n) {
    cluster <- object$partition[i]
    cluster_indices <- if (cluster == 1) object$I1 else object$I2
    obs_consensus[i] <- mean(object$consensus_matrix[i, cluster_indices])
  }

  low_consensus <- which(obs_consensus < 0.7)
  if (length(low_consensus) > 0) {
    cat(sprintf("\nObservations with low consensus (< 0.7): %d\n",
                length(low_consensus)))
    cat("  Indices:", head(low_consensus, 10))
    if (length(low_consensus) > 10) cat(" ...")
    cat("\n")
  }

  invisible(object)
}
