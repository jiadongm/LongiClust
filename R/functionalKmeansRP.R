#' Functional K-Means Clustering with Random Projections
#'
#' @description
#' Implements the functional K-means algorithm using random projections to find
#' the optimal projection for clustering multivariate functional data.
#'
#' This function follows the algorithm from the manuscript:
#' 1. For each random projection Ψ*_b in Q_B:
#'    - Project the multivariate functional data onto Ψ*_b
#'    - Perform k-means clustering on the projected data
#'    - Compute the tightness criterion T̂(Ψ*_b)
#' 2. Select the projection Ψ̂_QB that minimizes T̂(Ψ)
#'
#' @param X_list List of length d0. Each element is an n x grid_length matrix
#'        representing one dimension of the multivariate functional data
#' @param Q_B Object of class 'random_projections' from generateRP()
#' @param k_clusters Integer. Number of clusters (default: 2)
#' @param nstart Integer. Number of random starts for k-means (default: 25)
#' @param iter_max Integer. Maximum iterations for k-means (default: 100)
#' @param verbose Logical. Print progress messages (default: FALSE)
#'
#' @return List containing:
#'   \item{optimal_projection}{Index of the optimal projection in Q_B}
#'   \item{optimal_Psi}{The optimal projection object}
#'   \item{cluster_labels}{Vector of cluster assignments (1 to k)}
#'   \item{projected_data}{n x d0 matrix of data projected onto optimal Psi}
#'   \item{tightness_values}{Vector of tightness values for each projection}
#'   \item{min_tightness}{Minimum tightness value achieved}
#'   \item{centers}{Cluster centers in projected space}
#'   \item{within_ss}{Within-cluster sum of squares for each cluster}
#'   \item{total_within_ss}{Total within-cluster sum of squares}
#'
#' @details
#' The tightness criterion is defined as:
#' \deqn{T̂(Î_1, Î_2 | Ψ) = \sum_{j=1}^{d_0} \frac{1}{\hat{\sigma}^2(\psi_j)}
#'       \frac{1}{n} \sum_{k=1}^2 \sum_{i \in Î_k} [X_{i,j}(\psi_j) - \bar{X}_{k,j}(\psi_j)]^2}
#'
#' where \eqn{\hat{\sigma}^2(\psi_j)} is the sample variance of the j-th projected dimension.
#'
#' @examples
#' \dontrun{
#' # Generate simulation data
#' data1 <- simulate_functional_data(
#'   within_var = 0.3, between_sep = 2.5,
#'   start1_mean = 2, start1_sd = 0.3,
#'   end1_mean = 0, end1_sd = 0.3,
#'   start2_mean = 0, start2_sd = 0.3,
#'   end2_mean = 4, end2_sd = 0.3
#' )
#'
#' # For univariate functional data (d0 = 1)
#' X_list <- list(data1$all_curves)
#'
#' # Generate random projections
#' Q_B <- generateRP(B = 100, d0 = 1, r = 10, basis = "bspline")
#'
#' # Run functional K-means
#' result <- functionalKmeansRP(X_list, Q_B, k_clusters = 2)
#'
#' # Check clustering accuracy
#' table(Predicted = result$cluster_labels, True = data1$true_labels)
#' }
#'
#' @export
functionalKmeansRP <- function(X_list,
                               Q_B,
                               k_clusters = 2,
                               nstart = 25,
                               iter_max = 100,
                               verbose = FALSE) {

  # Validate inputs
  if (!inherits(Q_B, "random_projections")) {
    stop("Q_B must be an object of class 'random_projections' from generateRP()")
  }

  d0 <- attr(Q_B, "d0")
  B <- attr(Q_B, "B")

  if (length(X_list) != d0) {
    stop(paste("X_list must have length d0 =", d0))
  }

  n <- nrow(X_list[[1]])

  # Validate that all dimensions have the same number of observations
  for (j in 1:d0) {
    if (nrow(X_list[[j]]) != n) {
      stop("All elements of X_list must have the same number of rows (observations)")
    }
  }

  if (verbose) {
    cat(sprintf("Starting Functional K-Means with Random Projections\n"))
    cat(sprintf("  Number of observations: %d\n", n))
    cat(sprintf("  Data dimension: %d\n", d0))
    cat(sprintf("  Number of projections to evaluate: %d\n", B))
    cat(sprintf("  Number of clusters: %d\n\n", k_clusters))
  }

  # Initialize storage for tightness values
  tightness_values <- numeric(B)
  cluster_assignments <- vector("list", B)
  projected_data_list <- vector("list", B)
  kmeans_results <- vector("list", B)

  # Loop through each projection in Q_B
  for (b in 1:B) {

    if (verbose && b %% 10 == 0) {
      cat(sprintf("Processing projection %d / %d\n", b, B))
    }

    # Get the b-th projection
    Psi_star_b <- Q_B[[b]]

    # Project the multivariate functional data onto Psi_star_b
    # Result is n x d0 matrix
    X_projected <- project_multivariate_data(X_list, Psi_star_b)

    # Store projected data
    projected_data_list[[b]] <- X_projected

    # Perform k-means clustering on the projected data
    kmeans_result <- stats::kmeans(X_projected,
                                   centers = k_clusters,
                                   nstart = nstart,
                                   iter.max = iter_max)

    # Store k-means results
    kmeans_results[[b]] <- kmeans_result
    cluster_assignments[[b]] <- kmeans_result$cluster

    # Compute tightness criterion T̂(Î_1, Î_2 | Ψ*_b)
    tightness <- compute_tightness(X_projected,
                                   kmeans_result$cluster,
                                   k_clusters)

    tightness_values[b] <- tightness
  }

  # Find the optimal projection (minimum tightness)
  optimal_idx <- which.min(tightness_values)
  min_tightness <- tightness_values[optimal_idx]

  if (verbose) {
    cat(sprintf("\nOptimal projection found: #%d\n", optimal_idx))
    cat(sprintf("Minimum tightness: %.6f\n", min_tightness))
  }

  # Extract results for optimal projection
  optimal_Psi <- Q_B[[optimal_idx]]
  optimal_clusters <- cluster_assignments[[optimal_idx]]
  optimal_projected_data <- projected_data_list[[optimal_idx]]
  optimal_kmeans <- kmeans_results[[optimal_idx]]

  # Create return object
  result <- list(
    optimal_projection = optimal_idx,
    optimal_Psi = optimal_Psi,
    cluster_labels = optimal_clusters,
    projected_data = optimal_projected_data,
    tightness_values = tightness_values,
    min_tightness = min_tightness,
    centers = optimal_kmeans$centers,
    within_ss = optimal_kmeans$withinss,
    total_within_ss = optimal_kmeans$tot.withinss,
    kmeans_result = optimal_kmeans,
    Q_B = Q_B,
    k_clusters = k_clusters,
    n = n,
    d0 = d0,
    B = B
  )

  class(result) <- c("functionalKmeansRP", "list")

  return(result)
}


#' Compute Tightness Criterion for Functional K-Means
#'
#' @description
#' Computes the weighted tightness criterion from equation (4.2) in the manuscript:
#' \deqn{T̂(Î_1, Î_2 | Ψ) = \sum_{j=1}^{d_0} \frac{1}{\hat{\sigma}^2(\psi_j)}
#'       \frac{1}{n} \sum_{k=1}^K \sum_{i \in Î_k} [X_{i,j}(\psi_j) - \bar{X}_{k,j}(\psi_j)]^2}
#'
#' @param X_projected Matrix. n x d0 matrix of projected data
#' @param cluster_labels Vector. Cluster assignments (length n)
#' @param k_clusters Integer. Number of clusters
#'
#' @return Numeric. The tightness value
#'
#' @details
#' The tightness criterion standardizes each dimension by its variance before
#' computing the within-cluster sum of squares. This prevents dimensions with
#' large variance from dominating the clustering.
#'
#' @keywords internal
compute_tightness <- function(X_projected, cluster_labels, k_clusters) {

  n <- nrow(X_projected)
  d0 <- ncol(X_projected)

  # Compute overall variance for each dimension (standardization)
  sigma_sq <- apply(X_projected, 2, stats::var)

  # Avoid division by zero
  sigma_sq[sigma_sq < 1e-10] <- 1e-10

  # Initialize tightness
  tightness <- 0

  # Loop over each dimension
  for (j in 1:d0) {

    # Loop over each cluster
    dimension_tightness <- 0

    for (k in 1:k_clusters) {

      # Get indices of observations in cluster k
      cluster_k_idx <- which(cluster_labels == k)
      nk <- length(cluster_k_idx)

      if (nk > 0) {
        # Get data for cluster k, dimension j
        X_kj <- X_projected[cluster_k_idx, j]

        # Compute cluster mean
        X_bar_kj <- mean(X_kj)

        # Compute within-cluster sum of squares
        ss_kj <- sum((X_kj - X_bar_kj)^2)

        dimension_tightness <- dimension_tightness + ss_kj
      }
    }

    # Add weighted contribution from dimension j
    tightness <- tightness + (1 / sigma_sq[j]) * (1 / n) * dimension_tightness
  }

  return(tightness)
}


#' Print Method for functionalKmeansRP
#'
#' @param x Object of class functionalKmeansRP
#' @param ... Additional arguments (unused)
#'
#' @export
print.functionalKmeansRP <- function(x, ...) {
  cat("Functional K-Means Clustering with Random Projections\n")
  cat("======================================================\n\n")

  cat("Clustering Summary:\n")
  cat(sprintf("  Number of observations: %d\n", x$n))
  cat(sprintf("  Data dimension (d0): %d\n", x$d0))
  cat(sprintf("  Number of clusters: %d\n", x$k_clusters))
  cat(sprintf("  Number of projections evaluated: %d\n\n", x$B))

  cat("Optimal Projection:\n")
  cat(sprintf("  Selected projection: #%d\n", x$optimal_projection))
  cat(sprintf("  Minimum tightness: %.6f\n\n", x$min_tightness))

  cat("Cluster Sizes:\n")
  cluster_table <- table(x$cluster_labels)
  for (k in 1:x$k_clusters) {
    cat(sprintf("  Cluster %d: %d observations\n", k, cluster_table[k]))
  }

  cat(sprintf("\nTotal within-cluster sum of squares: %.4f\n", x$total_within_ss))

  cat("\nUse summary() for more details.\n")
  cat("Use plot() to visualize clustering results.\n")
}


#' Summary Method for functionalKmeansRP
#'
#' @param object Object of class functionalKmeansRP
#' @param ... Additional arguments (unused)
#'
#' @export
summary.functionalKmeansRP <- function(object, ...) {

  cat("Functional K-Means Clustering with Random Projections\n")
  cat("======================================================\n\n")

  cat("Data Information:\n")
  cat(sprintf("  Number of observations: %d\n", object$n))
  cat(sprintf("  Data dimension (d0): %d\n", object$d0))
  cat(sprintf("  Number of clusters: %d\n\n", object$k_clusters))

  cat("Algorithm Configuration:\n")
  cat(sprintf("  Number of random projections (B): %d\n", object$B))
  cat(sprintf("  Basis type: %s\n", attr(object$Q_B, "basis_type")))
  cat(sprintf("  Basis truncation (r): %d\n\n", attr(object$Q_B, "r")))

  cat("Optimal Projection:\n")
  cat(sprintf("  Selected projection: #%d out of %d\n",
              object$optimal_projection, object$B))
  cat(sprintf("  Minimum tightness: %.6f\n", object$min_tightness))
  cat(sprintf("  Mean tightness across all projections: %.6f\n",
              mean(object$tightness_values)))
  cat(sprintf("  SD of tightness values: %.6f\n\n",
              sd(object$tightness_values)))

  cat("Clustering Results:\n")
  cluster_table <- table(object$cluster_labels)
  for (k in 1:object$k_clusters) {
    cat(sprintf("  Cluster %d:\n", k))
    cat(sprintf("    Size: %d observations (%.1f%%)\n",
                cluster_table[k], 100 * cluster_table[k] / object$n))
    cat(sprintf("    Within-cluster SS: %.4f\n", object$within_ss[k]))
  }

  cat(sprintf("\nTotal within-cluster sum of squares: %.4f\n",
              object$total_within_ss))
  cat(sprintf("Between-cluster sum of squares: %.4f\n",
              object$kmeans_result$betweenss))

  # Silhouette coefficient if available
  if (requireNamespace("cluster", quietly = TRUE)) {
    sil <- cluster::silhouette(object$cluster_labels,
                               stats::dist(object$projected_data))
    avg_sil <- mean(sil[, "sil_width"])
    cat(sprintf("Average silhouette width: %.4f\n", avg_sil))
  }

  cat("\n")
  invisible(object)
}


#' Plot Method for functionalKmeansRP
#'
#' @param x Object of class functionalKmeansRP
#' @param type Character. Type of plot: "tightness", "clusters", "projection", or "all"
#' @param ... Additional arguments passed to plotting functions
#'
#' @export
plot.functionalKmeansRP <- function(x, type = "all", ...) {

  type <- match.arg(type, c("all", "tightness", "clusters", "projection"))

  if (type == "all") {
    # Create a multi-panel plot
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    if (x$d0 == 1) {
      par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
    } else if (x$d0 == 2) {
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    } else {
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    }

    plot_tightness(x, ...)
    plot_clusters(x, ...)

    if (x$d0 >= 2) {
      plot_projection_2d(x, dims = c(1, 2), ...)
    }

  } else if (type == "tightness") {
    plot_tightness(x, ...)
  } else if (type == "clusters") {
    plot_clusters(x, ...)
  } else if (type == "projection") {
    if (x$d0 == 1) {
      plot_clusters(x, ...)
    } else {
      plot_projection_2d(x, dims = c(1, 2), ...)
    }
  }

  invisible(x)
}


#' Plot Tightness Values
#'
#' @keywords internal
plot_tightness <- function(x, ...) {
  plot(1:x$B, x$tightness_values,
       type = "p", pch = 20, col = "gray60",
       xlab = "Projection Index",
       ylab = "Tightness",
       main = "Tightness Values for Each Projection",
       ...)

  # Highlight optimal projection
  points(x$optimal_projection, x$min_tightness,
         pch = 19, col = "red", cex = 1.5)

  # Add horizontal line for mean
  abline(h = mean(x$tightness_values), col = "blue", lty = 2)

  legend("topright",
         legend = c("Optimal", "Mean"),
         col = c("red", "blue"),
         pch = c(19, NA),
         lty = c(NA, 2),
         bty = "n")

  grid(col = "lightgray", lty = "dotted")
}


#' Plot Clusters in Projected Space
#'
#' @keywords internal
plot_clusters <- function(x, ...) {

  if (x$d0 == 1) {
    # 1D case: histogram or strip chart
    plot(x$projected_data[, 1], rep(0, x$n),
         col = x$cluster_labels + 1,
         pch = 19,
         xlab = "Projected Data (Dimension 1)",
         ylab = "",
         main = "Cluster Assignments (1D Projection)",
         yaxt = "n",
         ...)

    # Add cluster centers
    for (k in 1:x$k_clusters) {
      abline(v = x$centers[k], col = k + 1, lwd = 2, lty = 2)
    }

    legend("topright",
           legend = paste("Cluster", 1:x$k_clusters),
           col = (1:x$k_clusters) + 1,
           pch = 19,
           bty = "n")

  } else {
    # Multidimensional case: scatterplot of first two dimensions
    plot_projection_2d(x, dims = c(1, 2), ...)
  }
}


#' Plot 2D Projection
#'
#' @keywords internal
plot_projection_2d <- function(x, dims = c(1, 2), ...) {

  if (x$d0 < 2) {
    warning("Cannot create 2D plot with d0 < 2")
    return(invisible(NULL))
  }

  dim1 <- dims[1]
  dim2 <- dims[2]

  plot(x$projected_data[, dim1], x$projected_data[, dim2],
       col = x$cluster_labels + 1,
       pch = 19,
       xlab = sprintf("Projected Dimension %d", dim1),
       ylab = sprintf("Projected Dimension %d", dim2),
       main = "Cluster Assignments (2D Projection)",
       ...)

  # Add cluster centers
  points(x$centers[, dim1], x$centers[, dim2],
         col = 1:x$k_clusters + 1,
         pch = 4,
         cex = 2,
         lwd = 3)

  legend("topright",
         legend = c(paste("Cluster", 1:x$k_clusters), "Centers"),
         col = c((1:x$k_clusters) + 1, "black"),
         pch = c(rep(19, x$k_clusters), 4),
         bty = "n")

  grid(col = "lightgray", lty = "dotted")
}


#' Evaluate Clustering Accuracy
#'
#' @description
#' Computes clustering accuracy metrics when true labels are known
#'
#' @param predicted Vector. Predicted cluster labels
#' @param true Vector. True cluster labels
#'
#' @return List containing:
#'   \item{accuracy}{Overall accuracy (best permutation)}
#'   \item{confusion_matrix}{Confusion matrix}
#'   \item{adjusted_rand_index}{Adjusted Rand Index (if mclust available)}
#'
#' @export
evaluate_clustering <- function(predicted, true) {

  # Create confusion matrix
  conf_mat <- table(Predicted = predicted, True = true)

  # Find best permutation for accuracy
  k <- length(unique(true))
  perms <- gtools::permutations(k, k)

  best_accuracy <- 0
  best_perm <- NULL

  for (i in 1:nrow(perms)) {
    perm <- perms[i, ]
    mapped_pred <- perm[predicted]
    accuracy <- mean(mapped_pred == true)

    if (accuracy > best_accuracy) {
      best_accuracy <- accuracy
      best_perm <- perm
    }
  }

  # Adjusted Rand Index
  ari <- NULL
  if (requireNamespace("mclust", quietly = TRUE)) {
    ari <- mclust::adjustedRandIndex(predicted, true)
  }

  result <- list(
    accuracy = best_accuracy,
    confusion_matrix = conf_mat,
    adjusted_rand_index = ari,
    best_permutation = best_perm
  )

  return(result)
}
