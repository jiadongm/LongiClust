#' Simulate Multivariate Functional Data with Two Clusters
#'
#' @description
#' Generates synthetic multivariate functional data for two clusters with
#' controlled separation and within-group variation. This function is designed
#' for testing clustering algorithms on longitudinal omics data, where data
#' objects are vector-valued functions defined on a common domain.
#'
#' @param n1 Integer. Number of observations in cluster 1 (default: 50)
#' @param n2 Integer. Number of observations in cluster 2 (default: 50)
#' @param within_var Numeric. Controls within-group variation via eigenvalue
#'   scaling. Higher values increase variability within each cluster (default: 1.0)
#' @param between_sep Numeric. Controls separation between cluster means.
#'   Higher values increase distance between cluster centroids (default: 1.0)
#' @param start1_mean Numeric. Mean starting level for cluster 1 (default: 0)
#' @param start1_sd Numeric. Standard deviation of starting levels in cluster 1 (default: 0.5)
#' @param end1_mean Numeric. Mean ending level for cluster 1 (default: 0)
#' @param end1_sd Numeric. Standard deviation of ending levels in cluster 1 (default: 0.5)
#' @param start2_mean Numeric. Mean starting level for cluster 2 (default: 0)
#' @param start2_sd Numeric. Standard deviation of starting levels in cluster 2 (default: 0.5)
#' @param end2_mean Numeric. Mean ending level for cluster 2 (default: 0)
#' @param end2_sd Numeric. Standard deviation of ending levels in cluster 2 (default: 0.5)
#' @param n_basis Integer. Number of basis functions for functional expansion (default: 50)
#' @param n_grid Integer. Number of grid points for discretization (default: 50)
#'
#' @details
#' The function generates functional data using a basis expansion approach:
#'
#' **Data Generation Model:**
#' Each functional observation X_i(t) is generated as:
#' \deqn{X_i(t) = \sum_{j=1}^{n_{basis}} (\sqrt{\theta_j} Z_{ij} + \mu_{kj}) \phi_j(t) + L_i(t)}
#'
#' where:
#' \itemize{
#'   \item \eqn{\phi_j(t) = \sqrt{2} \sin(\pi j t)} are orthonormal basis functions
#'   \item \eqn{\theta_j = within_var \cdot j^{-2}} are eigenvalues controlling variance
#'   \item \eqn{Z_{ij} \sim N(0,1)} are random coefficients
#'   \item \eqn{\mu_{kj}} are cluster-specific mean coefficients
#'   \item \eqn{L_i(t)} is a linear trend from random start to end level
#' }
#'
#' **Cluster Separation:**
#' The between_sep parameter scales the difference between cluster means:
#' \deqn{\mu_k = \bar{\mu} + between_sep \cdot (\mu_k^{base} - \bar{\mu})}
#'
#' where \eqn{\mu_k^{base}} are predefined cluster-specific coefficients and
#' \eqn{\bar{\mu}} is their average.
#'
#' **Within-Cluster Variation:**
#' The within_var parameter controls the eigenvalues, affecting how much
#' individual curves deviate from their cluster mean. Individual-specific
#' linear trends add additional within-cluster heterogeneity.
#'
#' **Use Cases:**
#' This simulation is particularly useful for:
#' \itemize{
#'   \item Testing clustering algorithm performance under various difficulty levels
#'   \item Evaluating sensitivity to within-cluster variation
#'   \item Assessing robustness to cluster overlap
#'   \item Mimicking longitudinal gene expression or biomarker trajectories
#' }
#'
#' @return A list with the following components:
#'   \item{curves_pop1}{Matrix (n1 x n_grid) of functional observations for cluster 1}
#'   \item{curves_pop2}{Matrix (n2 x n_grid) of functional observations for cluster 2}
#'   \item{all_curves}{Matrix ((n1+n2) x n_grid) of all functional observations}
#'   \item{true_labels}{Vector of length (n1+n2) with true cluster assignments}
#'   \item{t_grid}{Vector of length n_grid with evaluation points}
#'   \item{mean1}{Vector of length n_grid with cluster 1 mean function}
#'   \item{mean2}{Vector of length n_grid with cluster 2 mean function}
#'   \item{within_sd1}{Numeric. Average standard deviation within cluster 1}
#'   \item{within_sd2}{Numeric. Average standard deviation within cluster 2}
#'   \item{mean_diff}{Numeric. RMSE of difference between cluster means}
#'   \item{separation_ratio}{Numeric. Ratio of between-cluster distance to
#'     average within-cluster standard deviation. Higher values indicate
#'     better-separated clusters}
#'   \item{params}{List of all input parameters}
#'
#' @examples
#' \dontrun{
#' # Scenario 1: Well-separated clusters with low within-group variation
#' # (Easy clustering problem)
#' data1 <- simulate_functional_data(
#'   n1 = 50, n2 = 50,
#'   within_var = 0.3,
#'   between_sep = 2.5,
#'   start1_mean = 2, start1_sd = 0.3,
#'   end1_mean = 0, end1_sd = 0.3,
#'   start2_mean = 0, start2_sd = 0.3,
#'   end2_mean = 4, end2_sd = 0.3
#' )
#' cat("Separation ratio:", data1$separation_ratio, "\n")
#'
#' # Scenario 2: High within-group variation with good separation
#' # (Moderate difficulty)
#' data2 <- simulate_functional_data(
#'   within_var = 2.0,
#'   between_sep = 2.5,
#'   start1_mean = 2, start1_sd = 0.5,
#'   end1_mean = 0, end1_sd = 0.5,
#'   start2_mean = 0, start2_sd = 0.5,
#'   end2_mean = 4, end2_sd = 0.5
#' )
#'
#' # Scenario 3: Overlapping clusters with low within-group variation
#' # (Difficult due to overlap)
#' data3 <- simulate_functional_data(
#'   within_var = 0.3,
#'   between_sep = 1.0,
#'   start1_mean = 2, start1_sd = 0.3,
#'   end1_mean = 0, end1_sd = 0.3,
#'   start2_mean = 0, start2_sd = 0.3,
#'   end2_mean = 1, end2_sd = 0.3
#' )
#'
#' # Scenario 4: High within-variation and low separation
#' # (Very difficult clustering problem)
#' data4 <- simulate_functional_data(
#'   within_var = 2.0,
#'   between_sep = 1.0
#' )
#'
#' # Multivariate functional data (d0 = 3 variables)
#' data_mv <- simulate_multivariate_functional_data(
#'   d0 = 3,
#'   n1 = 40, n2 = 60,
#'   within_var = 1.0,
#'   between_sep = 1.5
#' )
#' }
#'
#' @seealso
#' \code{\link{plot_functional_data}} for visualization
#' \code{\link{simulate_multivariate_functional_data}} for multivariate version
#'
#' @references
#' Delaigle, A., Hall, P., & Bathia, N. (2012). Componentwise classification
#' and clustering of functional data. Biometrika, 99(2), 299-313.
#'
#' Jacques, J., & Preda, C. (2014). Model-based clustering for multivariate
#' functional data. Computational Statistics & Data Analysis, 71, 92-106.
#'
#' @export
simulate_functional_data <- function(n1 = 50, n2 = 50,
                                     within_var = 1.0,      # Controls within-group variation
                                     between_sep = 1.0,     # Controls between-group separation
                                     start1_mean = 0, start1_sd = 0.5,  # Pop 1 starting level
                                     end1_mean = 0, end1_sd = 0.5,      # Pop 1 ending level
                                     start2_mean = 0, start2_sd = 0.5,  # Pop 2 starting level
                                     end2_mean = 0, end2_sd = 0.5,      # Pop 2 ending level
                                     n_basis = 50,
                                     n_grid = 50) {

  # Time grid
  t_grid <- seq(0, 1, length.out = n_grid)

  # Basis functions
  phi <- function(j, t) {
    sqrt(2) * sin(pi * j * t)
  }

  # Eigenvalues (controlling variance) - scaled by within_var
  theta <- sapply(1:n_basis, function(j) within_var * j^(-2))

  # Base mean coefficients (scaled by between_sep for group separation)
  mu1_base <- rep(0, n_basis)
  mu1_base[1:6] <- c(0, -0.30, 0.60, -0.30, 0.60, -0.30)

  mu2_base <- rep(0, n_basis)
  mu2_base[1:6] <- c(0, 0.45, 0.45, 0.09, 0.84, 0.60)

  # Scale the difference between groups
  mu_mean <- (mu1_base + mu2_base) / 2
  mu1 <- mu_mean + between_sep * (mu1_base - mu_mean)
  mu2 <- mu_mean + between_sep * (mu2_base - mu_mean)

  # Function to generate curves with group-specific random start/end levels
  generate_curves <- function(n_curves, mu_coef, theta, n_basis, n_grid, t_grid,
                              start_mean, start_sd, end_mean, end_sd) {
    curves <- matrix(0, nrow = n_curves, ncol = n_grid)

    for (i in 1:n_curves) {
      # Generate random coefficients Z_ij ~ N(0, 1)
      Z <- rnorm(n_basis)

      # Compute the base curve
      curve <- rep(0, n_grid)
      for (j in 1:n_basis) {
        coef <- sqrt(theta[j]) * Z[j] + mu_coef[j]
        curve <- curve + coef * phi(j, t_grid)
      }

      # Random starting and ending levels for this individual (group-specific)
      start_level <- rnorm(1, mean = start_mean, sd = start_sd)
      end_level <- rnorm(1, mean = end_mean, sd = end_sd)

      # Create linear trend from start to end
      linear_trend <- seq(start_level, end_level, length.out = n_grid)

      # Add the linear trend to the curve
      curves[i, ] <- curve + linear_trend
    }
    return(curves)
  }

  # Generate curves for each population with different start/end parameters
  curves_pop1 <- generate_curves(n1, mu1, theta, n_basis, n_grid, t_grid,
                                 start1_mean, start1_sd, end1_mean, end1_sd)
  curves_pop2 <- generate_curves(n2, mu2, theta, n_basis, n_grid, t_grid,
                                 start2_mean, start2_sd, end2_mean, end2_sd)

  # Calculate separation metrics
  mean1 <- colMeans(curves_pop1)
  mean2 <- colMeans(curves_pop2)
  mean_diff <- sqrt(mean((mean1 - mean2)^2))  # RMSE of mean difference

  within_sd1 <- mean(apply(curves_pop1, 2, sd))
  within_sd2 <- mean(apply(curves_pop2, 2, sd))
  avg_within_sd <- (within_sd1 + within_sd2) / 2

  separation_ratio <- mean_diff / avg_within_sd

  return(list(
    curves_pop1 = curves_pop1,
    curves_pop2 = curves_pop2,
    all_curves = rbind(curves_pop1, curves_pop2),
    true_labels = c(rep(1, n1), rep(2, n2)),
    t_grid = t_grid,
    mean1 = mean1,
    mean2 = mean2,
    within_sd1 = within_sd1,
    within_sd2 = within_sd2,
    mean_diff = mean_diff,
    separation_ratio = separation_ratio,
    params = list(within_var = within_var,
                  between_sep = between_sep,
                  start1_mean = start1_mean, start1_sd = start1_sd,
                  end1_mean = end1_mean, end1_sd = end1_sd,
                  start2_mean = start2_mean, start2_sd = start2_sd,
                  end2_mean = end2_mean, end2_sd = end2_sd,
                  n1 = n1, n2 = n2)
  ))
}

