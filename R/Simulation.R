# Simulate functional data with controlled variation and group separation


# Master function to control simulation parameters
#' Title
#'
#' @param n1
#' @param n2
#' @param within_var
#' @param between_sep
#' @param start1_mean
#' @param start1_sd
#' @param end1_mean
#' @param end1_sd
#' @param start2_mean
#' @param start2_sd
#' @param end2_mean
#' @param end2_sd
#' @param n_basis
#' @param n_grid
#'
#' @returns
#' @export
#'
#' @examples
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

