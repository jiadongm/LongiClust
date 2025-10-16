# Simulate functional data with controlled variation and group separation
set.seed(123)

# Master function to control simulation parameters
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

# Plotting function
plot_functional_data <- function(data, show_metrics = TRUE) {
  par(mar = c(4, 4, 4, 1))
  
  ylim_range <- range(data$all_curves)
  
  matplot(data$t_grid, t(data$curves_pop1), type = "l", lty = 1, 
          col = rgb(0, 0, 1, 0.2),
          xlab = "Time (t)", ylab = "Expression Level X(t)", 
          main = sprintf("Functional Data"),
                         #\n Pop1: Start=%.1f±%.1f, End=%.1f±%.1f | Pop2: Start=%.1f±%.1f, End=%.1f±%.1f",
                         #data$params$start1_mean, data$params$start1_sd,
                         #data$params$end1_mean, data$params$end1_sd,
                         #data$params$start2_mean, data$params$start2_sd,
                         #data$params$end2_mean, data$params$end2_sd),
          ylim = ylim_range, cex.main = 0.8)
  
  matlines(data$t_grid, t(data$curves_pop2), type = "l", lty = 1, 
           col = rgb(1, 0, 0, 0.2))
  
  lines(data$t_grid, data$mean1, col = "blue", lwd = 3)
  lines(data$t_grid, data$mean2, col = "red", lwd = 3)

  
  legend("topright", 
         legend = c("Pop 1", "Pop 2"), 
         col = c("blue", "red"), 
         lwd = c(2, 2),
         bty = "n", cex = 0.8)
  
  if (show_metrics) {
    metrics_text <- sprintf(
      "Separation Ratio: %.2f\nAvg Within SD: %.2f",
      data$separation_ratio, 
      (data$within_sd1 + data$within_sd2)/2
    )
    mtext(metrics_text, side = 3, line = -3, adj = 0.02, cex = 0.7, col = "gray30")
  }
}

# Compare different scenarios
par(mfrow = c(2, 2), mar = c(4, 4, 4, 1))

# Scenario 1: Low within variation & high separation between groups

data1 <- simulate_functional_data(
  within_var = 0.3, between_sep = 2.5,
  start1_mean = 2, start1_sd = 0.3,
  end1_mean = 0, end1_sd = 0.3,
  start2_mean = 0, start2_sd = 0.3,
  end2_mean = 4, end2_sd = 0.3
)
plot_functional_data(data1)


# Scenario 2: High within variation & high separation between groups

data2 <- simulate_functional_data(
  within_var = 2.0, between_sep = 2.5,
  start1_mean = 2, start1_sd = 0.5,
  end1_mean = 0, end1_sd = 0.5,
  start2_mean = 0, start2_sd = 0.5,
  end2_mean = 4, end2_sd = 0.5
)
plot_functional_data(data2)


# Scenario 3: Low within variation & low separation between groups

data3 <- simulate_functional_data(
  within_var = 0.3, between_sep = 1,
  start1_mean = 2, start1_sd = 0.3,
  end1_mean = 0, end1_sd = 0.3,
  start2_mean = 0, start2_sd = 0.3,
  end2_mean = 1, end2_sd = 0.3
)

plot_functional_data(data3)


# Scenario 4: High within variation & low separation between groups

data4 <- simulate_functional_data(
  within_var = 2.0, between_sep =1,
  start1_mean = 2, start1_sd = 0.5,
  end1_mean = 0, end1_sd = 0.5,
  start2_mean = 0, start2_sd = 0.5,
  end2_mean = 1, end2_sd = 0.5
)

plot_functional_data(data4)
