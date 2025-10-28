#' Visualize Simulated Functional Data
#'
#' @param data List. Output from simulate_functional_data()
#' @param show_metrics Logical. Whether to display separation metrics on plot (default: TRUE)
#' @param alpha Numeric. Transparency for individual curves (default: 0.2)
#' @param col_pop1 Character. Color for cluster 1 (default: "blue")
#' @param col_pop2 Character. Color for cluster 2 (default: "red")
#' @param lwd_mean Numeric. Line width for mean functions (default: 3)
#' @param ... Additional arguments passed to matplot()
#'
#' @details
#' Creates a visualization showing:
#' \itemize{
#'   \item Individual trajectories for each cluster (semi-transparent)
#'   \item Mean trajectories for each cluster (solid lines)
#'   \item Optional separation metrics overlay
#' }
#'
#' The separation ratio (mean difference / average within-cluster SD) provides
#' a standardized measure of clustering difficulty, similar to Cohen's d in
#' two-sample comparisons.
#'
#' @examples
#' \dontrun{
#' # Generate and visualize data
#' data <- simulate_functional_data(
#'   within_var = 1.0,
#'   between_sep = 2.0
#' )
#' plot_functional_data(data)
#'
#' # Multiple scenarios in one figure
#' par(mfrow = c(2, 2))
#'
#' data1 <- simulate_functional_data(within_var = 0.3, between_sep = 2.5)
#' plot_functional_data(data1, main = "Easy: Low var, high sep")
#'
#' data2 <- simulate_functional_data(within_var = 2.0, between_sep = 2.5)
#' plot_functional_data(data2, main = "Moderate: High var, high sep")
#'
#' data3 <- simulate_functional_data(within_var = 0.3, between_sep = 1.0)
#' plot_functional_data(data3, main = "Moderate: Low var, low sep")
#'
#' data4 <- simulate_functional_data(within_var = 2.0, between_sep = 1.0)
#' plot_functional_data(data4, main = "Hard: High var, low sep")
#' }
#'
#' @export
plot_functional_data <- function(data, show_metrics = TRUE) {
  #par(mar = c(4, 4, 4, 1))

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

