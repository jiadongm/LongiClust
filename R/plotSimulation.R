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

