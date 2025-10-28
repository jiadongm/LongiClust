# Enhanced Visualization: Cluster Predicted Profiles with Correct Predictions
# =============================================================================

# Source files
source("R/Simulation.R")
source("R/plotSimulation.R")
source("R/generateRP.R")
source("R/functionalKmeansRP.R")

# Set seed for reproducibility
set.seed(123)
par(mfrow = c(1,1))

# Generate data
data1 <- simulate_functional_data(
  within_var = 2.0, between_sep = 2.5,
  start1_mean = 2, start1_sd = 0.5,
  end1_mean = 0, end1_sd = 0.5,
  start2_mean = 0, start2_sd = 0.5,
  end2_mean = 4, end2_sd = 0.5
)

cat("\n=== Functional K-Means Clustering Results ===\n\n")

# Show original data
plot_functional_data(data1)

# Cluster
X_list <- list(data1$all_curves)
Q_B <- generateRP(B = 100, d0 = 1, r = 10, basis = "bspline", grid = data1$t_grid)
result <- functionalKmeansRP(X_list, Q_B, k_clusters = 2, verbose = TRUE)

cat("\n")
print(result)

# Compute confusion matrix to determine label mapping
confusion <- table(Predicted = result$cluster_labels, True = data1$true_labels)
cat("\n\nConfusion Matrix:\n")
print(confusion)

# Determine if we need to flip labels
# Check which mapping gives higher accuracy
acc_direct <- (confusion[1,1] + confusion[2,2]) / sum(confusion)
acc_flipped <- (confusion[1,2] + confusion[2,1]) / sum(confusion)

if (acc_flipped > acc_direct) {
  # Flip predicted labels
  predicted_labels_aligned <- 3 - result$cluster_labels  # 1->2, 2->1
  cat("\nNote: Flipped cluster labels for alignment with true labels\n")
} else {
  predicted_labels_aligned <- result$cluster_labels
}

# Calculate accuracy
accuracy <- max(acc_direct, acc_flipped) * 100
cat(sprintf("\nClustering Accuracy: %.1f%%\n", accuracy))

# Identify correct and incorrect predictions
correct_predictions <- (predicted_labels_aligned == data1$true_labels)
n_correct <- sum(correct_predictions)
n_incorrect <- sum(!correct_predictions)

cat(sprintf("Correct predictions: %d / %d\n", n_correct, length(correct_predictions)))
cat(sprintf("Incorrect predictions: %d / %d\n\n", n_incorrect, length(correct_predictions)))


# =============================================================================
# ENHANCED VISUALIZATION
# =============================================================================

# Create comprehensive figure
par(mfrow = c(2, 3), mar = c(4, 4, 3, 2), oma = c(0, 0, 3, 0))

# Color schemes
col_true1 <- rgb(0, 0, 1, 0.6)      # Blue for true cluster 1
col_true2 <- rgb(1, 0, 0, 0.6)      # Red for true cluster 2
col_correct <- rgb(0, 0.8, 0, 0.8)  # Green for correct predictions
col_incorrect <- rgb(1, 0.5, 0, 0.8) # Orange for incorrect predictions

# ------------------------------------------------------------------
# PLOT 1: Original Data (True Labels)
# ------------------------------------------------------------------
plot(data1$t_grid, data1$all_curves[1, ], type = "n",
     ylim = range(data1$all_curves),
     xlab = "Time", ylab = "Value",
     main = "Original Data (True Labels)")

# Plot population 1
for (i in which(data1$true_labels == 1)) {
  lines(data1$t_grid, data1$all_curves[i, ], col = col_true1, lwd = 1)
}

# Plot population 2
for (i in which(data1$true_labels == 2)) {
  lines(data1$t_grid, data1$all_curves[i, ], col = col_true2, lwd = 1)
}

# Add mean curves
lines(data1$t_grid, data1$mean1, col = "black", lwd = 1)
lines(data1$t_grid, data1$mean2, col = "black", lwd = 2)

legend("topright",
       legend = c("True Pop 1", "True Pop 2", "Mean 1", "Mean 2"),
       col = c(col_true1, col_true2, "blue", "red"),
       lwd = c(1, 1, 3, 3),
       bty = "n", cex = 0.8)

grid(col = "lightgray", lty = "dotted")


# ------------------------------------------------------------------
# PLOT 2: Predicted Clusters (with correct/incorrect overlay)
# ------------------------------------------------------------------
plot(data1$t_grid, data1$all_curves[1, ], type = "n",
     ylim = range(data1$all_curves),
     xlab = "Time", ylab = "Value",
     main = "Predicted Clusters")

# Plot predicted cluster 1
for (i in which(predicted_labels_aligned == 1)) {
  if (correct_predictions[i]) {
    lines(data1$t_grid, data1$all_curves[i, ], col = col_true1, lwd = 1.5)
  } else {
    lines(data1$t_grid, data1$all_curves[i, ], col = col_incorrect, lwd = 1.5, lty = 2)
  }
}

# Plot predicted cluster 2
for (i in which(predicted_labels_aligned == 2)) {
  if (correct_predictions[i]) {
    lines(data1$t_grid, data1$all_curves[i, ], col = col_true2, lwd = 1.5)
  } else {
    lines(data1$t_grid, data1$all_curves[i, ], col = col_incorrect, lwd = 1.5, lty = 2)
  }
}

# Add mean curves for predicted clusters
pred_mean1 <- colMeans(data1$all_curves[predicted_labels_aligned == 1, , drop = FALSE])
pred_mean2 <- colMeans(data1$all_curves[predicted_labels_aligned == 2, , drop = FALSE])
lines(data1$t_grid, pred_mean1, col = "black", lwd = 1)
lines(data1$t_grid, pred_mean2, col = "black", lwd = 2)

legend("topright",
       legend = c("Pred Cluster 1", "Pred Cluster 2", "Misclassified"),
       col = c(col_true1, col_true2, col_incorrect),
       lwd = c(1.5, 1.5, 1.5),
       lty = c(1, 1, 2),
       bty = "n", cex = 0.8)

grid(col = "lightgray", lty = "dotted")


# ------------------------------------------------------------------
# PLOT 3: Correctly Predicted Only
# ------------------------------------------------------------------
plot(data1$t_grid, data1$all_curves[1, ], type = "n",
     ylim = range(data1$all_curves),
     xlab = "Time", ylab = "Value",
     main = sprintf("Correct Predictions (%d/%d = %.1f%%)",
                    n_correct, length(correct_predictions), accuracy))

# Plot only correctly predicted curves
correct_idx <- which(correct_predictions)

for (i in correct_idx) {
  if (data1$true_labels[i] == 1) {
    lines(data1$t_grid, data1$all_curves[i, ], col = col_correct, lwd = 1.5)
  } else {
    lines(data1$t_grid, data1$all_curves[i, ], col = rgb(0, 0.5, 0, 0.6), lwd = 1.5)
  }
}

# Add mean curves
lines(data1$t_grid, data1$mean1, col = "black", lwd = 3, lty = 1)
lines(data1$t_grid, data1$mean2, col = "black", lwd = 3, lty = 2)

legend("topright",
       legend = c("Correct Pop 1", "Correct Pop 2", "True Means"),
       col = c(col_correct, rgb(0, 0.5, 0, 0.6), "darkgreen"),
       lwd = c(1.5, 1.5, 3),
       bty = "n", cex = 0.8)

grid(col = "lightgray", lty = "dotted")


# ------------------------------------------------------------------
# PLOT 4: Incorrectly Predicted (if any)
# ------------------------------------------------------------------
if (n_incorrect > 0) {
  plot(data1$t_grid, data1$all_curves[1, ], type = "n",
       ylim = range(data1$all_curves),
       xlab = "Time", ylab = "Value",
       main = sprintf("Incorrect Predictions (%d/%d)",
                      n_incorrect, length(correct_predictions)))

  # Plot only incorrectly predicted curves
  incorrect_idx <- which(!correct_predictions)

  for (i in incorrect_idx) {
    if (data1$true_labels[i] == 1) {
      # Should be in cluster 1 but predicted as 2
      lines(data1$t_grid, data1$all_curves[i, ], col = "orange", lwd = 2)
      points(data1$t_grid, data1$all_curves[i, ], pch = 1, col = "blue", cex = 0.5)
    } else {
      # Should be in cluster 2 but predicted as 1
      lines(data1$t_grid, data1$all_curves[i, ], col = "darkorange", lwd = 2)
      points(data1$t_grid, data1$all_curves[i, ], pch = 1, col = "red", cex = 0.5)
    }
  }

  # Add reference means
  lines(data1$t_grid, data1$mean1, col = "black", lwd = 2, lty = 1)
  lines(data1$t_grid, data1$mean2, col = "black", lwd = 2, lty = 2)

  legend("topright",
         legend = c("Misclass. from Pop 1", "Misclass. from Pop 2", "True Means"),
         col = c("orange", "darkorange", "gray"),
         lwd = c(2, 2, 2),
         bty = "n", cex = 0.8)

  grid(col = "lightgray", lty = "dotted")
} else {
  plot.new()
  text(0.5, 0.5, "Perfect Clustering!\nNo Misclassifications",
       cex = 1.5, col = "darkgreen", font = 2)
}


# ------------------------------------------------------------------
# PLOT 5: Projected Data Space
# ------------------------------------------------------------------
# Show the 1D projection space where clustering was performed
proj_data <- result$projected_data[, 1]

# Create histogram/density plot
hist(proj_data[predicted_labels_aligned == 1],
     breaks = 15, col = rgb(0, 0, 1, 0.5),
     xlim = range(proj_data),
     xlab = "Projected Value",
     main = "Projected Data Distribution",
     border = "white")
hist(proj_data[predicted_labels_aligned == 2],
     breaks = 15, col = rgb(1, 0, 0, 0.5),
     add = TRUE, border = "white")

# Add cluster centers
abline(v = result$centers[1], col = "blue", lwd = 3, lty = 2)
abline(v = result$centers[2], col = "red", lwd = 3, lty = 2)

# Mark misclassified points
if (n_incorrect > 0) {
  incorrect_idx <- which(!correct_predictions)
  for (i in incorrect_idx) {
    abline(v = proj_data[i], col = "orange", lwd = 1, lty = 3)
  }
}

legend("topright",
       legend = c("Cluster 1", "Cluster 2", "Centers",
                  if(n_incorrect > 0) "Misclassified"),
       col = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5), "gray",
               if(n_incorrect > 0) "orange"),
       lwd = c(10, 10, 3, if(n_incorrect > 0) 1),
       lty = c(1, 1, 2, if(n_incorrect > 0) 3),
       bty = "n", cex = 0.8)

grid(col = "lightgray", lty = "dotted")


# ------------------------------------------------------------------
# PLOT 6: Mean Profiles Comparison
# ------------------------------------------------------------------
plot(data1$t_grid, data1$mean1, type = "l", col = "blue", lwd = 3,
     ylim = range(c(data1$mean1, data1$mean2, pred_mean1, pred_mean2)),
     xlab = "Time", ylab = "Mean Value",
     main = "True vs Predicted Mean Profiles")

lines(data1$t_grid, data1$mean2, col = "red", lwd = 3)
lines(data1$t_grid, pred_mean1, col = "blue", lwd = 2, lty = 2)
lines(data1$t_grid, pred_mean2, col = "red", lwd = 2, lty = 2)

legend("topright",
       legend = c("True Mean 1", "True Mean 2",
                  "Predicted Mean 1", "Predicted Mean 2"),
       col = c("blue", "red", "blue", "red"),
       lwd = c(3, 3, 2, 2),
       lty = c(1, 1, 2, 2),
       bty = "n", cex = 0.8)

grid(col = "lightgray", lty = "dotted")

# Add RMSE of mean profiles
rmse1 <- sqrt(mean((data1$mean1 - pred_mean1)^2))
rmse2 <- sqrt(mean((data1$mean2 - pred_mean2)^2))
text(mean(data1$t_grid), min(data1$mean1, data1$mean2),
     sprintf("RMSE: Pop1=%.3f, Pop2=%.3f", rmse1, rmse2),
     cex = 0.8, col = "darkgray")


# Add overall title
title(main = sprintf("Functional K-Means Clustering Results (Accuracy: %.1f%%)", accuracy),
      outer = TRUE, cex.main = 1.3, font.main = 2)


# =============================================================================
# DETAILED BREAKDOWN BY INDIVIDUAL
# =============================================================================

cat("\n\n=== Detailed Classification Breakdown ===\n\n")

# Create detailed table
breakdown <- data.frame(
  ID = 1:length(data1$true_labels),
  True_Label = data1$true_labels,
  Predicted_Label = predicted_labels_aligned,
  Correct = correct_predictions,
  Projected_Value = proj_data
)

# Show first few and any misclassified
cat("First 10 observations:\n")
print(head(breakdown, 10))

if (n_incorrect > 0) {
  cat("\n\nMisclassified observations:\n")
  print(breakdown[!correct_predictions, ])

  # Analyze misclassifications
  cat("\n\nMisclassification Analysis:\n")
  misclass_from_1 <- sum(data1$true_labels[!correct_predictions] == 1)
  misclass_from_2 <- sum(data1$true_labels[!correct_predictions] == 2)

  cat(sprintf("  - From Population 1: %d (%.1f%% of Pop 1)\n",
              misclass_from_1,
              100 * misclass_from_1 / sum(data1$true_labels == 1)))
  cat(sprintf("  - From Population 2: %d (%.1f%% of Pop 2)\n",
              misclass_from_2,
              100 * misclass_from_2 / sum(data1$true_labels == 2)))
}


# =============================================================================
# ADDITIONAL METRICS
# =============================================================================

cat("\n\n=== Additional Clustering Metrics ===\n\n")

# Sensitivity and Specificity (treating Pop 1 as "positive")
TP <- sum(predicted_labels_aligned == 1 & data1$true_labels == 1)
FP <- sum(predicted_labels_aligned == 1 & data1$true_labels == 2)
TN <- sum(predicted_labels_aligned == 2 & data1$true_labels == 2)
FN <- sum(predicted_labels_aligned == 2 & data1$true_labels == 1)

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
ppv <- TP / (TP + FP)
npv <- TN / (TN + FN)

cat(sprintf("Sensitivity (Pop 1): %.2f%%\n", sensitivity * 100))
cat(sprintf("Specificity (Pop 2): %.2f%%\n", specificity * 100))
cat(sprintf("Positive Predictive Value: %.2f%%\n", ppv * 100))
cat(sprintf("Negative Predictive Value: %.2f%%\n", npv * 100))

# Clustering quality metrics
cat(sprintf("\nTightness: %.4f\n", result$min_tightness))
cat(sprintf("Total within-cluster SS: %.2f\n", result$total_within_ss))
cat(sprintf("Between-cluster SS: %.2f\n", result$kmeans_result$betweenss))
cat(sprintf("Total SS: %.2f\n", result$kmeans_result$totss))
cat(sprintf("Between SS / Total SS: %.2f%%\n",
            100 * result$kmeans_result$betweenss / result$kmeans_result$totss))

cat("\n=== Analysis Complete ===\n\n")
