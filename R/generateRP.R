#' Visualize Basis Support Regions (B-splines)
#'
#' @param rp Object of class random_projections from generateRP()
#' @param threshold Numeric. Values below this threshold are considered zero (default: 1e-6)
#' @param col Character. Color for support regions (default: "lightblue")
#' @param alpha Numeric. Transparency for support regions (default: 0.3)
#' @param ... Additional arguments passed to plot()
#'
#' @details
#' This function visualizes the support regions (non-zero intervals) of each
#' basis function, which is particularly useful for B-splines. The plot shows
#' which time intervals each basis function covers, helping users understand
#' how projections can achieve local feature selection.
#'
#' @examples
#' \dontrun{
#'   # Visualize support regions for B-splines
#'   Q_B <- generateRP(B = 10, d0 = 2, r = 8, basis = "bspline")
#'   plot_basis_support(Q_B)
#'
#'   # With custom knots
#'   knots <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
#'   Q_B <- generateRP(B = 10, d0 = 2, r = 8, basis = "bspline",
#'                    bspline_knots = knots)
#'   plot_basis_support(Q_B)
#' }
#'
#'
#' @export
plot_basis_support <- function(rp, threshold = 1e-6,
                               col = "lightblue", alpha = 0.3, ...) {

  if (!inherits(rp, "random_projections")) {
    stop("Input must be an object of class 'random_projections'")
  }

  basis_functions <- get_basis_functions(rp)
  grid <- attr(rp, "grid")
  r <- attr(rp, "r")
  basis_type <- attr(rp, "basis_type")

  # Calculate support regions
  support_regions <- vector("list", r)
  for (i in 1:r) {
    non_zero_idx <- which(abs(basis_functions[i, ]) > threshold)
    if (length(non_zero_idx) > 0) {
      support_regions[[i]] <- c(grid[min(non_zero_idx)],
                                grid[max(non_zero_idx)])
    } else {
      support_regions[[i]] <- c(NA, NA)
    }
  }

  # Create the plot
  plot(1, type = "n", xlim = range(grid), ylim = c(0.5, r + 0.5),
       xlab = "Time/Position", ylab = "Basis Function Index",
       main = paste("Support Regions -", tools::toTitleCase(basis_type), "Basis"),
       yaxt = "n", ...)

  axis(2, at = 1:r, labels = 1:r, las = 1)

  # Add support regions as rectangles
  col_with_alpha <- adjustcolor(col, alpha.f = alpha)
  for (i in 1:r) {
    if (!any(is.na(support_regions[[i]]))) {
      rect(support_regions[[i]][1], i - 0.4,
           support_regions[[i]][2], i + 0.4,
           col = col_with_alpha, border = col)

      # Add center point
      center <- mean(support_regions[[i]])
      points(center, i, pch = 19, col = col, cex = 0.8)
    }
  }

  # Add grid lines
  abline(h = 1:r, col = "gray90", lty = 2)
  grid(nx = 10, ny = 0, col = "gray80", lty = "dotted")

  # Add legend
  legend("topright",
         legend = c("Support region", "Center"),
         fill = c(col_with_alpha, NA),
         border = c(col, NA),
         pch = c(NA, 19),
         col = c(NA, col),
         bty = "n")

  # Add summary text
  support_widths <- sapply(support_regions, function(x) {
    if (!any(is.na(x))) x[2] - x[1] else NA
  })

  mtext(sprintf("Mean support width: %.3f", mean(support_widths, na.rm = TRUE)),
        side = 3, line = 0.5, cex = 0.8, col = "gray40")
}#' Generate Random Projections for Functional Data Clustering
#'
#' @param B Integer. Number of random projections to generate
#' @param d0 Integer. Dimension of multivariate functional data
#' @param r Integer. Number of basis functions to use for each projection
#' @param basis Character or matrix. Basis type to use:
#'   \itemize{
#'     \item "bspline" - B-spline basis (default, recommended for omics data)
#'     \item "fourier" - Fourier basis
#'     \item matrix - Custom r x grid_length matrix of basis functions
#'   }
#' @param grid Vector. Grid points where basis functions are evaluated
#' @param bspline_degree Integer. Degree of B-splines if basis = "bspline" (default 3)
#' @param bspline_knots Vector. Optional knot positions for B-splines. If NULL,
#'        uses equally spaced internal knots. Knots should be in range of grid.
#'
#' @return List with B elements, each containing:
#'   \item{Psi}{List of d0 projection functions (each as coefficient vector)}
#'   \item{coefficients}{d0 x r matrix of normalized coefficients}
#'   \item{projection_id}{Integer identifying the projection}
#'
#' @details
#' B-splines are the default and recommended basis for longitudinal omics data because:
#' \itemize{
#'   \item Local support: Each basis function is zero outside its support interval
#'   \item Feature selection: Projections can be zero-valued over certain intervals
#'   \item Flexibility: Can focus on specific time windows of biological interest
#' }
#'
#' Custom knot placement allows you to concentrate basis functions in regions
#' where you expect important biological changes.
#'
#' @examples
#' \dontrun{
#' # Generate projections with B-spline basis (default)
#' projections <- generateRP(B = 100, d0 = 3, r = 10)
#'
#' # B-spline with custom knot placement (focus on early time points)
#' grid <- seq(0, 1, length.out = 100)
#' knots <- c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 1)  # More knots early
#' projections <- generateRP(B = 100, d0 = 3, r = 10,
#'                          basis = "bspline", bspline_knots = knots, grid = grid)
#'
#' # Fourier basis for periodic patterns
#' projections <- generateRP(B = 100, d0 = 3, r = 10, basis = "fourier")
#'
#' # Custom basis
#' custom_basis <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
#' projections <- generateRP(B = 100, d0 = 3, r = 10,
#'                          basis = custom_basis, grid = grid)
#' }
generateRP <- function(B, d0, r,
                       basis = "bspline",
                       grid = seq(0, 1, length.out = 100),
                       bspline_degree = 3,
                       bspline_knots = NULL) {

  # Validate inputs
  if (B <= 0 || d0 <= 0 || r <= 0) {
    stop("B, d0, and r must be positive integers")
  }

  # Handle basis specification
  if (is.character(basis)) {
    basis <- tolower(basis)

    if (basis == "fourier") {
      basis_functions <- create_fourier_basis(r, grid)
      basis_type <- "fourier"
    } else if (basis == "bspline") {
      basis_functions <- create_bspline_basis(r, grid,
                                              degree = bspline_degree,
                                              knots = bspline_knots)
      basis_type <- "bspline"
    } else {
      stop("basis must be 'fourier', 'bspline', or a matrix")
    }

  } else if (is.matrix(basis)) {
    # User provided custom basis
    basis_functions <- basis
    basis_type <- "custom"

    # Validate custom basis dimensions
    if (nrow(basis_functions) != r) {
      stop(paste("Custom basis matrix must have r =", r, "rows"))
    }
    if (ncol(basis_functions) != length(grid)) {
      stop(paste("Custom basis matrix must have", length(grid),
                 "columns to match grid length"))
    }

  } else {
    stop("basis must be a character string ('fourier' or 'bspline') or a matrix")
  }

  # VECTORIZED APPROACH: Generate all B*d0 random coefficient vectors at once
  # Draw all random coefficients in one call: B*d0 rows, r columns
  all_coefficients <- matrix(stats::rnorm(B * d0 * r), nrow = B * d0, ncol = r)

  # Normalize each row to unit length using vectorized operations
  row_norms <- sqrt(rowSums(all_coefficients^2))
  all_coefficients <- all_coefficients / row_norms

  # Reshape into 3D array: B x d0 x r for easier indexing
  coef_array <- array(all_coefficients, dim = c(d0, B, r))
  coef_array <- aperm(coef_array, c(2, 1, 3))  # Permute to B x d0 x r

  # Pre-allocate the entire list structure
  Q_B <- vector("list", B)

  # Build projection objects efficiently
  for (b in 1:B) {
    # Extract coefficient matrix for this projection (d0 x r)
    # Handle d0 = 1 case where subsetting drops dimension
    if (d0 == 1) {
      coef_matrix <- matrix(coef_array[b, , ], nrow = 1, ncol = r)
    } else {
      coef_matrix <- coef_array[b, , ]
    }

    # Create list of d0 projection functions
    Psi_b <- lapply(1:d0, function(j) {
      list(
        coefficients = coef_matrix[j, ],
        basis = basis_functions,
        grid = grid
      )
    })

    # Store this projection
    Q_B[[b]] <- list(
      Psi = Psi_b,
      coefficients = coef_matrix,
      projection_id = b
    )
  }

  class(Q_B) <- c("random_projections", "list")
  attr(Q_B, "d0") <- d0
  attr(Q_B, "r") <- r
  attr(Q_B, "B") <- B
  attr(Q_B, "basis_type") <- basis_type
  attr(Q_B, "basis_functions") <- basis_functions
  attr(Q_B, "grid") <- grid
  attr(Q_B, "bspline_degree") <- if (basis_type == "bspline") bspline_degree else NULL
  attr(Q_B, "bspline_knots") <- if (basis_type == "bspline") bspline_knots else NULL

  return(Q_B)
}


#' Create Fourier Basis Functions
#'
#' @param r Integer. Number of basis functions
#' @param grid Vector. Grid points for evaluation
#'
#' @return Matrix of r x length(grid) with basis functions in rows
#'
create_fourier_basis <- function(r, grid) {

  n_grid <- length(grid)
  basis_matrix <- matrix(0, nrow = r, ncol = n_grid)

  # Constant function
  if (r >= 1) {
    basis_matrix[1, ] <- 1 / sqrt(grid[n_grid] - grid[1])
  }

  # Sine and cosine pairs
  k <- 1
  idx <- 2
  while (idx <= r) {
    period <- 2 * pi * k / (grid[n_grid] - grid[1])

    # Sine function
    if (idx <= r) {
      basis_matrix[idx, ] <- sqrt(2 / (grid[n_grid] - grid[1])) *
        sin(period * grid)
      idx <- idx + 1
    }

    # Cosine function
    if (idx <= r) {
      basis_matrix[idx, ] <- sqrt(2 / (grid[n_grid] - grid[1])) *
        cos(period * grid)
      idx <- idx + 1
    }

    k <- k + 1
  }

  return(basis_matrix)
}


#' Create B-spline Basis Functions
#'
#' @param r Integer. Number of basis functions
#' @param grid Vector. Grid points for evaluation
#' @param degree Integer. Degree of B-splines (default 3 for cubic)
#' @param knots Vector. Optional internal knot positions. If NULL, uses equally
#'        spaced knots. Should be within range of grid.
#'
#' @return Matrix of r x length(grid) with basis functions in rows
#'
#' @details
#' B-splines have local support, meaning each basis function is non-zero only
#' over a subset of the domain. This property is ideal for:
#' \itemize{
#'   \item Feature selection in time-series data
#'   \item Identifying important time windows
#'   \item Creating projections that focus on specific intervals
#' }
#'
#' Custom knot placement allows concentration of basis functions in regions
#' of interest (e.g., early response periods in omics experiments).
#'
#' @examples
#' \dontrun{
#'   grid <- seq(0, 1, length.out = 100)
#'
#'   # Equally spaced knots (default)
#'   basis1 <- create_bspline_basis(r = 10, grid = grid)
#'
#'   # Custom knots - more density early
#'   knots <- c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1)
#'   basis2 <- create_bspline_basis(r = 10, grid = grid, knots = knots)
#' }
#'
create_bspline_basis <- function(r, grid, degree = 3, knots = NULL) {

  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is required for B-spline basis")
  }

  # Validate inputs
  if (r < degree + 1) {
    stop(paste("r must be at least degree + 1 =", degree + 1))
  }

  grid_range <- range(grid)

  # Handle knot specification
  if (!is.null(knots)) {
    # User-provided knots
    if (any(knots < grid_range[1] | knots > grid_range[2])) {
      warning("Some knots are outside grid range. They will be clipped.")
      knots <- pmax(grid_range[1], pmin(grid_range[2], knots))
    }

    # Ensure knots include boundaries
    if (knots[1] != grid_range[1]) {
      knots <- c(grid_range[1], knots)
    }
    if (knots[length(knots)] != grid_range[2]) {
      knots <- c(knots, grid_range[2])
    }

    # Remove duplicates and sort
    knots <- sort(unique(knots))

    # Calculate number of internal knots needed
    n_internal_knots <- r - degree - 1
    n_provided_internal <- length(knots) - 2  # Exclude boundaries

    if (n_provided_internal > n_internal_knots) {
      warning(paste("Too many knots provided. Using first",
                    n_internal_knots, "internal knots."))
      internal_knots <- knots[2:(n_internal_knots + 1)]
    } else if (n_provided_internal < n_internal_knots) {
      # Add equally spaced knots to reach desired number
      n_to_add <- n_internal_knots - n_provided_internal
      additional_knots <- seq(grid_range[1], grid_range[2],
                              length.out = n_to_add + 2)[2:(n_to_add + 1)]
      internal_knots <- sort(c(knots[2:(length(knots)-1)], additional_knots))
    } else {
      internal_knots <- knots[2:(length(knots)-1)]
    }

    # Create B-spline basis with specified knots
    basis_eval <- splines::bs(grid, knots = internal_knots,
                              degree = degree,
                              Boundary.knots = grid_range,
                              intercept = TRUE)
  } else {
    # Default: equally spaced internal knots
    basis_eval <- splines::bs(grid, df = r, degree = degree,
                              intercept = TRUE)
  }

  # Transpose so basis functions are in rows
  basis_matrix <- t(basis_eval)

  # Normalize each basis function to have unit L2 norm
  grid_diff <- diff(grid)
  for (i in 1:nrow(basis_matrix)) {
    # Approximate L2 norm using trapezoidal rule
    norm_sq <- sum(basis_matrix[i, -1]^2 * grid_diff +
                     basis_matrix[i, -length(grid)]^2 * grid_diff) / 2
    if (norm_sq > 1e-10) {  # Avoid division by zero
      basis_matrix[i, ] <- basis_matrix[i, ] / sqrt(norm_sq)
    }
  }

  return(basis_matrix)
}


#' Evaluate Projection Function on Functional Data
#'
#' @param X Matrix. n x grid_length matrix of functional observations
#' @param psi List. Projection function from generate_random_projections
#'
#' @return Vector of n projected values
#'
project_functional_data <- function(X, psi) {

  # Extract basis and coefficients
  chi <- psi$basis  # r x grid_length
  a <- psi$coefficients  # length r
  grid <- psi$grid

  # Compute psi(x) = sum_m a_m * chi_m(x)
  psi_values <- as.vector(a %*% chi)  # length grid_length

  # Compute X_i(psi) = integral X_i(t) * psi(t) dt for each i
  # Using trapezoidal rule
  n <- nrow(X)
  grid_diff <- diff(grid)

  projections <- numeric(n)
  for (i in 1:n) {
    integrand <- X[i, ] * psi_values
    # Trapezoidal rule
    projections[i] <- sum((integrand[-1] + integrand[-length(integrand)]) *
                            grid_diff) / 2
  }

  return(projections)
}


#' Project Multivariate Functional Data onto Psi
#'
#' @param X_list List of length d0. Each element is n x grid_length matrix
#' @param Psi_star List. One element from Q_B (contains d0 projection functions)
#'
#' @return Matrix of n x d0 projected values
#'
project_multivariate_data <- function(X_list, Psi_star) {

  d0 <- length(X_list)
  n <- nrow(X_list[[1]])

  # Initialize projection matrix
  X_projected <- matrix(0, nrow = n, ncol = d0)

  # Project each dimension
  for (j in 1:d0) {
    X_projected[, j] <- project_functional_data(X_list[[j]], Psi_star$Psi[[j]])
  }

  return(X_projected)
}


#' Print method for random_projections
#'
#' @param x Object of class random_projections
#' @param ... Additional arguments (unused)
#'
print.random_projections <- function(x, ...) {
  B <- attr(x, "B")
  d0 <- attr(x, "d0")
  r <- attr(x, "r")
  basis_type <- attr(x, "basis_type")
  grid <- attr(x, "grid")

  cat("Random Projections for Functional Data Clustering\n")
  cat("=================================================\n")
  cat(sprintf("Number of projections (B): %d\n", B))
  cat(sprintf("Data dimension (d0): %d\n", d0))
  cat(sprintf("Basis truncation (r): %d\n", r))
  cat(sprintf("Basis type: %s\n", basis_type))

  if (basis_type == "bspline") {
    degree <- attr(x, "bspline_degree")
    knots <- attr(x, "bspline_knots")
    cat(sprintf("  B-spline degree: %d\n", degree))
    if (!is.null(knots)) {
      cat(sprintf("  Custom knots: %d specified\n", length(knots)))
    } else {
      cat("  Knots: equally spaced (default)\n")
    }
  }

  cat(sprintf("Grid length: %d points from %.3f to %.3f\n",
              length(grid), min(grid), max(grid)))
  cat("\nEach projection consists of", d0, "orthonormal functions\n")
  cat("represented as linear combinations of", r, "basis functions.\n")

  if (basis_type == "bspline") {
    cat("\nB-splines provide local support for feature selection:\n")
    cat("  - Each basis function is non-zero only over a subset of the domain\n")
    cat("  - Projections can focus on specific time windows\n")
  }

  cat("\nUse get_basis_functions() to extract the basis functions.\n")
  cat("Use plot() to visualize basis functions and/or projections.\n")
}


#' Plot Method for Random Projections
#'
#' @param x Object of class random_projections from generateRP()
#' @param type Character. Type of plot: "both" (default), "basis", or "projections"
#' @param which_projection Integer. Which projection to visualize (default: 1)
#' @param which_dimension Integer vector. Which dimensions to show for projections
#'        (default: all if d0 <= 4, otherwise 1:4)
#' @param which_basis Integer vector. Which basis functions to plot (default: all if r <= 10)
#' @param col_basis Character vector or function. Colors for basis functions
#' @param col_proj Character vector or function. Colors for projection functions
#' @param lwd Numeric. Line width (default: 2)
#' @param ... Additional arguments passed to plotting functions
#'
#' @examples
#' \dontrun{
#'   Q_B <- generateRP(B = 50, d0 = 3, r = 8, basis = "fourier")
#'
#'   # Plot both basis and projection functions
#'   plot(Q_B)
#'
#'   # Plot only basis functions
#'   plot(Q_B, type = "basis")
#'
#'   # Plot only projection functions
#'   plot(Q_B, type = "projections")
#'
#'   # Customize which projection and dimensions to show
#'   plot(Q_B, which_projection = 5, which_dimension = 1:2)
#' }
#'
#' @export
plot.random_projections <- function(x,
                                    type = c("both", "basis", "projections"),
                                    which_projection = 1,
                                    which_dimension = NULL,
                                    which_basis = NULL,
                                    col_basis = NULL,
                                    col_proj = NULL,
                                    lwd = 2,
                                    ...) {

  type <- match.arg(type)

  # Validate which_projection
  B <- attr(x, "B")
  if (which_projection < 1 || which_projection > B) {
    stop(paste("which_projection must be between 1 and", B))
  }

  if (type == "basis") {
    # Plot only basis functions
    plot_basis_functions(x, which_basis = which_basis,
                         col = col_basis, lwd = lwd, ...)

  } else if (type == "projections") {
    # Plot only projection functions
    plot_projection_functions(x, which_projection = which_projection,
                              which_dimension = which_dimension,
                              col = col_proj, lwd = lwd, ...)

  } else {  # "both"
    # Save original par settings
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    # Set up 1x2 layout
    graphics::par(mfrow = c(1, 2))

    # Plot basis functions on left
    plot_basis_functions(x, which_basis = which_basis,
                         col = col_basis, lwd = lwd,
                         legend = FALSE, ...)

    # Plot projection functions on right
    plot_projection_functions(x, which_projection = which_projection,
                              which_dimension = which_dimension,
                              col = col_proj, lwd = lwd, ...)
  }

  invisible(x)
}


#' Extract Basis Functions from Random Projections Object
#'
#' @param rp Object of class random_projections from generateRP()
#'
#' @return Matrix of r x grid_length with basis functions in rows
#'
#' @examples
#' \dontrun{
#'   Q_B <- generateRP(B = 10, d0 = 2, r = 5)
#'   basis <- get_basis_functions(Q_B)
#'   dim(basis)  # 5 x 100
#' }
#'
#' @export
get_basis_functions <- function(rp) {
  if (!inherits(rp, "random_projections")) {
    stop("Input must be an object of class 'random_projections'")
  }

  basis_functions <- attr(rp, "basis_functions")

  if (is.null(basis_functions)) {
    stop("No basis functions found in the random_projections object")
  }

  return(basis_functions)
}


#' Plot Basis Functions
#'
#' @param rp Object of class random_projections from generateRP()
#' @param which_basis Integer vector. Which basis functions to plot (default: all)
#' @param col Character vector or function. Colors for the basis functions
#' @param lty Integer vector. Line types for the basis functions
#' @param lwd Numeric. Line width (default: 2)
#' @param main Character. Main title for the plot
#' @param xlab Character. X-axis label
#' @param ylab Character. Y-axis label
#' @param legend Logical. Whether to add a legend (default: TRUE if r <= 10)
#' @param ... Additional arguments passed to plot()
#'
#' @examples
#' \dontrun{
#'   # Plot Fourier basis
#'   Q_fourier <- generateRP(B = 10, d0 = 2, r = 8, basis = "fourier")
#'   plot_basis_functions(Q_fourier)
#'
#'   # Plot only first 5 B-spline basis functions
#'   Q_bspline <- generateRP(B = 10, d0 = 2, r = 12, basis = "bspline")
#'   plot_basis_functions(Q_bspline, which_basis = 1:5)
#' }
#'
#' @export
plot_basis_functions <- function(rp, which_basis = NULL,
                                 col = NULL, lty = 1, lwd = 2,
                                 main = NULL, xlab = "t", ylab = "Basis function value",
                                 legend = NULL, ...) {

  if (!inherits(rp, "random_projections")) {
    stop("Input must be an object of class 'random_projections'")
  }

  # Extract basis functions and grid
  basis_functions <- get_basis_functions(rp)
  grid <- attr(rp, "grid")
  r <- attr(rp, "r")
  basis_type <- attr(rp, "basis_type")

  # Determine which basis functions to plot
  if (is.null(which_basis)) {
    which_basis <- 1:r
  } else {
    if (any(which_basis < 1 | which_basis > r)) {
      stop(paste("which_basis must be between 1 and", r))
    }
  }

  n_plot <- length(which_basis)

  # Set up colors
  if (is.null(col)) {
    if (n_plot <= 8) {
      col <- 1:n_plot
    } else {
      col <- grDevices::rainbow(n_plot)
    }
  }

  # Recycle colors and line types if needed
  col <- rep_len(col, n_plot)
  lty <- rep_len(lty, n_plot)

  # Set default title
  if (is.null(main)) {
    main <- paste(tools::toTitleCase(basis_type), "Basis Functions")
  }

  # Determine legend display
  if (is.null(legend)) {
    legend <- (n_plot <= 10)
  }

  # Calculate y-axis limits
  y_range <- range(basis_functions[which_basis, ])

  # Create the plot
  plot(grid, basis_functions[which_basis[1], ], type = "l",
       col = col[1], lty = lty[1], lwd = lwd,
       ylim = y_range, main = main, xlab = xlab, ylab = ylab, ...)

  # Add remaining basis functions
  if (n_plot > 1) {
    for (i in 2:n_plot) {
      graphics::lines(grid, basis_functions[which_basis[i], ],
            col = col[i], lty = lty[i], lwd = lwd)
    }
  }

  # Add legend if requested
  if (legend) {
    legend("topright",
           legend = paste("Basis", which_basis),
           col = col, lty = lty, lwd = lwd,
           bty = "n", cex = 0.8)
  }

  # Add grid
  grid(col = "lightgray", lty = "dotted")
}


#' Plot Projection Functions
#'
#' @param rp Object of class random_projections from generateRP()
#' @param which_projection Integer. Which projection to visualize (default: 1)
#' @param which_dimension Integer vector. Which dimensions to show
#'        (default: all if d0 <= 4, otherwise 1:4)
#' @param col Character vector or function. Colors for projection functions
#' @param lty Integer vector. Line types
#' @param lwd Numeric. Line width (default: 2)
#' @param main Character. Main title
#' @param xlab Character. X-axis label
#' @param ylab Character. Y-axis label
#' @param legend Logical. Whether to add a legend (default: TRUE)
#' @param show_components Logical. If TRUE, also show individual basis components
#'        (default: FALSE)
#' @param ... Additional arguments passed to plot()
#'
#' @examples
#' \dontrun{
#'   Q_B <- generateRP(B = 50, d0 = 3, r = 8, basis = "fourier")
#'   plot_projection_functions(Q_B, which_projection = 1)
#'
#'   # Show only first 2 dimensions
#'   plot_projection_functions(Q_B, which_dimension = 1:2)
#'
#'   # Show with basis components
#'   plot_projection_functions(Q_B, show_components = TRUE)
#' }
#'
#' @export
plot_projection_functions <- function(rp, which_projection = 1,
                                      which_dimension = NULL,
                                      col = NULL, lty = 1, lwd = 2,
                                      main = NULL, xlab = "t",
                                      ylab = "Projection function value",
                                      legend = TRUE,
                                      show_components = FALSE,
                                      ...) {

  if (!inherits(rp, "random_projections")) {
    stop("Input must be an object of class 'random_projections'")
  }

  # Extract attributes
  B <- attr(rp, "B")
  d0 <- attr(rp, "d0")
  r <- attr(rp, "r")
  grid <- attr(rp, "grid")
  basis_functions <- get_basis_functions(rp)
  basis_type <- attr(rp, "basis_type")

  # Validate which_projection
  if (which_projection < 1 || which_projection > B) {
    stop(paste("which_projection must be between 1 and", B))
  }

  # Extract the selected projection
  projection <- rp[[which_projection]]
  coef_matrix <- projection$coefficients  # d0 x r matrix

  # Ensure coef_matrix is always a matrix (handles d0 = 1 case)
  if (!is.matrix(coef_matrix)) {
    coef_matrix <- matrix(coef_matrix, nrow = 1)
  }

  # Determine which dimensions to plot
  if (is.null(which_dimension)) {
    if (d0 <= 4) {
      which_dimension <- 1:d0
    } else {
      which_dimension <- 1:4
      message(paste("d0 =", d0, "> 4, plotting only first 4 dimensions.",
                    "Use which_dimension to specify others."))
    }
  } else {
    if (any(which_dimension < 1 | which_dimension > d0)) {
      stop(paste("which_dimension must be between 1 and", d0))
    }
  }

  n_plot <- length(which_dimension)

  # Set up colors
  if (is.null(col)) {
    if (n_plot <= 8) {
      col <- 1:n_plot
    } else {
      col <- grDevices::rainbow(n_plot)
    }
  }
  col <- rep_len(col, n_plot)
  lty <- rep_len(lty, n_plot)

  # Compute projection functions: psi_j(t) = sum_m a_{j,m} * chi_m(t)
  # Result is d0 x grid_length matrix
  projection_functions <- coef_matrix %*% basis_functions

  # Ensure projection_functions is a matrix (handles d0 = 1 case)
  if (!is.matrix(projection_functions)) {
    projection_functions <- matrix(projection_functions, nrow = 1)
  }

  # Set default title
  if (is.null(main)) {
    main <- sprintf("Projection Functions (Projection #%d)", which_projection)
  }

  # Calculate y-axis limits
  y_range <- range(projection_functions[which_dimension, ])

  # Create the plot
  plot(grid, projection_functions[which_dimension[1], ],
       type = "l", col = col[1], lty = lty[1], lwd = lwd,
       ylim = y_range, main = main, xlab = xlab, ylab = ylab, ...)

  # Add remaining projection functions
  if (n_plot > 1) {
    for (i in 2:n_plot) {
      graphics::lines(grid, projection_functions[which_dimension[i], ],
            col = col[i], lty = lty[i], lwd = lwd)
    }
  }

  # Optionally show individual basis components
  if (show_components && n_plot == 1) {
    # Only show components if plotting a single dimension
    j <- which_dimension[1]
    for (m in 1:r) {
      if (abs(coef_matrix[j, m]) > 1e-6) {  # Only show non-negligible components
        component <- coef_matrix[j, m] * basis_functions[m, ]
        graphics::lines(grid, component, col = "gray70", lty = 2, lwd = 1)
      }
    }
  }

  # Add legend
  if (legend) {
    legend_labels <- paste("Dimension", which_dimension)
    legend("topright",
           legend = legend_labels,
           col = col, lty = lty, lwd = lwd,
           bty = "n", cex = 0.8)
  }

  # Add grid
  grid(col = "lightgray", lty = "dotted")

  # Add zero line for reference
  abline(h = 0, col = "gray50", lty = 3)
}


#' Get Projection Functions as Matrix
#'
#' @param rp Object of class random_projections from generateRP()
#' @param which_projection Integer. Which projection to extract (default: 1)
#'
#' @return Matrix of d0 x grid_length with projection functions in rows
#'
#' @examples
#' \dontrun{
#'   Q_B <- generateRP(B = 10, d0 = 3, r = 5)
#'   proj_funcs <- get_projection_functions(Q_B, which_projection = 1)
#'   dim(proj_funcs)  # 3 x 100
#' }
#'
#' @export
get_projection_functions <- function(rp, which_projection = 1) {

  if (!inherits(rp, "random_projections")) {
    stop("Input must be an object of class 'random_projections'")
  }

  B <- attr(rp, "B")
  if (which_projection < 1 || which_projection > B) {
    stop(paste("which_projection must be between 1 and", B))
  }

  # Extract basis and coefficients
  basis_functions <- get_basis_functions(rp)
  coef_matrix <- rp[[which_projection]]$coefficients

  # Ensure coef_matrix is a matrix (handles d0 = 1 case)
  if (!is.matrix(coef_matrix)) {
    coef_matrix <- matrix(coef_matrix, nrow = 1)
  }

  # Compute projection functions
  projection_functions <- coef_matrix %*% basis_functions

  # Ensure result is a matrix (handles d0 = 1 case)
  if (!is.matrix(projection_functions)) {
    projection_functions <- matrix(projection_functions, nrow = 1)
  }

  return(projection_functions)
}


#' Summary Method for Random Projections
#'
#' @param object Object of class random_projections
#' @param ... Additional arguments (unused)
#'
#' @export
summary.random_projections <- function(object, ...) {
  B <- attr(object, "B")
  d0 <- attr(object, "d0")
  r <- attr(object, "r")
  basis_type <- attr(object, "basis_type")
  grid <- attr(object, "grid")
  basis_functions <- attr(object, "basis_functions")

  cat("Summary of Random Projections\n")
  cat("==============================\n\n")

  cat("Projection Configuration:\n")
  cat(sprintf("  Number of projections: %d\n", B))
  cat(sprintf("  Data dimension: %d\n", d0))
  cat(sprintf("  Basis functions used: %d\n", r))
  cat(sprintf("  Basis type: %s\n", basis_type))

  if (basis_type == "bspline") {
    degree <- attr(object, "bspline_degree")
    knots <- attr(object, "bspline_knots")
    cat(sprintf("  B-spline degree: %d\n", degree))
    if (!is.null(knots)) {
      cat(sprintf("  Custom knots: %d positions\n", length(knots)))
      cat("  Knot locations:", paste(round(knots, 3), collapse = ", "), "\n")
    } else {
      cat("  Knots: equally spaced internal knots\n")
    }
  }

  cat(sprintf("  Grid: [%.3f, %.3f] with %d points\n",
              min(grid), max(grid), length(grid)))

  cat("\nBasis Function Statistics:\n")
  basis_norms <- apply(basis_functions, 1, function(x) {
    grid_diff <- diff(grid)
    sqrt(sum((x[-1]^2 + x[-length(x)]^2) * grid_diff / 2))
  })

  cat(sprintf("  L2 norm range: [%.4f, %.4f]\n",
              min(basis_norms), max(basis_norms)))
  cat(sprintf("  Mean L2 norm: %.4f\n", mean(basis_norms)))

  # Calculate support (non-zero regions) for B-splines
  if (basis_type == "bspline") {
    supports <- apply(basis_functions, 1, function(x) {
      non_zero <- which(abs(x) > 1e-10)
      if (length(non_zero) > 0) {
        c(grid[min(non_zero)], grid[max(non_zero)])
      } else {
        c(NA, NA)
      }
    })
    support_widths <- supports[2, ] - supports[1, ]
    cat(sprintf("  Basis support width: [%.3f, %.3f] (mean: %.3f)\n",
                min(support_widths, na.rm = TRUE),
                max(support_widths, na.rm = TRUE),
                mean(support_widths, na.rm = TRUE)))
  }

  cat("\nCoefficient Statistics (across all projections):\n")
  all_coefs <- sapply(object, function(x) as.vector(x$coefficients))
  cat(sprintf("  Mean coefficient: %.4f\n", mean(all_coefs)))
  cat(sprintf("  SD of coefficients: %.4f\n", stats::sd(all_coefs)))
  cat(sprintf("  Coefficient range: [%.4f, %.4f]\n",
              min(all_coefs), max(all_coefs)))

  invisible(object)
}
