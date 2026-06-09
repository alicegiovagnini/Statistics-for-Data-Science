# STATISTICS FOR DATA SCIENCE — PROJECT A.Y. 2025/26
# Paper: "A General Procedure to Combine Estimators"
#        Lavancier & Rochet (2015), Computational Statistics & Data Analysis
#
# FILE 1/4 — SETUP AND CORE THEORY (Sections 2-3 of the paper)
#
# This file is sourced by each of the three application scripts:
#   01_section_4_1.R   — averaging of mean and median
#   03_section_4_3.R   — Boolean model estimation
#   04_section_4_4.R   — quantile estimation under misspecification
#
# It loads the packages, sets up parallelism and defines the averaging
# procedure (the shared "theory"). It produces no output on its own.

# PART 0: Setup

# Fix the global RNG seed so the whole script is reproducible run-to-run.
set.seed(42)

# Packages the script depends on: plotting (ggplot2/gridExtra), table
# printing (knitr), multicore Monte Carlo (parallel), the QP solver used as a
# fallback for the convex averaging weights (quadprog) and showtext, which lets
# the plots use the Cambria font of the presentation slides.
required_packages <- c("ggplot2", "gridExtra", "knitr", "parallel", "quadprog",
                       "showtext")
# Install any package that is not already available, quietly, before loading.
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# Attach the packages so their functions are on the search path.
library(ggplot2)     # grammar-of-graphics plots
library(gridExtra)   # arrange several ggplot objects on one page
library(knitr)       # kable() for clean console tables
library(parallel)    # mclapply() for parallel Monte Carlo replicates
library(quadprog)    # solve.QP() fallback for constrained weights
library(showtext)    # render custom (Cambria) fonts in the plots

# Use every CPU core but one (leaves the machine responsive); on Windows
# fork-based mclapply is unavailable, so set n_cores <- 1L there.
n_cores <- max(1L, detectCores() - 1L)
# Report the chosen degree of parallelism so the user can sanity-check it.
cat(sprintf("Available cores for parallelization: %d\n", n_cores))

# Paper-style cell formatter: combine an MSE value and its standard deviation
# into a single string "mse (sd)", as in the tables of Lavancier & Rochet (2015).
# `digits` controls the rounding of both numbers.
mse_sd <- function(mse, sd, digits = 2) {
  sprintf("%.*f (%.*f)", digits, mse, digits, sd)
}


# PRESENTATION STYLE (shared by all plots)
#
# The plots are exported into slides that use the Cambria font and a navy/gold
# colour scheme; the helpers below keep every figure consistent with them.

# Register Cambria from the Microsoft PowerPoint bundle (the font is not in the
# usual macOS font folders). If it cannot be loaded, fall back to a generic
# serif so the scripts still run on any machine.
PRES_FONT <- "serif"
local({
  ppt <- "/Applications/Microsoft PowerPoint.app/Contents/Resources/DFonts"
  reg <- file.path(ppt, "Cambria.ttc")
  if (file.exists(reg)) {
    tryCatch({
      font_add("Cambria",
               regular    = reg,
               bold       = file.path(ppt, "Cambriab.ttf"),
               italic     = file.path(ppt, "Cambriai.ttf"),
               bolditalic = file.path(ppt, "Cambriaz.ttf"))
      showtext_auto()              # render text with showtext on every device
      showtext_opts(dpi = 96)      # match the RStudio on-screen graphics device
      PRES_FONT <<- "Cambria"
    }, error = function(e)
      message("Cambria could not be loaded, using serif instead: ",
              conditionMessage(e)))
  }
})

# Slide colour palette, sampled from the four header boxes of the deck:
# navy, red, green and gold. A fifth muted purple covers the only plot that has
# five series (Section 4.4); the other plots never reach it.
PRES_PALETTE <- c("#1E245A", "#D2483C", "#427836", "#EAA842", "#8C5AA8")

# Assign palette colours to a vector of estimator labels, in order. Returns a
# named vector suitable for scale_*_manual().
estimator_colors <- function(levels) {
  cols        <- PRES_PALETTE[seq_along(levels)]
  names(cols) <- levels
  cols
}

# Shared ggplot theme: Cambria font, white background, light horizontal grid
# only, and large text (axis titles, tick labels and legend) so the figures stay
# very legible once placed on a slide with plenty of room.
theme_pres <- function(base_size = 17) {
  theme_minimal(base_size = base_size, base_family = PRES_FONT) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88"),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.title         = element_text(size = rel(1.15), face = "bold"),
      axis.title         = element_text(size = rel(1.30), color = "grey20"),
      axis.text          = element_text(size = rel(1.10), color = "grey25"),
      legend.title       = element_text(size = rel(1.30), face = "bold"),
      legend.text        = element_text(size = rel(1.20)),
      legend.key.size    = unit(1.3, "lines")
    )
}


# PART 1 — AVERAGING PROCEDURE (Section 2 of the paper)
#
# Given k estimators T1,...,Tk of a parameter theta, the averaging estimator is:
#   theta_hat = lambda_hat' * T
# with optimal weights (eq. 9):
#   lambda* = Sigma^{-1} 1 / (1' Sigma^{-1} 1)
# subject to sum(lambda) = 1 (Lambda_max, univariate case).
#
# For d parameters simultaneously (multivariate, eq. 11):
#   theta_hat = (J' Sigma^{-1} J)^{-1} J' Sigma^{-1} T
# where J assigns estimators to parameters.

# Optimal averaging weights (univariate, Lambda_max, eq. 9)
#' @param Sigma_hat  k x k estimated MSE matrix
#' @return  numeric vector of length k (weights summing to 1)
averaging_weights <- function(Sigma_hat) {
  # Column of ones, one entry per estimator (the "1" vector in eq. 9).
  ones        <- rep(1, nrow(Sigma_hat))
  # Solve Sigma * y = 1 instead of forming Sigma^{-1} explicitly:
  # more stable and gives the numerator Sigma^{-1} 1 directly.
  Sigma_inv_1 <- solve(Sigma_hat, ones)
  # Normalise by 1' Sigma^{-1} 1 so the weights sum to one (the constraint).
  Sigma_inv_1 / as.numeric(t(ones) %*% Sigma_inv_1)
}

# Averaging estimator (univariate)
#' @param T_vec      k-vector of point estimates
#' @param Sigma_hat  k x k estimated MSE matrix
#' @return  scalar averaging estimate
averaging_estimator <- function(T_vec, Sigma_hat) {
  # Combined estimate = weighted sum of the individual estimators.
  as.numeric(averaging_weights(Sigma_hat) %*% T_vec)
}

# Asymptotic variance of the averaging estimator (Proposition 3.3)
#' @param Sigma_hat  k x k estimated MSE matrix
#' @return  scalar estimated variance
averaging_variance <- function(Sigma_hat) {
  # Recompute the optimal weights for this Sigma.
  w <- averaging_weights(Sigma_hat)
  # Var(theta_hat) = w' Sigma w; used later to build confidence intervals.
  as.numeric(t(w) %*% Sigma_hat %*% w)
}

# Multivariate averaging estimator (eq. 11)
#' @param T_vec      k-vector of point estimates
#' @param Sigma_hat  k x k estimated MSE matrix
#' @param J          k x d assignment matrix (J[i,j]=1 if T_i estimates theta_j)
#' @return  d-vector of averaging estimates
averaging_multivariate <- function(T_vec, Sigma_hat, J) {
  # Invert the MSE matrix once; reused in both A and b below.
  Sigma_inv <- solve(Sigma_hat)
  # A = J' Sigma^{-1} J : the d x d "information" matrix across parameters.
  A         <- t(J) %*% Sigma_inv %*% J
  # b = J' Sigma^{-1} T : projects the estimators onto the parameter space.
  b_vec     <- t(J) %*% Sigma_inv %*% T_vec
  # Solve A theta = b to get the d combined estimates (eq. 11).
  solve(A, b_vec)
}

# Convex averaging: min lambda' Sigma lambda  s.t. sum=1, lambda>=0
#' Implements the iterative algorithm from Section 2.3, with quadprog fallback.
#' @param T_vec      k-vector of point estimates
#' @param Sigma_hat  k x k estimated MSE matrix
#' @return  list with: theta (estimate), weights, support (active indices)
averaging_convex <- function(T_vec, Sigma_hat) {
  # Number of estimators to combine.
  k      <- length(T_vec)
  # "Active set": indices currently allowed a non-zero weight; start with all.
  active <- seq_len(k)

  # Active-set loop: at most k iterations, dropping one estimator each time
  # a weight would go negative, until all remaining weights are non-negative.
  for (iter in seq_len(k)) {
    # Restrict Sigma and the ones-vector to the currently active estimators.
    S_sub  <- Sigma_hat[active, active, drop = FALSE]
    ones   <- rep(1, length(active))

    # Unconstrained optimal direction Sigma^{-1} 1 on the active block;
    # tryCatch guards against a singular sub-matrix.
    Sinv_1 <- tryCatch(solve(S_sub, ones), error = function(e) NULL)
    # Bail out of the loop if the sub-problem could not be solved.
    if (is.null(Sinv_1)) break

    # Normalising constant 1' Sigma^{-1} 1.
    denom <- as.numeric(t(ones) %*% Sinv_1)
    # If it is not positive/finite the solution is invalid; stop and fall back.
    if (!is.finite(denom) || denom <= 0) break

    # Candidate weights on the active set (the unconstrained eq.-9 solution).
    w_sub <- Sinv_1 / denom

    # If every weight is (numerically) non-negative, the constraint lambda>=0
    # is satisfied and we have the convex solution.
    if (all(w_sub >= -1e-10)) {
      # Expand back to the full length-k weight vector, zeros off the support.
      w_full         <- rep(0, k)
      w_full[active] <- pmax(w_sub, 0)
      # Renormalise to defend against tiny negative round-off being clipped.
      w_full         <- w_full / sum(w_full)
      # Return estimate, weights and which estimators are actually used.
      return(list(theta   = as.numeric(w_full %*% T_vec),
                  weights = w_full,
                  support = active))
    }

    # Otherwise drop the estimator with the most negative weight and retry.
    remove_idx <- which.min(w_sub)
    active     <- active[-remove_idx]
    # If we removed everything, give up on the active-set method.
    if (length(active) == 0) break
  }

  # Fallback: solve the constrained QP directly when the active-set loop fails.
  tryCatch({
    # Symmetrise Sigma and add a tiny ridge so the QP matrix is positive definite.
    D_mat <- (Sigma_hat + t(Sigma_hat)) / 2 + diag(1e-8, k)
    # Constraints: first column enforces sum(lambda)=1 (equality), the diagonal
    # block enforces lambda >= 0 (inequalities).
    A_mat <- cbind(rep(1, k), diag(k))
    # Right-hand sides matching the constraints above.
    b_vec <- c(1, rep(0, k))
    # Minimise lambda' D lambda with linear term 0; meq=1 marks the first
    # constraint as an equality (the sum-to-one condition).
    sol   <- quadprog::solve.QP(D_mat, rep(0, k), A_mat, b_vec, meq = 1)
    # Clip negatives from numerical noise and renormalise to sum one.
    w_full <- pmax(sol$solution, 0)
    w_full <- w_full / sum(w_full)
    # Same return shape as the active-set branch.
    list(theta   = as.numeric(w_full %*% T_vec),
         weights = w_full,
         support = which(w_full > 1e-10))
  }, error = function(e) {
    # Last resort: if even the QP fails, keep the single lowest-MSE estimator.
    best          <- which.min(diag(Sigma_hat))
    w_full        <- rep(0, k)
    w_full[best]  <- 1
    list(theta = T_vec[best], weights = w_full, support = best)
  })
}
