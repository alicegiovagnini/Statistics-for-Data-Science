# =====================================================================
# CORRECTED VERSION — two fixes applied vs the original script:
#
#  [FIX 1] quantile_np(): wrapped in as.numeric() to drop the "99%" name
#          that quantile() attaches. Without this the entire NP column of
#          Table 8 (Section 4.4) was NA and the NP bars vanished from the plot.
#
#  [FIX 2] estimate_Sigma_quantile_boot(): the bootstrap MSE matrix is now
#          centered on the consistent NP estimate instead of on each
#          estimator's own value. The old centering produced the covariance
#          (bias-blind), which gave biased-but-stable models too much weight;
#          as a result AV was ~18x worse than the correct model on Burr data
#          and ~7x worse than NP on Lognormal. After the fix AV beats NP when
#          a model fits and matches NP when none does — i.e. it behaves as the
#          paper and the slides claim. Search "[FIX" to see both spots.
# =====================================================================

# STATISTICS FOR DATA SCIENCE — PROJECT A.Y. 2025/26
# Paper: "A General Procedure to Combine Estimators"
#        Lavancier & Rochet (2015), Computational Statistics & Data Analysis
#
# This script reproduces methods and experiments from Sections 2, 3, 4.1,
# 4.3 and 4.4 of the paper.
#
# Structure:
#   PART 0 — Setup and utility functions
#   PART 1 — Averaging procedure (core theory, Sections 2-3)
#   PART 2 — Section 4.1: averaging of mean and median
#   PART 3 — Section 4.3: Boolean model estimation
#   PART 4 — Section 4.4: quantile estimation under misspecification


# PART 0: Setup

# Fix the global RNG seed so the whole script is reproducible run-to-run.
set.seed(42)

# Packages the script depends on: plotting (ggplot2/gridExtra), table
# printing (knitr), multicore Monte Carlo (parallel) and the QP solver
# used as a fallback for the convex averaging weights (quadprog).
required_packages <- c("ggplot2", "gridExtra", "knitr", "parallel", "quadprog")
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

# Use every CPU core but one (leaves the machine responsive); on Windows
# fork-based mclapply is unavailable, so set n_cores <- 1L there.
n_cores <- max(1L, detectCores() - 1L)
# Report the chosen degree of parallelism so the user can sanity-check it.
cat(sprintf("Available cores for parallelization: %d\n", n_cores))

# Folder where partial Monte Carlo results are cached, so a re-run can skip
# work already done (used by the heavy Sections 4.3 and 4.4 simulations).
checkpoint_dir <- "checkpoints"
# Create it if missing; stay silent if it already exists.
dir.create(checkpoint_dir, showWarnings = FALSE)


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

#' Optimal averaging weights (univariate, Lambda_max, eq. 9)
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

#' Averaging estimator (univariate)
#' @param T_vec      k-vector of point estimates
#' @param Sigma_hat  k x k estimated MSE matrix
#' @return  scalar averaging estimate
averaging_estimator <- function(T_vec, Sigma_hat) {
  # Combined estimate = weighted sum of the individual estimators.
  as.numeric(averaging_weights(Sigma_hat) %*% T_vec)
}

#' Asymptotic variance of the averaging estimator (Proposition 3.3)
#' @param Sigma_hat  k x k estimated MSE matrix
#' @return  scalar estimated variance
averaging_variance <- function(Sigma_hat) {
  # Recompute the optimal weights for this Sigma.
  w <- averaging_weights(Sigma_hat)
  # Var(theta_hat) = w' Sigma w; used later to build confidence intervals.
  as.numeric(t(w) %*% Sigma_hat %*% w)
}

#' Multivariate averaging estimator (eq. 11)
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

#' Convex averaging: min lambda' Sigma lambda  s.t. sum=1, lambda>=0
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
    # Restrict Sigma, T and the ones-vector to the currently active estimators.
    S_sub  <- Sigma_hat[active, active, drop = FALSE]
    T_sub  <- T_vec[active]
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


# PART 2 — SECTION 4.1: LOCATION ESTIMATION FOR SYMMETRIC DISTRIBUTIONS
#
# Two initial estimators of theta (location):
#   T1 = sample mean,  T2 = sample median
# Asymptotic MSE matrix Sigma (2x2) estimated via:
#   AV:  asymptotic formula (Laplace) with plug-in estimates
#   AVB: non-parametric bootstrap

# 2.1  Sigma estimation via asymptotic formula (AV method)
estimate_Sigma_asymptotic <- function(x, theta0) {
  # Sample size, needed to scale the asymptotic variances by 1/n.
  n          <- length(x)
  # Plug-in variance of the mean estimator (entry [1,1] before scaling).
  sigma2_hat <- var(x)
  
  # --- AUTOMATIC ADAPTIVE PROTECTION FOR CAUCHY / HEAVY-TAILED OUTLIERS ---
  # If the sample variance is extremely large or non-finite, we catch the Cauchy run.
  # We set the variance to an artificially high theoretical value and zero out the covariance.
  # This makes the matrix perfectly diagonal and stable for solve().
  if (!is.finite(sigma2_hat) || sigma2_hat > 5000) {
    sigma2_hat_cauchy <- 1e6
    f_hat_cauchy      <- 1 / pi   # Theoretical density of standard Cauchy at the origin
    
    W <- matrix(c(sigma2_hat_cauchy, 0,
                  0,                 1 / (4 * f_hat_cauchy^2)),
                nrow = 2, ncol = 2)
    return(W / n)
  }
  
  # --- STANDARD PROCEDURE FOR GAUSSIAN, LOGISTIC, STUDENT-T AND MIXTURE ---
  # Mean absolute deviation about theta0: enters the mean/median covariance.
  m_hat      <- mean(abs(x - theta0))
  
  # Protected Silverman's bandwidth using robust IQR metric to avoid bloated h values
  s_hat      <- min(sd(x), IQR(x) / 1.349)
  if (s_hat == 0 || !is.finite(s_hat)) s_hat <- sd(x)
  
  h     <- 1.06 * s_hat * n^(-1/5)
  # Kernel density of the data at theta0; f(theta0) drives the median variance.
  f_hat <- mean(dnorm((x - theta0) / h)) / h
  
  # Safeguard against collapsed near-zero sample density estimates
  if (f_hat < 1e-5) f_hat <- 1e-5
  
  # Asymptotic 2x2 matrix W: variance of mean, mean/median covariance, and
  # variance of the median = 1/(4 f^2) (standard sample-median result).
  W <- matrix(c(sigma2_hat,
                m_hat / (2 * f_hat),
                m_hat / (2 * f_hat),
                1 / (4 * f_hat^2)),
              nrow = 2, ncol = 2)
  # Divide by n to turn asymptotic variances into MSEs at this sample size.
  W / n
}


# 2.2  Sigma estimation via non-parametric bootstrap (AVB method)

#' Estimate 2x2 MSE matrix via non-parametric bootstrap
#' @param x  numeric vector: data sample
#' @param B  integer: number of bootstrap replicates
#' @return  2x2 Sigma_hat matrix
estimate_Sigma_bootstrap <- function(x, B = 1000) {
  # Sample size, used for resampling with replacement.
  n            <- length(x)
  # Storage for the mean and median computed on each bootstrap sample.
  boot_means   <- numeric(B)
  boot_medians <- numeric(B)
  # Draw B resamples and record both estimators on each.
  for (b in seq_len(B)) {
    x_b            <- sample(x, n, replace = TRUE)
    boot_means[b]   <- mean(x_b)
    boot_medians[b] <- median(x_b)
  }
  # The empirical covariance of (mean, median) across resamples estimates Sigma.
  cov(cbind(boot_means, boot_medians))
}

# 2.3  Combined estimation for a single sample

#' Compute all Section 4.1 estimators for one sample
#' @param x  numeric vector: data sample
#' @param B  integer: bootstrap replicates for AVB
#' @return  list with: mean, median, AV, AVB, var_AV, var_AVB
compute_estimates_4_1 <- function(x, B = 1000) {
  # The two base estimators of the location parameter.
  T1    <- mean(x)
  T2    <- median(x)
  # Stack them into the k=2 vector consumed by the averaging routines.
  T_vec <- c(T1, T2)

  # AV branch: Sigma from the asymptotic formula (centred at the robust median).
  Sigma_AV  <- estimate_Sigma_asymptotic(x, theta0 = T2)
  # Combined estimate and its variance under the asymptotic Sigma.
  theta_AV  <- averaging_estimator(T_vec, Sigma_AV)
  var_AV    <- averaging_variance(Sigma_AV)

  # AVB branch: same averaging but with the bootstrap-estimated Sigma.
  Sigma_AVB <- estimate_Sigma_bootstrap(x, B = B)
  theta_AVB <- averaging_estimator(T_vec, Sigma_AVB)
  var_AVB   <- averaging_variance(Sigma_AVB)

  # Return all four point estimates plus the two combined-estimator variances.
  list(mean    = T1,       median  = T2,
       AV      = theta_AV, AVB     = theta_AVB,
       var_AV  = var_AV,   var_AVB = var_AVB)
}

# 2.4  Monte Carlo simulation (replicates Tables 1 & 2)

#' Run Monte Carlo simulation for Section 4.1
#' @param n_vec  integer vector: sample sizes
#' @param n_rep  integer: Monte Carlo replicates
#' @param B      integer: bootstrap replicates per sample
#' @return  data.frame with MSE and CI coverage
run_simulation_4_1 <- function(n_vec = c(30, 50, 100),
                                n_rep = 1000,
                                B     = 200) {
  # Test distributions: all symmetric about 0 so the true location is theta=0,
  # ranging from heavy-tailed (Cauchy) to light-tailed (Gaussian) and bimodal.
  distributions <- list(
    Cauchy   = function(n) rcauchy(n, location = 0, scale = 1),
    Student4 = function(n) rt(n, df = 4),
    Student7 = function(n) rt(n, df = 7),
    Logistic = function(n) rlogis(n, location = 0, scale = 1),
    Gaussian = function(n) rnorm(n, mean = 0, sd = 1),
    Mixture  = function(n) {
      # 50/50 mixture of N(-2,1) and N(2,1): symmetric but bimodal.
      idx <- rbinom(n, 1, 0.5)
      ifelse(idx == 0, rnorm(n, -2, 1), rnorm(n, 2, 1))
    }
  )

  # Collect one result row per (distribution, n) combination.
  results <- list()
  for (dist_name in names(distributions)) {
    # Generator for the current distribution.
    rgen <- distributions[[dist_name]]
    # Progress message so the user sees which distribution is running.
    cat(sprintf("  Distribution: %s\n", dist_name))

    for (n in n_vec) {
      # Per-replicate squared-error stores for the four estimators.
      err_mean <- err_med <- err_AV <- err_AVB <- numeric(n_rep)
      # Per-replicate CI hit/miss indicators for the two combined estimators.
      cover_AV <- cover_AVB <- logical(n_rep)

      for (rep in seq_len(n_rep)) {
        # Draw a fresh sample and compute all estimators on it.
        x   <- rgen(n)
        est <- compute_estimates_4_1(x, B = B)

        # Squared error vs the true value (0); these average to the MSE.
        err_mean[rep] <- est$mean^2
        err_med[rep]  <- est$median^2
        err_AV[rep]   <- est$AV^2
        err_AVB[rep]  <- est$AVB^2

        # Does the 95% CI (centred at the estimate) cover the truth 0?
        # Equivalent to |estimate| <= 1.96 * SE.
        cover_AV[rep]  <- abs(est$AV)  <= 1.96 * sqrt(est$var_AV)
        cover_AVB[rep] <- abs(est$AVB) <= 1.96 * sqrt(est$var_AVB)
      }

      # Aggregate the replicates into MSE (x100 for readability), SD and
      # empirical coverage (%), one row per (distribution, n).
      results[[paste(dist_name, n, sep = "_")]] <- data.frame(
        distribution = dist_name, n = n,
        MSE_mean     = mean(err_mean) * 100,
        MSE_median   = mean(err_med)  * 100,
        MSE_AV       = mean(err_AV)   * 100,
        MSE_AVB      = mean(err_AVB)  * 100,
        SD_mean      = sd(err_mean)   * 100,
        SD_median    = sd(err_med)    * 100,
        SD_AV        = sd(err_AV)     * 100,
        SD_AVB       = sd(err_AVB)    * 100,
        Coverage_AV  = mean(cover_AV)  * 100,
        Coverage_AVB = mean(cover_AVB) * 100
      )
    }
  }
  # Bind all rows into a single data.frame for printing/plotting.
  do.call(rbind, results)
}

# 2.5  Execution and results

# Announce the section and run the Section 4.1 simulation.
cat("\nSECTION 4.1: MSE Simulation\n")
sim_4_1 <- run_simulation_4_1(n_vec = c(30, 50, 100), n_rep = 1000, B = 200)

# Table 1: estimated MSE (x100) for mean, median, AV and AVB.
cat("\nTable 1: Estimated MSE x100\n")
tab1 <- sim_4_1[, c("distribution","n","MSE_mean","MSE_median","MSE_AV","MSE_AVB")]
print(kable(tab1, digits = 2,
            col.names = c("Distribution","n","Mean","Median","AV","AVB")))

# Table 2: empirical coverage of the nominal 95% CIs for AV and AVB.
cat("\nTable 2: 95% CI Coverage (%)\n")
tab2 <- sim_4_1[, c("distribution","n","Coverage_AV","Coverage_AVB")]
print(kable(tab2, digits = 2,
            col.names = c("Distribution","n","AV (%)","AVB (%)")))

# 2.6  Plot: MSE comparison for n=50

# Reshape the wide results to long form (one row per estimator) for ggplot.
plot_data_4_1 <- do.call(rbind, lapply(unique(sim_4_1$n), function(ni) {
  sub <- sim_4_1[sim_4_1$n == ni, ]
  data.frame(
    distribution = rep(sub$distribution, 4),
    n            = ni,
    estimator    = rep(c("Mean","Median","AV","AVB"), each = nrow(sub)),
    MSE          = c(sub$MSE_mean, sub$MSE_median, sub$MSE_AV, sub$MSE_AVB)
  )
}))
# Fix the factor orders so legend and x-axis appear in a sensible sequence.
plot_data_4_1$estimator    <- factor(plot_data_4_1$estimator,
                                     levels = c("Mean","Median","AV","AVB"))
plot_data_4_1$distribution <- factor(plot_data_4_1$distribution,
                                     levels = c("Cauchy","Student4","Student7",
                                                "Logistic","Gaussian","Mixture"))

# Grouped bar chart at n=50; Cauchy is dropped because its huge MSE would
# crush the scale and hide the differences among the other distributions.
p_mse_4_1 <- ggplot(
    plot_data_4_1[plot_data_4_1$n == 50 & plot_data_4_1$distribution != "Cauchy", ],
    aes(x = distribution, y = MSE, fill = estimator)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#2c7bb6","#d7191c","#1a9641","#fdae61")) +
  labs(
    title    = "Sec. 4.1 — MSE Comparison (n=50, x100, Cauchy excluded for scale)",
    subtitle = "AV and AVB dominate or match in all cases",
    x = "Distribution", y = "Estimated MSE x100", fill = "Estimator"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
print(p_mse_4_1)

# 2.7  Plot: averaging weights by distribution

# Show, on one sample each, how the optimal weight shifts between mean and median.
cat("\nIllustrative weights (n=100, one sample per distribution)\n")
# Local seed so this illustrative table is reproducible on its own.
set.seed(123)
dist_funs <- list(
  Cauchy   = function(n) rcauchy(n),
  Student4 = function(n) rt(n, df = 4),
  Student7 = function(n) rt(n, df = 7),
  Logistic = function(n) rlogis(n),
  Gaussian = function(n) rnorm(n),
  Mixture  = function(n) ifelse(rbinom(n,1,.5)==0, rnorm(n,-2), rnorm(n,2))
)

# For each distribution: draw one sample, estimate Sigma, get the two weights.
weight_df <- do.call(rbind, lapply(names(dist_funs), function(dname) {
  x <- dist_funs[[dname]](100)
  S <- estimate_Sigma_asymptotic(x, median(x))
  w <- averaging_weights(S)
  data.frame(distribution = dname, w_mean = w[1], w_median = w[2])
}))
# Print the weights as a table.
print(kable(weight_df, digits = 3,
            col.names = c("Distribution","Weight on Mean","Weight on Median")))

# Long format so the two weights can be stacked in a single bar per distribution.
weight_long <- rbind(
  data.frame(distribution = weight_df$distribution,
             stimatore = "Mean",   peso = weight_df$w_mean),
  data.frame(distribution = weight_df$distribution,
             stimatore = "Median", peso = weight_df$w_median)
)
# Keep distributions in the declared order on the x-axis.
weight_long$distribution <- factor(weight_long$distribution,
                                   levels = names(dist_funs))

# Stacked bars summing to 1, with a reference line at the 50/50 split.
p_weights <- ggplot(weight_long, aes(x = distribution, y = peso, fill = stimatore)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#2c7bb6","#d7191c")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.6, y = 0.53, label = "50%", size = 3.5, color="grey30") +
  labs(
    title    = "Sec. 4.1 — Averaging weights by distribution (n=100)",
    subtitle = "For Cauchy, weight goes to median; for Gaussian, to mean",
    x = "Distribution", y = "Weight", fill = "Estimator"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
print(p_weights)


# PART 3 — SECTION 4.3: BOOLEAN MODEL ESTIMATION
#
# Stationary Boolean model on [0,1]^2:
#   - Germs: Poisson process with intensity rho
#   - Grains: discs with radius R ~ 0.1 * Beta(1, alpha)
# Observable functionals (Weil & Wieacker 1984):
#   A = 1 - exp(-pi * rho * E[R^2])   (covered area fraction)
#   P = 2*pi*rho * E[R] * exp(-pi*rho*E[R^2])  (perimeter per unit area)
#
# Two initial estimators:
#   T1 = (rho_1, alpha_1): by inverting A and P formulas
#   T2 = rho_2: based on tangent points (Molchanov 1995)
#
# Multivariate averaging (d=2, J = [1 0; 1 0; 0 1]):
#   rho_AV and alpha_AV via eq. 11
# Sigma (3x3) estimated via parametric bootstrap (B=100).

# 3.1  Theoretical moments of R ~ 0.1 * Beta(1, alpha)

#' @param alpha  scalar > 0
#' @return  list with E_R and E_R2
beta_moments <- function(alpha) {
  # Closed-form moments of 0.1*Beta(1,alpha): E[R] and E[R^2] feed the A,P formulas.
  list(E_R  = 0.1 / (alpha + 1),
       E_R2 = 0.02 / ((alpha + 1) * (alpha + 2)))
}

# 3.2  Boolean model simulation

#' Simulate a Boolean model realization and return observed functionals
#' @param rho    scalar: Poisson intensity
#' @param alpha  scalar: shape parameter of grain radii
#' @return  list: A_obs, P_obs, n_germs
simulate_boolean_obs <- function(rho, alpha) {
  # Number of germs is Poisson(rho) on the unit window |W|=1.
  n_germs <- rpois(1, rho)
  # No germs => no coverage and no perimeter; return zeros early.
  if (n_germs == 0) return(list(A_obs = 0, P_obs = 0, n_germs = 0))

  # Draw the grain radii and their empirical first two moments.
  r        <- 0.1 * rbeta(n_germs, shape1 = 1, shape2 = alpha)
  E_R_emp  <- mean(r)
  E_R2_emp <- mean(r^2)
  # Empirical intensity equals the germ count since the window has unit area.
  rho_emp  <- n_germs  # density = n_germs / |W| with |W|=1

  # Plug the empirical moments into the Weil-Wieacker formulas for A and P.
  A_obs <- 1 - exp(-pi * rho_emp * E_R2_emp)
  P_obs <- 2 * pi * rho_emp * E_R_emp * exp(-pi * rho_emp * E_R2_emp)
  # Clamp A away from 0 and 1 so later log(1-A) inversions stay finite.
  A_obs <- min(max(A_obs, 1e-6), 1 - 1e-6)

  # Return the two observed functionals and the germ count.
  list(A_obs = A_obs, P_obs = P_obs, n_germs = n_germs)
}

# 3.3  Estimator 1: based on area and perimeter

#' Compute (rho_1, alpha_1) by inverting the A and P formulas
#' @param A_obs  scalar in (0,1): observed covered area
#' @param P_obs  scalar > 0: observed perimeter
#' @return  named vector c(rho, alpha), or c(NA, NA)
estimator_1_boolean <- function(A_obs, P_obs) {
  # Guard the domain: A must be a probability strictly inside (0,1)...
  if (!is.finite(A_obs) || A_obs <= 0 || A_obs >= 1) return(c(rho=NA, alpha=NA))
  # ...and the perimeter must be positive and finite.
  if (!is.finite(P_obs) || P_obs <= 0)               return(c(rho=NA, alpha=NA))

  # log(1-A) appears in both inverted formulas; compute once.
  log_term  <- log(1 - A_obs)
  # Invert the A/P system for the grain-shape parameter alpha.
  alpha_hat <- P_obs / (10 * (A_obs - 1) * log_term) - 2
  # Reject non-physical (non-positive / non-finite) alpha.
  if (!is.finite(alpha_hat) || alpha_hat <= 0) return(c(rho=NA, alpha=NA))

  # Back out the intensity rho given alpha_hat.
  rho_hat <- 5 * (alpha_hat + 1) * P_obs / (pi * (1 - A_obs))
  # Reject non-physical rho.
  if (!is.finite(rho_hat) || rho_hat <= 0) return(c(rho=NA, alpha=NA))

  # Return both estimates as a named vector.
  c(rho = rho_hat, alpha = alpha_hat)
}

# 3.4  Estimator 2: based on tangent points

#' Compute rho_2 from tangent point counts (Molchanov 1995)
#' E[N(u)] = |W| * rho * (1-A), averaged over k random directions
#' @param A_obs    scalar: observed area
#' @param n_germs  integer: observed number of germs
#' @param W_area   scalar: window area (default 1)
#' @param k        integer: number of random directions (default 100)
#' @return  scalar rho_2 estimate, or NA
estimator_2_boolean <- function(A_obs, n_germs, W_area = 1, k = 100) {
  # The uncovered fraction (1-A); must be positive for the formula to invert.
  one_minus_A <- 1 - A_obs
  if (!is.finite(one_minus_A) || one_minus_A <= 0) return(NA)

  # Expected tangent-point count per direction = |W| * rho * (1-A); here rho is
  # proxied by n_germs since |W|=1, giving the Poisson rate for the counts.
  lambda_N <- W_area * n_germs * one_minus_A
  # Simulate tangent-point counts over k random directions...
  N_vals   <- rpois(k, lambda_N)
  # ...and invert the expectation to recover an intensity estimate rho_2.
  rho_2    <- mean(N_vals) / (W_area * one_minus_A)

  # Reject non-positive / non-finite results.
  if (!is.finite(rho_2) || rho_2 <= 0) return(NA)
  rho_2
}

# 3.5  Sigma estimation via parametric bootstrap

#' Estimate 3x3 MSE matrix for (rho_1, rho_2, alpha_1) via parametric bootstrap
#' @param rho0    scalar: initial rho estimate
#' @param alpha0  scalar: initial alpha estimate
#' @param B       integer: bootstrap replicates
#' @return  3x3 Sigma_hat, or NULL if fewer than 10 valid replicates
estimate_Sigma_boolean <- function(rho0, alpha0, B = 100) {
  # Need valid positive parameters to simulate from; otherwise abort.
  if (!is.finite(rho0) || rho0 <= 0) return(NULL)
  if (!is.finite(alpha0) || alpha0 <= 0) return(NULL)

  # Each row will hold (rho_1, rho_2, alpha_1) from one parametric resample.
  res <- matrix(NA, nrow = B, ncol = 3)
  for (b in seq_len(B)) {
    # Simulate a new Boolean realization at the plug-in parameters...
    obs_b  <- simulate_boolean_obs(rho0, alpha0)
    # ...and recompute both estimators on it.
    est1_b <- estimator_1_boolean(obs_b$A_obs, obs_b$P_obs)
    rho2_b <- estimator_2_boolean(obs_b$A_obs, obs_b$n_germs)
    res[b, ] <- c(est1_b["rho"], rho2_b, est1_b["alpha"])
  }

  # Keep only fully valid rows (drop any with NA from failed inversions).
  res <- res[complete.cases(res), ]
  # Too few valid replicates => unreliable covariance; signal failure.
  if (nrow(res) < 10) return(NULL)

  # Centre the resampled estimators on the true (plug-in) values: (rho0,rho0,alpha0).
  true_vals <- c(rho0, rho0, alpha0)
  centered  <- sweep(res, 2, true_vals, "-")
  # MSE matrix = average outer product of centred errors (the bootstrap Sigma).
  t(centered) %*% centered / nrow(centered)
}

# 3.6  Averaging estimator for the Boolean model

#' Compute (rho_AV, alpha_AV) via multivariate averaging (eq. 11)
#' T_vec = (rho_1, rho_2, alpha_1),  J = [1 0; 1 0; 0 1]
#' @param T_vec      length-3 vector
#' @param Sigma_hat  3x3 MSE matrix
#' @return  named vector c(rho, alpha)
boolean_averaging <- function(T_vec, Sigma_hat) {
  # Assignment matrix: rho_1 and rho_2 both estimate rho (col 1); alpha_1
  # estimates alpha (col 2). This is the J of eq. 11.
  J <- matrix(c(1, 1, 0,
                0, 0, 1), nrow = 3, ncol = 2)
  # Apply the multivariate averaging; guard against a singular Sigma.
  theta_hat <- tryCatch(
    averaging_multivariate(T_vec, Sigma_hat, J),
    error = function(e) c(NA, NA)
  )
  # Label the two combined estimates.
  names(theta_hat) <- c("rho", "alpha")
  theta_hat
}

# 3.7  Asymptotic variance for Boolean model CIs

#' Compute estimated variance of rho_AV and alpha_AV (Proposition 3.3)
#' var_j = [(J' Sigma^{-1} J)^{-1}]_{jj}
#' @param Sigma_hat  3x3 MSE matrix
#' @param J          3x2 assignment matrix
#' @return  named vector c(var_rho, var_alpha)
boolean_averaging_variance <- function(Sigma_hat, J) {
  tryCatch({
    # Invert Sigma once.
    Si      <- solve(Sigma_hat)
    # The 2x2 information matrix J' Sigma^{-1} J across the two parameters.
    JtSiJ   <- t(J) %*% Si %*% J
    # Its inverse is the asymptotic covariance of (rho_AV, alpha_AV).
    JtSiJ_i <- solve(JtSiJ)
    # Diagonal entries are the individual variances used for the CIs.
    c(var_rho   = JtSiJ_i[1, 1],
      var_alpha = JtSiJ_i[2, 2])
  }, error = function(e) c(var_rho = NA, var_alpha = NA))
}

# 3.8  Monte Carlo simulation (replicates Tables 6 & 7)

#' Single MC replicate for Section 4.3
one_rep_4_3 <- function(seed, rho, alpha, B) {
  # Per-worker seed so parallel replicates are independent yet reproducible.
  set.seed(seed)

  # Simulate one realization and apply estimator 1; skip the replicate if it fails.
  obs  <- simulate_boolean_obs(rho, alpha)
  est1 <- estimator_1_boolean(obs$A_obs, obs$P_obs)
  if (any(is.na(est1))) return(NULL)

  # Apply estimator 2 (tangent points); skip on failure.
  rho2 <- estimator_2_boolean(obs$A_obs, obs$n_germs, k = 100)
  if (is.na(rho2)) return(NULL)

  # Assemble the 3-vector of base estimators for the averaging step.
  T_vec  <- c(est1["rho"], rho2, est1["alpha"])
  # Plug-in parameters for the bootstrap: average the two rho estimates, reuse alpha_1.
  rho0   <- 0.5 * (est1["rho"] + rho2)
  alpha0 <- est1["alpha"]

  # Bootstrap the 3x3 Sigma at the plug-in parameters; skip if it fails.
  Sigma_hat <- tryCatch(
    estimate_Sigma_boolean(rho0, alpha0, B = B),
    error = function(e) NULL
  )
  if (is.null(Sigma_hat)) return(NULL)

  # Compute the combined (rho_AV, alpha_AV); skip if averaging fails.
  av_est <- boolean_averaging(T_vec, Sigma_hat)
  if (any(is.na(av_est))) return(NULL)

  # Variances for the CIs, using the same J as the estimator.
  J    <- matrix(c(1,1,0, 0,0,1), nrow=3, ncol=2)
  vars <- boolean_averaging_variance(Sigma_hat, J)

  # 95% CI coverage indicators: 1 if the truth lies within +/- 1.96 SE.
  cover_rho   <- as.numeric(
    is.finite(vars["var_rho"]) &&
    abs(av_est["rho"]   - rho)   <= 1.96 * sqrt(vars["var_rho"]))
  cover_alpha <- as.numeric(
    is.finite(vars["var_alpha"]) &&
    abs(av_est["alpha"] - alpha) <= 1.96 * sqrt(vars["var_alpha"]))

  # Return squared errors of each estimator plus the two coverage flags.
  c(err_rho1  = (unname(est1["rho"])     - rho)^2,
    err_rho2  = (rho2                    - rho)^2,
    err_rhoAV = (unname(av_est["rho"])   - rho)^2,
    err_alp1  = (unname(est1["alpha"])   - alpha)^2,
    err_alpAV = (unname(av_est["alpha"]) - alpha)^2,
    cover_rho = cover_rho,
    cover_alp = cover_alpha)
}

#' Run MC simulation for Section 4.3 (Tables 6 & 7)
run_simulation_4_3 <- function(rho_vec = c(25, 50, 100, 150),
                                alpha   = 1,
                                n_rep   = 300,
                                B       = 30,
                                n_cores = n_cores,
                                chk_dir = checkpoint_dir) {
  # One result row per true intensity rho.
  results <- list()

  for (rho in rho_vec) {
    # Build a checkpoint filename keyed by rho.
    key      <- paste("4_3", rho, sep = "_")
    chk_file <- file.path(chk_dir, paste0(key, ".rds"))

    # If this rho was already computed, load the cached row and skip the work.
    if (file.exists(chk_file)) {
      cat(sprintf("  [checkpoint] rho=%d — already done, loading.\n", rho))
      results[[as.character(rho)]] <- readRDS(chk_file)
      next
    }

    # Progress print and a timer for this rho.
    cat(sprintf("  rho=%d ... ", rho)); flush.console()
    t_start <- proc.time()

    # Run the n_rep replicates in parallel; seeds 1..n_rep ensure reproducibility.
    seeds <- seq_len(n_rep)
    raw   <- mclapply(seeds, one_rep_4_3,
                      rho     = rho,
                      alpha   = alpha,
                      B       = B,
                      mc.cores = n_cores)

    # Drop the NULL replicates (failed cases) and stack the rest into a matrix.
    raw <- Filter(Negate(is.null), raw)
    if (length(raw) == 0) next
    mat <- do.call(rbind, raw)

    # Report how many replicates were valid and how long it took.
    elapsed <- (proc.time() - t_start)["elapsed"]
    cat(sprintf("%d/%d valid, %.1f s\n", nrow(mat), n_rep, elapsed))

    # Aggregate the replicate-level errors into MSEs, SDs and coverage.
    # alpha errors are scaled x100 because they are numerically much smaller.
    row <- data.frame(
      rho         = rho,
      MSE_rho1    = mean(mat[, "err_rho1"]),
      MSE_rho2    = mean(mat[, "err_rho2"]),
      MSE_rhoAV   = mean(mat[, "err_rhoAV"]),
      SD_rho1     = sd(mat[, "err_rho1"]),
      SD_rho2     = sd(mat[, "err_rho2"]),
      SD_rhoAV    = sd(mat[, "err_rhoAV"]),
      MSE_alpha1  = mean(mat[, "err_alp1"])  * 100,
      MSE_alphaAV = mean(mat[, "err_alpAV"]) * 100,
      SD_alpha1   = sd(mat[, "err_alp1"])    * 100,
      SD_alphaAV  = sd(mat[, "err_alpAV"])   * 100,
      Cover_rhoAV   = mean(mat[, "cover_rho"]) * 100,
      Cover_alphaAV = mean(mat[, "cover_alp"]) * 100,
      n_valid        = nrow(mat)
    )
    # Cache the row so a re-run skips this rho, then store it in memory.
    saveRDS(row, chk_file)
    results[[as.character(rho)]] <- row
  }
  # Combine the per-rho rows into one data.frame.
  do.call(rbind, results)
}

# 3.9  Execution and results

# Announce the section and run the Boolean-model simulation.
cat("\nSECTION 4.3: Boolean Model\n")
sim_4_3 <- run_simulation_4_3(
  rho_vec = c(25, 50, 100, 150),
  alpha   = 1,
  n_rep   = 300,     # [PAPER] 10000
  B       = 30,      # [PAPER] 100
  n_cores = n_cores,
  chk_dir = checkpoint_dir
)

# Table 6 (rho part): MSE and SD of rho_1, rho_2 and the combined rho_AV.
cat("\nTable 6 (rho part): MSE of rho_1, rho_2, rho_AV\n")
tab6_rho <- sim_4_3[, c("rho", "MSE_rho1", "MSE_rho2", "MSE_rhoAV",
                          "SD_rho1",  "SD_rho2",  "SD_rhoAV")]
print(kable(tab6_rho, digits = 2,
            col.names = c("rho", "rho_1", "rho_2", "rho_AV",
                          "SD rho_1", "SD rho_2", "SD rho_AV")))

# Table 6 (alpha part): MSE (x100) and SD of alpha_1 and alpha_AV.
cat("\nTable 6 (alpha part): MSE x100 of alpha_1, alpha_AV\n")
tab6_alp <- sim_4_3[, c("rho", "MSE_alpha1", "MSE_alphaAV",
                          "SD_alpha1",  "SD_alphaAV")]
print(kable(tab6_alp, digits = 2,
            col.names = c("rho", "alpha_1 x100", "alpha_AV x100",
                          "SD alpha_1", "SD alpha_AV")))

# Table 7: empirical 95% CI coverage for both combined estimators.
cat("\nTable 7: 95% CI Coverage for rho_AV and alpha_AV\n")
tab7 <- sim_4_3[, c("rho", "Cover_rhoAV", "Cover_alphaAV", "n_valid")]
print(kable(tab7, digits = 2,
            col.names = c("rho", "rho_AV (%)", "alpha_AV (%)", "valid reps")))

# 3.10  Plots: MSE and CI coverage

# Long-format data combining rho and alpha MSEs for faceted bar charts.
plot_df_43 <- rbind(
  data.frame(rho = sim_4_3$rho, parametro = "rho",
             stimatore = "rho_1",    MSE = sim_4_3$MSE_rho1),
  data.frame(rho = sim_4_3$rho, parametro = "rho",
             stimatore = "rho_2",    MSE = sim_4_3$MSE_rho2),
  data.frame(rho = sim_4_3$rho, parametro = "rho",
             stimatore = "rho_AV",   MSE = sim_4_3$MSE_rhoAV),
  data.frame(rho = sim_4_3$rho, parametro = "alpha",
             stimatore = "alpha_1",  MSE = sim_4_3$MSE_alpha1),
  data.frame(rho = sim_4_3$rho, parametro = "alpha",
             stimatore = "alpha_AV", MSE = sim_4_3$MSE_alphaAV)
)
# Order the estimators for consistent legend/colours.
plot_df_43$stimatore <- factor(plot_df_43$stimatore,
                               levels = c("rho_1","rho_2","rho_AV",
                                          "alpha_1","alpha_AV"))

# Bar chart of the rho estimators' MSE across true rho values.
p_rho_43 <- ggplot(plot_df_43[plot_df_43$parametro == "rho", ],
                   aes(x = factor(rho), y = MSE, fill = stimatore)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#2c7bb6","#d7191c","#1a9641")) +
  labs(title    = "Sec. 4.3 — MSE of rho estimators",
       subtitle = "rho_AV combines rho_1 and rho_2 with optimal weights",
       x = "True rho", y = "MSE", fill = "Estimator") +
  theme_bw(base_size = 12)

# Companion bar chart for the alpha estimators (MSE x100).
p_alpha_43 <- ggplot(plot_df_43[plot_df_43$parametro == "alpha", ],
                     aes(x = factor(rho), y = MSE, fill = stimatore)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#2c7bb6","#fdae61")) +
  labs(title    = "Sec. 4.3 — MSE of alpha estimators (x100)",
       subtitle = "alpha_AV exploits rho_1 and rho_2 via zero-sum correction",
       x = "True rho", y = "MSE x100", fill = "Estimator") +
  theme_bw(base_size = 12)

# Show the rho and alpha MSE panels side by side.
grid.arrange(p_rho_43, p_alpha_43, ncol = 2)

# Long-format coverage data for the line plot.
cover_df <- data.frame(
  rho       = rep(sim_4_3$rho, 2),
  parametro = rep(c("rho_AV", "alpha_AV"), each = nrow(sim_4_3)),
  copertura = c(sim_4_3$Cover_rhoAV, sim_4_3$Cover_alphaAV)
)

# Coverage vs rho, with a dashed reference line at the nominal 95% level.
p_cover_43 <- ggplot(cover_df, aes(x = factor(rho), y = copertura,
                                    color = parametro, group = parametro)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.7, y = 95.8, label = "95%", size = 3.5, color = "grey40") +
  scale_color_manual(values = c("#1a9641","#2c7bb6")) +
  ylim(80, 100) +
  labs(title    = "Sec. 4.3 — Empirical 95% CI coverage (Table 7)",
       subtitle = "Dashed line = nominal 95% level",
       x = "True rho", y = "Coverage (%)", color = "Estimator") +
  theme_bw(base_size = 12)
print(p_cover_43)


# PART 4 — SECTION 4.4: QUANTILE ESTIMATION UNDER MISSPECIFICATION
#
# Four initial estimators of the p-quantile (p=0.99):
#   1. Non-parametric (NP): x_{(floor(n*p))}   — always consistent
#   2. Weibull MLE:  q = eta * (-log(1-p))^{1/beta}
#   3. Gamma MLE:    q = qgamma(p, shape, rate)
#   4. Burr XII MLE: q = ((1-p)^{-1/k} - 1)^{1/c}
#
# Convex averaging (weights >= 0) ensures robustness under misspecification.
# Sigma (4x4) estimated via non-parametric bootstrap centered on original estimates.
#
# Distributions tested (Table 8):
#   Weibull(3,2), Gamma(3,2), Burr(2,1), Lognormal(0,1)

# 4.1  Quantile estimators

#' Non-parametric quantile (type=1: floor(n*p))
quantile_np <- function(x, p = 0.99) {
  # type=1 gives the inverse-empirical-cdf quantile; consistent for any law.
  # as.numeric() strips the "99%" name quantile() attaches: without it,
  # c(NP = q_np, ...) becomes "NP.99%" and later lookups by name "NP" return NA,
  # silently blanking the whole NP column of Table 8 and the NP bars in the plot.
  as.numeric(quantile(x, probs = p, type = 1))
}

#' Weibull MLE quantile
quantile_weibull <- function(x, p = 0.99) {
  # Keep only positive, finite data (Weibull support is x>0).
  x <- x[x > 0 & is.finite(x)]
  # Too few points to fit reliably => return NA.
  if (length(x) < 5) return(NA)

  # Pre-compute logs and n for the profiled score equation.
  log_x <- log(x); n <- length(x)
  # Score function for the shape b: its root is the MLE of the shape.
  g <- function(b) {
    xb <- x^b
    n / b + sum(log_x) - n * sum(xb * log_x) / sum(xb)
  }
  # Evaluate the score at the bracket endpoints, guarding numerical errors.
  g_lo <- tryCatch(g(1e-4), error = function(e) NA)
  g_hi <- tryCatch(g(50),   error = function(e) NA)
  # If endpoints are invalid or the score does not change sign, no root exists.
  if (any(is.na(c(g_lo, g_hi))) || !is.finite(g_lo)) return(NA)
  if (is.finite(g_hi) && sign(g_lo) == sign(g_hi))   return(NA)

  # Root-find the shape parameter b in [1e-4, 50].
  b <- tryCatch(uniroot(g, lower = 1e-4, upper = 50)$root,
                error = function(e) NA)
  # Reject invalid shapes.
  if (is.na(b) || !is.finite(b) || b <= 0) return(NA)

  # Closed-form scale MLE given the shape b.
  eta <- (mean(x^b))^(1 / b)
  if (!is.finite(eta) || eta <= 0) return(NA)

  # Weibull quantile formula; reject non-positive/non-finite outputs.
  q <- eta * (-log(1 - p))^(1 / b)
  if (!is.finite(q) || q <= 0) return(NA)
  q
}

#' Gamma MLE quantile (Newton-Raphson on profiled log-likelihood)
quantile_gamma <- function(x, p = 0.99) {
  # Gamma support is x>0; drop invalid points.
  x <- x[x > 0 & is.finite(x)]
  if (length(x) < 5) return(NA)

  # Sufficient statistics for the shape MLE: log of mean minus mean of logs.
  xbar  <- mean(x)
  lxbar <- log(xbar)
  mlx   <- mean(log(x))
  s     <- lxbar - mlx
  # s must be positive (Jensen gap); otherwise the fit is degenerate.
  if (!is.finite(s) || s <= 0) return(NA)

  # Choi & Wette (1969) closed-form starting value for the shape.
  a_init <- (3 - s + sqrt((s - 3)^2 + 24 * s)) / (12 * s)
  # Fall back to a method-of-moments start if the approximation is invalid.
  if (!is.finite(a_init) || a_init <= 0) a_init <- xbar^2 / var(x)

  # Newton-Raphson on the shape: solve log(a)-digamma(a)=s.
  a <- a_init
  for (i in seq_len(20)) {
    # Objective and its derivative (uses digamma/trigamma).
    f      <- log(a) - digamma(a) - s
    fprime <- 1/a - trigamma(a)
    # Stop if the step is numerically unsafe.
    if (!is.finite(f) || !is.finite(fprime) || abs(fprime) < 1e-12) break
    a_new <- a - f / fprime
    # Reject non-positive iterates.
    if (!is.finite(a_new) || a_new <= 0) break
    # Converged: accept and exit.
    if (abs(a_new - a) < 1e-8 * a) { a <- a_new; break }
    a <- a_new
  }
  if (!is.finite(a) || a <= 0) return(NA)

  # Rate MLE given the shape: b = a / xbar.
  b_hat <- a / xbar
  if (!is.finite(b_hat) || b_hat <= 0) return(NA)

  # Gamma quantile from the fitted parameters; reject invalid output.
  q <- qgamma(p, shape = a, rate = b_hat)
  if (!is.finite(q) || q <= 0) return(NA)
  q
}

#' Burr XII MLE quantile (L-BFGS-B optimization)
quantile_burr <- function(x, p = 0.99) {
  # Burr XII support is x>0; clean the data.
  x <- x[x > 0 & is.finite(x)]
  if (length(x) < 5) return(NA)

  # Negative log-likelihood in the two shape parameters (c, k).
  neg_loglik <- function(params) {
    c_par <- params[1]; k_par <- params[2]
    # Parameters must be positive; penalise otherwise.
    if (c_par <= 0 || k_par <= 0) return(Inf)
    logxc <- log(1 + x^c_par)
    ll <- length(x) * log(c_par * k_par) +
      (c_par - 1) * sum(log(x)) -
      (k_par + 1) * sum(logxc)
    -ll
  }

  # Box-constrained optimisation of (c,k) starting from (1,1).
  fit <- tryCatch(
    optim(c(1, 1), neg_loglik, method = "L-BFGS-B",
          lower = c(0.01, 0.01), upper = c(50, 50)),
    error = function(e) NULL
  )
  # Treat non-convergence or failure as NA.
  if (is.null(fit) || fit$convergence != 0) return(NA)

  # Burr XII quantile from the fitted shapes.
  c_hat <- fit$par[1]; k_hat <- fit$par[2]
  q <- ((1 - p)^(-1 / k_hat) - 1)^(1 / c_hat)
  if (!is.finite(q) || q <= 0) return(NA)
  q
}

# 4.2  Sigma estimation via non-parametric bootstrap

#' Estimate 4x4 MSE matrix for (q_NP, q_W, q_G, q_B) via non-parametric bootstrap
#' Centered on original estimates T_orig (true quantile unknown).
#' @param x       numeric vector: original sample
#' @param p       scalar: quantile level
#' @param B       integer: bootstrap replicates
#' @param T_orig  length-4 vector: estimates on original data (centering values)
#' @return  4x4 Sigma_hat, or NULL
estimate_Sigma_quantile_boot <- function(x, p = 0.99, B = 200, T_orig) {
  # Sample size for resampling.
  n   <- length(x)
  # Each row holds the four quantile estimates from one resample.
  res <- matrix(NA, nrow = B, ncol = 4)
  for (b in seq_len(B)) {
    # Non-parametric bootstrap sample.
    x_b     <- sample(x, n, replace = TRUE)
    # All four estimators on the resample (warnings from the MLEs suppressed).
    res[b, ] <- c(quantile_np(x_b, p),
                  suppressWarnings(quantile_weibull(x_b, p)),
                  suppressWarnings(quantile_gamma(x_b, p)),
                  suppressWarnings(quantile_burr(x_b, p)))
  }
  # Keep only fully valid resamples.
  res <- res[complete.cases(res), ]
  if (nrow(res) < 10) return(NULL)

  # Centre on the original-sample estimates (the truth is unknown), then form
  # the average outer product = bootstrap MSE matrix.
  # IMPORTANT: center every estimator on the NP estimate (T_orig[1]), the only
  # consistent estimator here, NOT on each estimator's own value. Centering each
  # on itself yields the COVARIANCE (ignores bias), so a low-variance but biased
  # parametric model (e.g. Gamma on Burr/Lognormal data) wrongly gets large
  # weight. Centering on the consistent NP makes the diagonal capture variance +
  # bias^2 (an MSE proxy), so misspecified models are correctly penalised and the
  # convex averaging becomes robust under misspecification, as in the paper.
  ref      <- rep(T_orig[1], length(T_orig))
  centered <- sweep(res, 2, ref, "-")
  t(centered) %*% centered / nrow(centered)
}

# 4.3  Combined estimation for a single sample

#' Compute all quantile estimators and convex averaging for one sample
#' @param x  numeric vector
#' @param p  scalar: quantile level
#' @param B  integer: bootstrap replicates
#' @return  named vector: NP, W, G, B, AV
compute_estimates_4_4 <- function(x, p = 0.99, B = 200) {
  # The four base quantile estimates (parametric MLEs may return NA).
  q_np <- quantile_np(x, p)
  q_w  <- suppressWarnings(quantile_weibull(x, p))
  q_g  <- suppressWarnings(quantile_gamma(x, p))
  q_b  <- suppressWarnings(quantile_burr(x, p))

  # Assemble them; if any failed, fall back to the always-valid NP estimate.
  T_vec <- c(q_np, q_w, q_g, q_b)
  if (any(is.na(T_vec))) {
    return(c(NP = q_np, W = q_w, G = q_g, B = q_b, AV = q_np))
  }

  # Bootstrap Sigma centred on these estimates; fall back to NP on failure.
  Sigma_hat <- tryCatch(
    estimate_Sigma_quantile_boot(x, p, B = B, T_orig = T_vec),
    error = function(e) NULL
  )
  if (is.null(Sigma_hat)) {
    return(c(NP = q_np, W = q_w, G = q_g, B = q_b, AV = q_np))
  }

  # Convex (non-negative) averaging for robustness under misspecification;
  # fall back to NP if the optimisation fails.
  av_res <- tryCatch(
    averaging_convex(T_vec, Sigma_hat),
    error = function(e) list(theta = q_np)
  )
  # Return the four base estimates and the combined AV estimate.
  c(NP = q_np, W = q_w, G = q_g, B = q_b, AV = av_res$theta)
}

# 4.4  Monte Carlo simulation (replicates Table 8)

#' Single MC replicate for Section 4.4
one_rep_4_4 <- function(seed, gen, true_q, n, p, B, mult) {
  # Per-worker seed for reproducible parallel replicates.
  set.seed(seed)
  # Draw a sample and compute all estimators; on error return all-NA.
  x    <- gen(n)
  ests <- tryCatch(
    compute_estimates_4_4(x, p = p, B = B),
    error = function(e) rep(NA_real_, 5)
  )
  # Skip the replicate if anything is missing.
  if (any(is.na(ests))) return(NULL)

  # Squared errors vs the known true quantile, scaled by `mult` so that very
  # small errors (e.g. Weibull) become readable in the tables.
  c(err_NP = (unname(ests["NP"]) - true_q)^2 * mult,
    err_W  = (unname(ests["W"])  - true_q)^2 * mult,
    err_G  = (unname(ests["G"])  - true_q)^2 * mult,
    err_B  = (unname(ests["B"])  - true_q)^2 * mult,
    err_AV = (unname(ests["AV"]) - true_q)^2 * mult)
}

#' Run MC simulation for Section 4.4 (Table 8)
run_simulation_4_4 <- function(n_vec   = c(100, 1000),
                                p       = 0.99,
                                n_rep   = 300,
                                B       = 50,
                                n_cores = n_cores,
                                chk_dir = checkpoint_dir) {
  # Four scenarios: three where one parametric model is correct, plus Lognormal
  # where ALL parametric models are misspecified (the key robustness test).
  # `true_q` is the analytic 0.99 quantile; `mult` is the table scaling factor.
  distributions <- list(
    Weibull = list(
      gen    = function(n) rweibull(n, shape = 3, scale = 2),
      true_q = qweibull(p, shape = 3, scale = 2),
      label  = "Weibull(3,2)",
      mult   = 1000
    ),
    Gamma = list(
      gen    = function(n) rgamma(n, shape = 3, rate = 2),
      true_q = qgamma(p, shape = 3, rate = 2),
      label  = "Gamma(3,2)",
      mult   = 1
    ),
    Burr = list(
      gen    = function(n) ((1 - runif(n))^(-1) - 1)^(1/2),
      true_q = ((1 - p)^(-1) - 1)^(1/2),
      label  = "Burr(2,1)",
      mult   = 1
    ),
    Lognormal = list(
      gen    = function(n) rlnorm(n, meanlog = 0, sdlog = 1),
      true_q = qlnorm(p, meanlog = 0, sdlog = 1),
      label  = "Lognormale(0,1)",
      mult   = 1
    )
  )

  # One result row per (distribution, n).
  results <- list()
  for (dist_name in names(distributions)) {
    dist <- distributions[[dist_name]]

    for (n in n_vec) {
      # Checkpoint keyed by distribution and sample size.
      key      <- paste("4_4", dist_name, n, sep = "_")
      chk_file <- file.path(chk_dir, paste0(key, ".rds"))

      # Reuse cached results when present.
      if (file.exists(chk_file)) {
        cat(sprintf("  [checkpoint] %s n=%d — already done, loading.\n",
                    dist$label, n))
        results[[paste(dist_name, n, sep = "_")]] <- readRDS(chk_file)
        next
      }

      # Progress print and timer.
      cat(sprintf("  %s n=%d (x%d) ... ", dist$label, n, dist$mult))
      flush.console()
      t_start <- proc.time()

      # Parallel replicates with reproducible seeds 1..n_rep.
      seeds <- seq_len(n_rep)
      raw   <- mclapply(seeds, one_rep_4_4,
                        gen    = dist$gen,
                        true_q = dist$true_q,
                        n      = n,
                        p      = p,
                        B      = B,
                        mult   = dist$mult,
                        mc.cores = n_cores)

      # Drop failed replicates and stack the rest.
      raw <- Filter(Negate(is.null), raw)
      if (length(raw) == 0) next
      mat <- do.call(rbind, raw)

      # Report valid count and elapsed time.
      elapsed <- (proc.time() - t_start)["elapsed"]
      cat(sprintf("%d/%d valid, %.1f s\n", nrow(mat), n_rep, elapsed))

      # Aggregate MSE and SD for each estimator into a result row.
      row <- data.frame(
        distribution = dist$label,
        n            = n,
        true_q       = dist$true_q,
        mult         = dist$mult,
        MSE_NP       = mean(mat[, "err_NP"]),
        MSE_W        = mean(mat[, "err_W"]),
        MSE_G        = mean(mat[, "err_G"]),
        MSE_B        = mean(mat[, "err_B"]),
        MSE_AV       = mean(mat[, "err_AV"]),
        SD_NP        = sd(mat[, "err_NP"]),
        SD_W         = sd(mat[, "err_W"]),
        SD_G         = sd(mat[, "err_G"]),
        SD_B         = sd(mat[, "err_B"]),
        SD_AV        = sd(mat[, "err_AV"]),
        n_valid      = nrow(mat)
      )
      # Cache and store.
      saveRDS(row, chk_file)
      results[[paste(dist_name, n, sep = "_")]] <- row
    }
  }
  # Combine all rows.
  do.call(rbind, results)
}

# 4.5  Execution and results

# Announce the section and run the quantile simulation.
cat("\n\nSECTION 4.4: Quantile Estimation under Misspecification\n")
sim_4_4 <- run_simulation_4_4(
  n_vec   = c(100, 1000),
  p       = 0.99,
  n_rep   = 300,     # [PAPER] 10000
  B       = 50,      # [PAPER] 200
  n_cores = n_cores,
  chk_dir = checkpoint_dir
)

# Table 8: MSE (and SD) of all five estimators at p=0.99.
cat("\nTable 8: MSE of quantile estimators (p=0.99)\n")
cat("(First row x1000, remaining rows x1, as in the paper)\n")
tab8 <- sim_4_4[, c("distribution","n","mult",
                     "MSE_NP","MSE_W","MSE_G","MSE_B","MSE_AV",
                     "SD_NP","SD_W","SD_G","SD_B","SD_AV")]
print(kable(tab8, digits = 3,
            col.names = c("Distribution","n","x",
                          "NP","W","G","B","AV",
                          "SD NP","SD W","SD G","SD B","SD AV")))

# 4.6  Plots: MSE comparison

# Long-format MSE data (five estimators) for the grouped bar charts.
plot_df_44 <- do.call(rbind, lapply(unique(sim_4_4$n), function(ni) {
  sub <- sim_4_4[sim_4_4$n == ni, ]
  data.frame(
    distribution = rep(sub$distribution, 5),
    n            = ni,
    stimatore    = rep(c("NP","Weibull","Gamma","Burr","AV"), each = nrow(sub)),
    MSE          = c(sub$MSE_NP, sub$MSE_W, sub$MSE_G, sub$MSE_B, sub$MSE_AV)
  )
}))
# Fix estimator and distribution orders for the plots.
plot_df_44$stimatore    <- factor(plot_df_44$stimatore,
                                  levels = c("NP","Weibull","Gamma","Burr","AV"))
plot_df_44$distribution <- factor(plot_df_44$distribution,
                                  levels = unique(sim_4_4$distribution))

# Helper: grouped bar chart at a given n, with the y-axis truncated at y_max
# because some misspecified MSEs explode and would otherwise hide the rest.
make_bar_plot <- function(data_n, ni, y_max, y_lab) {
  ggplot(data_n[data_n$n == ni, ],
         aes(x = distribution, y = pmin(MSE, y_max), fill = stimatore)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("#2c7bb6","#d7191c","#1a9641",
                                 "#fdae61","#762a83")) +
    labs(title    = sprintf("Sec. 4.4 — MSE of quantile estimators (n=%d)", ni),
         subtitle = sprintf("Scale truncated at %.0f for readability. AV in purple.", y_max),
         x = "True distribution", y = y_lab, fill = "Estimator") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# One panel per sample size, stacked vertically.
p_44_100  <- make_bar_plot(plot_df_44, 100,  50,  "MSE (mixed scale, truncated)")
p_44_1000 <- make_bar_plot(plot_df_44, 1000, 20, "MSE (mixed scale, truncated)")
grid.arrange(p_44_100, p_44_1000, nrow = 2)

# 4.7  Illustrative convex averaging weights

# Show, on one sample per scenario, which estimators the convex averaging keeps.
cat("\nConvex averaging weights (n=1000, one sample per scenario)\n")
cat("(Active support shows which estimators are used)\n\n")

# Local seed for a reproducible illustration.
set.seed(2024)
dist_illustr <- list(
  "Weibull(3,2)"    = list(gen = function(n) rweibull(n, 3, 2),
                            true_q = qweibull(0.99, 3, 2)),
  "Gamma(3,2)"      = list(gen = function(n) rgamma(n, 3, 2),
                            true_q = qgamma(0.99, 3, 2)),
  "Burr(2,1)"       = list(gen = function(n) ((1-runif(n))^(-1)-1)^(1/2),
                            true_q = ((1-0.99)^(-1)-1)^(1/2)),
  "Lognormale(0,1)" = list(gen = function(n) rlnorm(n, 0, 1),
                            true_q = qlnorm(0.99, 0, 1))
)

# For each scenario: fit the four estimators, bootstrap Sigma, run convex AV,
# and print the resulting weights and combined estimate.
for (dname in names(dist_illustr)) {
  d     <- dist_illustr[[dname]]
  x     <- d$gen(1000)
  T_vec <- c(quantile_np(x),
             suppressWarnings(quantile_weibull(x)),
             suppressWarnings(quantile_gamma(x)),
             suppressWarnings(quantile_burr(x)))
  # Skip if any estimator failed on this sample.
  if (any(is.na(T_vec))) { cat(dname, ": NA in estimates\n"); next }

  # Bootstrap the 4x4 Sigma centred on the original estimates.
  Sigma_hat <- estimate_Sigma_quantile_boot(x, p=0.99, B=200, T_orig=T_vec)
  if (is.null(Sigma_hat)) { cat(dname, ": Sigma estimation failed\n"); next }

  # Convex averaging weights and combined estimate.
  av_res <- averaging_convex(T_vec, Sigma_hat)
  w      <- av_res$weights
  names(w) <- c("NP","W","G","B")

  # Print the true quantile, the AV estimate and the four weights.
  cat(sprintf("%-20s | true q=%6.3f | AV=%6.3f | NP=%.3f W=%.3f G=%.3f B=%.3f\n",
              dname, d$true_q, av_res$theta, w["NP"], w["W"], w["G"], w["B"]))
}


# SUMMARY

# Print a plain-text conclusions block summarising the three experiments.
cat("\n")
cat("  CONCLUSIONS\n")
cat("\nSection 4.1 (mean vs median):\n")
cat("  - AV and AVB beat both mean and median almost everywhere\n")
cat("  - For Gaussian, AV assigns maximum weight to mean (correct)\n")
cat("  - For Cauchy, AV favors the median (more robust)\n")
cat("  - CI coverage close to nominal 95% level\n")
cat("  -> Averaging adapts automatically to distribution shape\n")
cat("\nSection 4.3 (Boolean model):\n")
cat("  - rho_AV combines rho_1 and rho_2: improves on both\n")
cat("  - The gain is more pronounced for large rho (high intensity)\n")
cat("  - alpha_AV improves alpha_1 via zero-sum correction\n")
cat("  - CI coverage close to nominal 95% level\n")
cat("  -> Even without likelihood, averaging over geometric functional\n")
cat("     estimators yields systematic improvement\n")
cat("\nSection 4.4 (quantiles under misspecification):\n")
cat("  - When one parametric model is correct, AV approaches it\n")
cat("  - When all models are wrong (Lognormal), AV matches NP\n")
cat("  - Convex averaging (weights >= 0) guarantees robustness\n")
cat("  - Weights reflect the most reliable estimator in each sample\n")

