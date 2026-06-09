# STATISTICS FOR DATA SCIENCE — PROJECT A.Y. 2025/26
# Paper: Lavancier & Rochet (2015)
#
# FILE 4/4 — SECTION 4.4: quantile estimation under misspecification
# Standalone: sources the shared setup/theory, then runs Section 4.4
# (Table 8, plots) and prints its own paper-vs-script comparison.

# Load the shared setup and averaging theory (packages, n_cores, averaging_*).
source("00_setup_theory.R")


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
                                n_cores = n_cores) {
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
      label  = "Lognormal(0,1)",
      mult   = 1
    )
  )

  # One result row per (distribution, n).
  results <- list()
  for (dist_name in names(distributions)) {
    dist <- distributions[[dist_name]]

    for (n in n_vec) {
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
      # Store the row.
      results[[paste(dist_name, n, sep = "_")]] <- row
    }
  }
  # Combine all rows.
  do.call(rbind, results)
}

# 4.5  Execution and results

# Announce the section and run the quantile simulation.
cat("\n\nSECTION 4.4: Quantile Estimation under Misspecification\n")
# Simulation parameters for Section 4.4 (collected for the comparison printout).
NVEC_44 <- c(100, 1000)   # sample sizes      [PAPER] 100, 1000
P_44    <- 0.99           # quantile level    [PAPER] 0.99
NREP_44 <- 500            # Monte Carlo reps  [PAPER] 10000 (reduced for runtime)
B_44    <- 100            # bootstrap reps    [PAPER] 200   (reduced for runtime)
sim_4_4 <- run_simulation_4_4(
  n_vec   = NVEC_44,
  p       = P_44,
  n_rep   = NREP_44,
  B       = B_44,
  n_cores = n_cores
)

# Table 8: MSE of all five estimators at p=0.99, with the standard deviation in
# parentheses in the same cell (paper layout).
cat("\nTable 8: MSE of quantile estimators (p=0.99, SD in parentheses)\n")
cat("(First row x1000, remaining rows x1, as in the paper)\n")
tab8 <- data.frame(
  Distribution = sim_4_4$distribution,
  n            = sim_4_4$n,
  x            = sim_4_4$mult,
  NP           = mse_sd(sim_4_4$MSE_NP, sim_4_4$SD_NP, 3),
  W            = mse_sd(sim_4_4$MSE_W,  sim_4_4$SD_W,  3),
  G            = mse_sd(sim_4_4$MSE_G,  sim_4_4$SD_G,  3),
  B            = mse_sd(sim_4_4$MSE_B,  sim_4_4$SD_B,  3),
  AV           = mse_sd(sim_4_4$MSE_AV, sim_4_4$SD_AV, 3),
  stringsAsFactors = FALSE
)
print(kable(tab8, align = "l"))

# 4.6  Plots: MSE comparison

# Long-format MSE data (five estimators) for the grouped bar charts.
plot_df_44 <- do.call(rbind, lapply(unique(sim_4_4$n), function(ni) {
  sub <- sim_4_4[sim_4_4$n == ni, ]
  data.frame(
    distribution = rep(sub$distribution, 5),
    n            = ni,
    estimator    = rep(c("NP","Weibull","Gamma","Burr","AV"), each = nrow(sub)),
    MSE          = c(sub$MSE_NP, sub$MSE_W, sub$MSE_G, sub$MSE_B, sub$MSE_AV)
  )
}))
# Fix estimator and distribution orders for the plots.
plot_df_44$estimator    <- factor(plot_df_44$estimator,
                                  levels = c("NP","Weibull","Gamma","Burr","AV"))
plot_df_44$distribution <- factor(plot_df_44$distribution,
                                  levels = unique(sim_4_4$distribution))

# Helper: grouped bar chart at a given n, with the y-axis truncated at y_max
# because some misspecified MSEs explode and would otherwise hide the rest.
make_bar_plot <- function(data_n, ni, y_max, y_lab) {
  ggplot(data_n[data_n$n == ni, ],
         aes(x = distribution, y = pmin(MSE, y_max), fill = estimator)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = estimator_colors(levels(data_n$estimator))) +
    labs(title = sprintf("n = %d", ni),
         x = "True distribution", y = y_lab, fill = "Estimator") +
    theme_pres() +
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
  "Lognormal(0,1)" = list(gen = function(n) rlnorm(n, 0, 1),
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


# 4.8  Note on the parameters vs the paper

cat("\nNote: the only deviation is the number of Monte Carlo replicates (and the\n")
cat("bootstrap B), reduced from the 10000 of the paper to keep the runtime\n")
cat("manageable. All model parameters are identical to the paper.\n")


# 4.9  Conclusions

cat("\n")
cat("  CONCLUSIONS\n")
cat("\nSection 4.4 (quantiles under misspecification):\n")
cat("  - When one parametric model is correct, AV approaches it\n")
cat("  - When all models are wrong (Lognormal), AV matches NP\n")
cat("  - Convex averaging (weights >= 0) guarantees robustness\n")
cat("  - Weights reflect the most reliable estimator in each sample\n")
