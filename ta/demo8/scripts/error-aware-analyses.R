# error-aware-analyses.R
# -----------------------------------------------------------------
# Helper used by Demo 8 Part 3.
# Monte Carlo: assume every tip's pavo values are noisy with a per-
# variable SD (from measurement-error.R). Perturb every tip with
# Gaussian noise on each iteration and re-run Parts 1 + 2.
# Returns distributions of:
#   - PC1 loadings
#   - per-variable Blomberg's K and Pagel's lambda
#   - PGLS slopes for the two named bivariate fits
# -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(phytools)
  library(caper)
})

# Run all stats one time on a noise-perturbed trait matrix.
.one_iteration <- function(traits, tree, error_sd,
                           pgls_form_1 = A    ~ m,
                           pgls_form_2 = m_dL ~ m_dS) {

  # Gaussian noise: SD per column = error_sd[col], applied independently to each tip.
  noise <- matrix(rnorm(prod(dim(traits)),
                        sd = rep(error_sd[colnames(traits)],
                                 each = nrow(traits))),
                  nrow = nrow(traits),
                  dimnames = dimnames(traits))

  perturbed <- traits + noise

  # --- PC1 loadings (ahistorical PCA on perturbed data)
  pc_loadings <- tryCatch(prcomp(perturbed, scale. = TRUE)$rotation[, 1],
                          error = function(e) rep(NA_real_, ncol(perturbed)))

  # --- per-variable Blomberg's K and Pagel's lambda (no permutation test in MC)
  sig <- t(sapply(colnames(perturbed), function(v) {
    x <- setNames(perturbed[, v], rownames(perturbed))
    k <- tryCatch(phytools::phylosig(tree, x, method = "K"),
                  error = function(e) NA_real_)
    l <- tryCatch(phytools::phylosig(tree, x, method = "lambda"),
                  error = function(e) list(lambda = NA_real_))
    c(K = as.numeric(k), lambda = as.numeric(l$lambda))
  }))

  # --- PGLS slopes
  pdf <- data.frame(species = rownames(perturbed), perturbed)
  cd  <- caper::comparative.data(phy = tree, data = pdf,
                                 names.col = "species", vcv = TRUE)

  m1 <- tryCatch(caper::pgls(pgls_form_1, data = cd, lambda = "ML"),
                 error = function(e) NULL)
  m2 <- tryCatch(caper::pgls(pgls_form_2, data = cd, lambda = "ML"),
                 error = function(e) NULL)

  slopes <- c(
    slope_1 = if (!is.null(m1)) unname(coef(m1)[2]) else NA_real_,
    slope_2 = if (!is.null(m2)) unname(coef(m2)[2]) else NA_real_
  )

  list(pc_loadings = pc_loadings, signal = sig, pgls_slopes = slopes)
}

# Public entry point. Loops .one_iteration() n_iter times and reduces.
#
# Args:
#   traits   — species × variable trait matrix (rows in tree-tip order)
#   tree     — pruned ape::phylo with rownames(traits) == tree$tip.label
#   error_sd — named numeric vector of per-variable SDs (names must match
#              colnames(traits))
#   n_iter   — number of MC iterations (default 200)
#   seed     — optional integer for reproducibility
#
# Returns a list with raw matrices and summary 95 % CIs.
mc_redo_analyses <- function(traits, tree, error_sd, n_iter = 200,
                             seed = NULL, progress = TRUE) {

  stopifnot(all(colnames(traits) %in% names(error_sd)))
  stopifnot(all(rownames(traits) == tree$tip.label))

  if (!is.null(seed)) set.seed(seed)

  iters <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    iters[[i]] <- .one_iteration(traits, tree, error_sd)
    if (progress && (i %% 25 == 0 || i == n_iter)) {
      message(sprintf("  MC iter %3d / %d", i, n_iter))
    }
  }

  loadings_mat <- do.call(rbind, lapply(iters, `[[`, "pc_loadings"))
  K_mat        <- do.call(rbind, lapply(iters, function(it) it$signal[, "K"]))
  lambda_mat   <- do.call(rbind, lapply(iters, function(it) it$signal[, "lambda"]))
  slopes_mat   <- do.call(rbind, lapply(iters, `[[`, "pgls_slopes"))

  ci <- function(x) stats::quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

  list(
    n_iter      = n_iter,
    error_sd    = error_sd,
    pc1_loadings = loadings_mat,
    K            = K_mat,
    lambda       = lambda_mat,
    pgls_slopes  = slopes_mat,
    summaries = list(
      pc1_ci    = apply(loadings_mat, 2, ci),
      K_ci      = apply(K_mat, 2, ci),
      lambda_ci = apply(lambda_mat, 2, ci),
      slopes_ci = apply(slopes_mat, 2, ci)
    )
  )
}
