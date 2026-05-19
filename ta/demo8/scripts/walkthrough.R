# Demo 8 — walkthrough.R
# Same chunks as demo8-worksheet.qmd, no prose.
# Paste section by section into the RStudio Console.

# ============================================================
# Part 1 — Standing diversity (ignore the tree)
# ============================================================

# ---- 1.1 Setup -----------------------------------------------------
library(pavo)
library(ape)
library(ggplot2)

set.seed(2026)
bundle <- "chaetodontidae-mini-bundle"   # your own family: e.g. "acanthuridae-mini-bundle"

# Load pavo pipeline + bundle-loading helpers
source("scripts/measurement-error.R")


# ---- 1.2 Inventory the bundle --------------------------------------
tree <- read_tree_auto(bundle_tree_path(bundle))
tree <- normalize_tree_tips(tree)   # strip specimen suffixes -> Genus_species

picked_species <- readLines(file.path(bundle, "picked_species.txt"))
species <- intersect(picked_species, tree$tip.label)
dropped <- setdiff(picked_species, tree$tip.label)
if (length(dropped))
  message("Dropping ", length(dropped),
          " species with no tip on the tree: ",
          paste(dropped, collapse = ", "))
species_to_tip <- setNames(species, species)

stopifnot(all(file.exists(file.path(bundle, "images", species))))
stopifnot(all(species_to_tip %in% tree$tip.label))

length(species)
head(species)

readLines(file.path(bundle, "README.txt"), n = 25)


# ---- 1.3 Run pavo on each species (slow: ~30-60 s) -----------------
readme <- readLines(file.path(bundle, "README.txt"))
err_line <- grep("^Error-species", readme, value = TRUE)[1]
error_species <- if (is.na(err_line) || grepl("\\(none\\)|NONE|none", err_line)) NA_character_ else
                 trimws(sub("\\s*\\(.*", "",
                            sub(".*?:\\s*", "", err_line)))
cat("Error-species in this bundle:",
    if (is.na(error_species)) "NONE" else error_species, "\n")

pick_image <- function(sp) {
  if (!is.na(error_species) && sp == error_species) {
    file.path(bundle, "images", sp, "img-1.png")
  } else {
    file.path(bundle, "images", sp, "exemplar.png")
  }
}

source("scripts/measurement-error.R")   # pavo_pipeline_one_image(), calc_euc_lum_dists()

results <- lapply(species, function(sp) {
  row <- pavo_pipeline_one_image(pick_image(sp))
  row$species <- sp
  row
})
adj_df <- do.call(rbind, results)

dim(adj_df)
colnames(adj_df)


# ---- 1.4 Species × pattern-stat matrix -----------------------------
keep <- c("m", "A", "Jc", "Jt", "m_dS", "m_dL")
traits <- adj_df[, c("species", keep)]
rownames(traits) <- traits$species
traits$species   <- NULL

head(traits, 4)


# ---- 1.5 PCA -------------------------------------------------------
pc <- prcomp(traits, scale. = TRUE)
summary(pc)

plot(pc, type = "l", main = "Scree: variance per axis")

# Image-point biplot helper (shared between Part 1 and Part 2)
abbr_sp <- function(x) {
  vapply(strsplit(x, "_"), function(p) {
    g <- substr(p[1], 1, 4)
    s <- if (length(p) >= 2) substr(p[2], 1, 4) else ""
    paste0(g, "_", s)
  }, character(1))
}

img_path_for <- function(sp_dir, bundle) {
  c1 <- file.path(bundle, "images", sp_dir, "exemplar.png")
  c2 <- file.path(bundle, "images", sp_dir, "img-1.png")
  if (file.exists(c1)) c1 else c2
}

plot_pca_images <- function(scores, loadings = NULL, var_pct, bundle,
                            score_to_dir = NULL, main = "PCA",
                            img_frac = 0.085) {
  scores <- as.matrix(scores[, 1:2, drop = FALSE])
  keys   <- rownames(scores)
  dirs   <- if (is.null(score_to_dir)) keys else score_to_dir[keys]
  labels <- abbr_sp(dirs)
  imgs   <- lapply(dirs, function(d) png::readPNG(img_path_for(d, bundle)))
  xr <- range(scores[, 1]); yr <- range(scores[, 2])
  pad <- 0.18
  xr2 <- xr + c(-pad, pad) * diff(xr)
  yr2 <- yr + c(-pad, pad) * diff(yr)
  plot(NA, xlim = xr2, ylim = yr2, asp = 1,
       xlab = sprintf("PC1  (%.1f%% variance)", var_pct[1]),
       ylab = sprintf("PC2  (%.1f%% variance)", var_pct[2]),
       main = main)
  img_w <- diff(xr2) * img_frac
  for (i in seq_along(imgs)) {
    h_w <- dim(imgs[[i]])[1] / dim(imgs[[i]])[2]
    img_h <- img_w * h_w
    rasterImage(imgs[[i]],
                scores[i, 1] - img_w/2, scores[i, 2] - img_h/2,
                scores[i, 1] + img_w/2, scores[i, 2] + img_h/2)
    text(scores[i, 1], scores[i, 2] - img_h/2 - diff(yr2) * 0.018,
         labels[i], cex = 0.6, family = "mono")
  }
  if (!is.null(loadings)) {
    L <- as.matrix(loadings[, 1:2, drop = FALSE])
    arrow_scale <- 0.85 * min(diff(xr), diff(yr)) / 2 /
                   max(sqrt(L[, 1]^2 + L[, 2]^2))
    arrows(0, 0, L[, 1] * arrow_scale, L[, 2] * arrow_scale,
           length = 0.1, col = "#c14a4a", lwd = 1.5)
    text(L[, 1] * arrow_scale * 1.18, L[, 2] * arrow_scale * 1.18,
         rownames(L), col = "#c14a4a", font = 2, cex = 0.95)
  }
  invisible(NULL)
}

var_pct <- round(100 * summary(pc)$importance[2, 1:2], 1)
plot_pca_images(pc$x, pc$rotation, var_pct, bundle,
                main = "Chaetodontidae pattern PCA · biplot (no tree)")
plot_pca_images(pc$x, NULL, var_pct, bundle,
                main = "Where each species sits on PC1 vs PC2")


# ---- 1.6 Bivariate exploration -------------------------------------
cor_mat <- cor(traits)
round(cor_mat, 2)

pairs(traits[, c("m", "A", "Jc", "Jt", "m_dS", "m_dL")],
      lower.panel = panel.smooth, pch = 19, cex = 0.7)

fit_1 <- lm(A ~ m, data = traits)
summary(fit_1)$coefficients
summary(fit_1)$r.squared

fit_2 <- lm(m_dL ~ m_dS, data = traits)
summary(fit_2)$coefficients
summary(fit_2)$r.squared


# ---- 1.8 Ecology — pavo patterns vs how the fish lives -------------
eco <- read.csv(file.path(bundle, "species_traits.csv"))

combined <- merge(
  data.frame(species = rownames(traits), traits),
  eco, by = "species", all.x = TRUE
)
combined$log_depth  <- log10(combined$DepthRangeDeep)
combined$log_length <- log10(combined$Length)
rownames(combined)  <- combined$species

fit_eco_1 <- lm(m ~ Troph, data = combined)
fit_eco_2 <- lm(A ~ log_depth, data = combined)
summary(fit_eco_1)$coefficients
summary(fit_eco_2)$coefficients

par(mfrow = c(1, 2))
plot(m ~ Troph, data = combined, pch = 19, col = "#0a66c2")
abline(fit_eco_1, col = "#d4a853", lwd = 2)
plot(A ~ log_depth, data = combined, pch = 19, col = "#0a66c2")
abline(fit_eco_2, col = "#d4a853", lwd = 2)
par(mfrow = c(1, 1))

pc_scores <- pc$x[combined$species, "PC1"]
ecology_cor <- sapply(c("Troph", "log_depth", "log_length", "Vulnerability"),
                     function(v) cor(pc_scores, combined[, v], use = "complete.obs"))
round(ecology_cor, 3)


# ---- 1.9 Save outputs for Parts 2 & 3 ------------------------------
dir.create("results", showWarnings = FALSE)
saveRDS(traits,         "results/part1-traits.rds")
saveRDS(pc,             "results/part1-pca.rds")
saveRDS(cor_mat,        "results/part1-cor.rds")
saveRDS(species_to_tip, "results/part1-species-to-tip.rds")
saveRDS(combined,       "results/part1-combined.rds")


# ============================================================
# Part 2 — Phylogenetic correction (bring the tree back)
# ============================================================

library(phytools)
library(caper)

# ---- 2.1 Match the trait matrix to the tree -----------------------
if (!exists("traits"))         traits         <- readRDS("results/part1-traits.rds")
if (!exists("species_to_tip")) species_to_tip <- readRDS("results/part1-species-to-tip.rds")
if (!exists("tree"))           tree           <- read.tree(file.path(bundle, "chaetodontidae-mini-tree.tre"))

rownames(traits) <- species_to_tip[rownames(traits)]
tree   <- keep.tip(tree, rownames(traits))
traits <- traits[tree$tip.label, ]
stopifnot(all(rownames(traits) == tree$tip.label))


# ---- 2.2 phyloPCA --------------------------------------------------
# The 6-var set is full-rank; no need to drop columns.
phy_pc <- phyl.pca(tree, traits, method = "lambda", mode = "corr")
summary(phy_pc)

tip_to_dir  <- setNames(names(species_to_tip), species_to_tip)
phy_var_pct <- round(100 * diag(phy_pc$Eval) / sum(diag(phy_pc$Eval)), 1)
plot_pca_images(
  scores       = phy_pc$S,
  loadings     = phy_pc$L,
  var_pct      = phy_var_pct[1:2],
  bundle       = bundle,
  score_to_dir = tip_to_dir,
  main         = "Chaetodontidae pattern phylo-PCA · biplot"
)


# ---- 2.3 Compare loadings (ahistorical vs phylo) -------------------
ahist_pc <- readRDS("results/part1-pca.rds")

common_vars <- intersect(rownames(ahist_pc$rotation), rownames(phy_pc$L))
load_df <- data.frame(
  variable        = common_vars,
  ahistorical_PC1 = ahist_pc$rotation[common_vars, "PC1"],
  phyloPCA_PC1    = phy_pc$L[common_vars, 1]
)
print(round(load_df[, -1], 3))

ggplot(load_df, aes(ahistorical_PC1, phyloPCA_PC1, label = variable)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, color = "#d4a853") +
  geom_text(vjust = -0.8, size = 3) +
  labs(x = "Loading on PC1 (ahistorical)",
       y = "Loading on PC1 (phyloPCA)") +
  theme_minimal() + coord_equal()


# ---- 2.4 PGLS ------------------------------------------------------
traits_df <- data.frame(species = rownames(traits), traits)
cd <- comparative.data(phy = tree, data = traits_df,
                       names.col = "species", vcv = TRUE)

pgls_1 <- pgls(A    ~ m,     data = cd, lambda = "ML")
pgls_2 <- pgls(m_dL ~ m_dS,  data = cd, lambda = "ML")

summary(pgls_1)
summary(pgls_2)

compare <- data.frame(
  fit = c("A ~ m", "m_dL ~ m_dS"),
  lm_slope    = c(coef(fit_1)[2], coef(fit_2)[2]),
  lm_p        = c(summary(fit_1)$coefficients[2, 4],
                  summary(fit_2)$coefficients[2, 4]),
  pgls_slope  = c(coef(pgls_1)[2], coef(pgls_2)[2]),
  pgls_p      = c(summary(pgls_1)$coefficients[2, 4],
                  summary(pgls_2)$coefficients[2, 4]),
  pgls_lambda = c(summary(pgls_1)$param["lambda"],
                  summary(pgls_2)$param["lambda"])
)
round(compare[, -1], 3)


# ---- 2.5 Per-variable phylogenetic signal --------------------------
signal <- t(sapply(colnames(traits), function(v) {
  x <- setNames(traits[, v], rownames(traits))
  k_res <- phylosig(tree, x, method = "K",      test = TRUE, nsim = 999)
  l_res <- phylosig(tree, x, method = "lambda", test = TRUE)
  c(K = k_res$K, K_p = k_res$P,
    lambda = as.numeric(l_res$lambda), lambda_p = l_res$P)
}))
signal <- as.data.frame(signal)
round(signal, 3)


# ---- 2.6 Compare to Frédérich K_mult = 0.10 (Chaetodontidae) -------
signal$variable <- rownames(signal)
ggplot(signal, aes(reorder(variable, K), K, fill = K_p < 0.05)) +
  geom_col() +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "#d4a853") +
  geom_text(aes(label = sprintf("p=%.2f", K_p)),
            hjust = -0.1, size = 3, color = "grey25") +
  scale_fill_manual(values = c(`TRUE` = "#0a66c2", `FALSE` = "grey70"),
                    name = "p < 0.05") +
  coord_flip(ylim = c(0, max(signal$K) * 1.25)) +
  labs(subtitle = "dashed line = Frédérich 2026 multivariate K = 0.10",
       x = NULL, y = "Blomberg's K") +
  theme_minimal()


# ---- 2.8 Ecology after phylogenetic correction ---------------------
if (!exists("combined")) combined <- readRDS("results/part1-combined.rds")
combined$species <- species_to_tip[combined$species]
combined <- combined[match(tree$tip.label, combined$species), ]

cd_eco <- comparative.data(phy = tree, data = combined,
                           names.col = "species", vcv = TRUE, na.omit = FALSE)
pgls_eco_1 <- pgls(m ~ Troph,     data = cd_eco, lambda = "ML")
pgls_eco_2 <- pgls(A ~ log_depth, data = cd_eco, lambda = "ML")
summary(pgls_eco_1)
summary(pgls_eco_2)

eco_vars <- c("Troph", "log_depth", "log_length", "Vulnerability")
signal_eco <- t(sapply(eco_vars, function(v) {
  x <- setNames(combined[, v], combined$species)
  x <- x[!is.na(x)]
  k_res <- phytools::phylosig(ape::keep.tip(tree, names(x)), x,
                              method = "K", test = TRUE, nsim = 999)
  l_res <- phytools::phylosig(ape::keep.tip(tree, names(x)), x,
                              method = "lambda", test = TRUE)
  c(K = k_res$K, K_p = k_res$P,
    lambda = as.numeric(l_res$lambda), lambda_p = l_res$P)
}))
signal_eco <- as.data.frame(signal_eco)
round(signal_eco, 3)

# Side-by-side K bar: color patterns vs ecology
signal$kind     <- "color pattern"
signal$variable <- rownames(signal)
signal_eco$kind <- "ecology"
signal_eco$variable <- rownames(signal_eco)
side_by_side <- rbind(
  signal[,     c("variable", "K", "K_p", "kind")],
  signal_eco[, c("variable", "K", "K_p", "kind")]
)
ggplot(side_by_side, aes(reorder(variable, K), K, fill = kind)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  scale_fill_manual(values = c(`color pattern` = "#0a66c2",
                               `ecology`       = "#d4a853")) +
  coord_flip() + theme_minimal() +
  labs(x = NULL, y = "Blomberg's K")


# ---- 2.9 Save outputs ----------------------------------------------
saveRDS(phy_pc,  "results/part2-phylopca.rds")
saveRDS(signal,  "results/part2-signal.rds")
saveRDS(list(pgls_1 = pgls_1, pgls_2 = pgls_2,
             pgls_eco_1 = pgls_eco_1, pgls_eco_2 = pgls_eco_2),
        "results/part2-pgls.rds")
saveRDS(compare,    "results/part2-compare-lm-pgls.rds")
saveRDS(signal_eco, "results/part2-signal-ecology.rds")


# ============================================================
# Part 3 — Measurement-error propagation
# ============================================================

# ---- 3.1 Reload Parts 1 + 2 if you've restarted R ------------------
if (!exists("traits"))         traits         <- readRDS("results/part1-traits.rds")
if (!exists("species_to_tip")) species_to_tip <- readRDS("results/part1-species-to-tip.rds")
if (!exists("signal"))         signal         <- readRDS("results/part2-signal.rds")
if (!exists("pgls_1") || !exists("pgls_2")) {
  .pgls_list <- readRDS("results/part2-pgls.rds")
  pgls_1 <- .pgls_list$pgls_1; pgls_2 <- .pgls_list$pgls_2
}
if (!exists("tree")) {
  tree <- ape::read.tree(file.path(bundle, "chaetodontidae-mini-tree.tre"))
  rownames(traits) <- species_to_tip[rownames(traits)]
  tree   <- ape::keep.tip(tree, rownames(traits))
  traits <- traits[tree$tip.label, ]
}

# Parse error-species name from README (re-parsed so Part 3 stands alone).
readme <- readLines(file.path(bundle, "README.txt"))
err_line <- grep("^Error-species", readme, value = TRUE)[1]
error_species <- if (is.na(err_line) || grepl("\\(none\\)|NONE|none", err_line)) NA_character_ else
                 trimws(sub("\\s*\\(.*", "", sub(".*?:\\s*", "", err_line)))
stopifnot(!is.na(error_species))


# ---- 3.1b Run pavo on the 5-image stack ----------------------------
source("scripts/measurement-error.R")

err <- pavo_error_stats(bundle, error_species, n_images = 5)
err$per_image

data.frame(mean = round(err$mean, 4),
           sd   = round(err$sd,   4),
           cv   = round(err$sd / abs(err$mean), 3))


# ---- 3.2 Strip plot of per-image variation -------------------------
long <- as.data.frame(err$per_image)
long$image <- rownames(long)
long <- tidyr::pivot_longer(long, cols = -image,
                            names_to = "variable", values_to = "value")

ggplot(long, aes(variable, value)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7, color = "#d4a853", size = 2.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, color = "#181b2b") +
  labs(x = NULL, y = "pavo value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


# ---- 3.4 Monte Carlo (n_iter=100 ~ 1.5-2 min; 200 ~ 3-5 min) -------
source("scripts/error-aware-analyses.R")

mc <- mc_redo_analyses(traits, tree, err$sd, n_iter = 100, seed = 2026)


# ---- 3.5 K distributions vs observed -------------------------------
K_long <- tidyr::pivot_longer(as.data.frame(mc$K),
                              cols = everything(),
                              names_to = "variable", values_to = "K")
obs <- signal[, c("variable", "K")]
obs$variable      <- factor(obs$variable,     levels = colnames(mc$K))
K_long$variable   <- factor(K_long$variable,  levels = colnames(mc$K))

ggplot(K_long, aes(variable, K)) +
  geom_boxplot(fill = "grey85", color = "grey40", outlier.size = 0.6) +
  geom_point(data = obs, aes(variable, K),
             color = "#d4a853", size = 3) +
  geom_hline(yintercept = 0.10, linetype = "dashed", color = "#0a66c2") +
  coord_flip() +
  labs(x = NULL, y = "K") +
  theme_minimal()


# ---- 3.6 PGLS slope CIs --------------------------------------------
slope_ci <- mc$summaries$slopes_ci
data.frame(
  fit       = c("A ~ m", "m_c ~ m_r"),
  pgls_obs  = c(coef(pgls_1)[2], coef(pgls_2)[2]),
  ci_low    = slope_ci["2.5%",  ],
  ci_median = slope_ci["50%",   ],
  ci_high   = slope_ci["97.5%", ]
)


# ---- 3.8 Save outputs ----------------------------------------------
saveRDS(err, "results/part3-error-stats.rds")
saveRDS(mc,  "results/part3-mc.rds")
