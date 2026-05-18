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


# ---- 1.2 Inventory the bundle --------------------------------------
tree <- read.tree(file.path(bundle, "chaetodontidae-mini-tree.tre"))

picked_species <- readLines(file.path(bundle, "picked_species.txt"))
picked_tips    <- readLines(file.path(bundle, "picked_tip_labels.txt"))

match_tip <- function(sp, tips) {
  hits <- grep(paste0("^", sp, "(?:[_0-9]|$)"), tips, value = TRUE, perl = TRUE)
  if (length(hits) == 0) NA_character_ else hits[1]
}
species_to_tip <- vapply(picked_species, match_tip, character(1), tips = picked_tips)

unmatched <- picked_species[is.na(species_to_tip)]
if (length(unmatched)) {
  message("Dropping ", length(unmatched), " species with no tree tip: ",
          paste(unmatched, collapse = ", "))
}
species <- picked_species[!is.na(species_to_tip)]
species_to_tip <- species_to_tip[!is.na(species_to_tip)]

stopifnot(all(file.exists(file.path(bundle, "images", species))))
stopifnot(all(species_to_tip %in% tree$tip.label))

length(species)
head(species_to_tip)

readLines(file.path(bundle, "README.txt"), n = 25)


# ---- 1.3 Run pavo on each species (slow: ~30-60 s) -----------------
readme <- readLines(file.path(bundle, "README.txt"))
err_line <- grep("^Error-species", readme, value = TRUE)[1]
error_species <- if (is.na(err_line) || grepl("NONE", err_line)) NA_character_ else
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

run_one <- function(sp) {
  img      <- getimg(pick_image(sp))
  classed  <- classify(img, kcols = 4, plotnew = FALSE)
  adj      <- adjacent(classed, xpts = 50, xscale = 1)
  adj$species <- sp
  adj
}

results <- lapply(species, run_one)
adj_df  <- do.call(rbind, results)

dim(adj_df)
colnames(adj_df)


# ---- 1.4 Species × pattern-stat matrix -----------------------------
keep <- c("p_1", "p_2", "p_3", "p_4", "Sc", "m", "m_r", "m_c", "A")
traits <- adj_df[, c("species", keep)]
rownames(traits) <- traits$species
traits$species   <- NULL

head(traits, 4)


# ---- 1.5 PCA -------------------------------------------------------
pc <- prcomp(traits, scale. = TRUE)
summary(pc)

plot(pc, type = "l", main = "Scree: variance per axis")

biplot(pc, cex = 0.75, main = "Chaetodontidae pattern PCA (no tree)")

pc_df <- as.data.frame(pc$x[, 1:2])
pc_df$species <- rownames(pc_df)

ggplot(pc_df, aes(PC1, PC2, label = species)) +
  geom_point(size = 3, alpha = 0.7, color = "#0a66c2") +
  geom_text(vjust = -0.8, size = 3) +
  labs(title = "Pattern axes (no tree)",
       x = sprintf("PC1 (%.1f%%)", summary(pc)$importance[2, 1] * 100),
       y = sprintf("PC2 (%.1f%%)", summary(pc)$importance[2, 2] * 100)) +
  theme_minimal()


# ---- 1.6 Bivariate exploration -------------------------------------
cor_mat <- cor(traits)
round(cor_mat, 2)

pairs(traits[, c("Sc", "m", "m_r", "m_c", "A")],
      lower.panel = panel.smooth, pch = 19, cex = 0.7)

fit_1 <- lm(A ~ m, data = traits)
summary(fit_1)$coefficients
summary(fit_1)$r.squared

fit_2 <- lm(m_c ~ m_r, data = traits)
summary(fit_2)$coefficients
summary(fit_2)$r.squared


# ---- 1.8 Save outputs for Parts 2 & 3 ------------------------------
dir.create("results", showWarnings = FALSE)
saveRDS(traits,         "results/part1-traits.rds")
saveRDS(pc,             "results/part1-pca.rds")
saveRDS(cor_mat,        "results/part1-cor.rds")
saveRDS(species_to_tip, "results/part1-species-to-tip.rds")


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
# Drop algebraically redundant columns: p_4 (Σp_i=1) and m_c (A=m_r/m_c).
traits_phy <- traits[, setdiff(colnames(traits), c("p_4", "m_c"))]
phy_pc <- phyl.pca(tree, traits_phy, method = "lambda", mode = "corr")
summary(phy_pc)
biplot(phy_pc, main = "Chaetodontidae pattern phylo-PCA")


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

pgls_1 <- pgls(A ~ m, data = cd, lambda = "ML")
pgls_2 <- pgls(m_c ~ m_r,   data = cd, lambda = "ML")

summary(pgls_1)
summary(pgls_2)

compare <- data.frame(
  fit = c("A ~ m", "m_c ~ m_r"),
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


# ---- 2.8 Save outputs ----------------------------------------------
saveRDS(phy_pc,  "results/part2-phylopca.rds")
saveRDS(signal,  "results/part2-signal.rds")
saveRDS(list(pgls_1 = pgls_1, pgls_2 = pgls_2),
        "results/part2-pgls.rds")
saveRDS(compare, "results/part2-compare-lm-pgls.rds")


# ============================================================
# Part 3 — Measurement-error propagation
# ============================================================

# ---- 3.1 Run pavo on the 5-image stack -----------------------------
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


# ---- 3.4 Monte Carlo (slow: ~3-5 min) ------------------------------
source("scripts/error-aware-analyses.R")

mc <- mc_redo_analyses(traits, tree, err$sd, n_iter = 200, seed = 2026)


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
