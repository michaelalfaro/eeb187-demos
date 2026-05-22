# helpers.R — Shared utilities for EEB 187 Project Recipe Book
# Source this at the top of each recipe: source("../../code/helpers.R")

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(nlme)
  library(ggplot2)
})

# ── theme_poster ─────────────────────────────────────────────────
theme_poster <- theme_minimal(base_size = 20) +
  theme(
    plot.title    = element_text(size = 28, face = "bold", hjust = 0,
                                margin = margin(b = 12)),
    plot.subtitle = element_text(size = 22, color = "gray30",
                                margin = margin(b = 16)),
    axis.title    = element_text(size = 22, face = "bold"),
    axis.text     = element_text(size = 18),
    legend.title  = element_text(size = 20, face = "bold"),
    legend.text   = element_text(size = 18),
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(color = "gray85", linewidth = 0.4),
    panel.border       = element_rect(color = "gray50", fill = NA,
                                      linewidth = 0.5),
    plot.margin        = margin(20, 20, 20, 20)
  )

# ── load_bundle ──────────────────────────────────────────────────
load_bundle <- function(bundle_dir) {
  if (!dir.exists(bundle_dir))
    stop("Bundle directory not found: ", bundle_dir)
  tree_files <- list.files(bundle_dir, pattern = "-mini-tree\\.tre$",
                           full.names = TRUE)
  if (length(tree_files) == 0)
    stop("No *-mini-tree.tre file in ", bundle_dir)
  tree <- tryCatch(ape::read.nexus(tree_files[1]),
                   error = function(e) ape::read.tree(tree_files[1]))
  if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
  sp_file <- file.path(bundle_dir, "picked_species.txt")
  species <- if (file.exists(sp_file)) {
    readLines(sp_file, warn = FALSE)
  } else {
    tree$tip.label
  }
  traits_file <- file.path(bundle_dir, "species_traits.csv")
  traits <- if (file.exists(traits_file)) {
    read.csv(traits_file, stringsAsFactors = FALSE)
  } else {
    NULL
  }
  list(tree = tree, species = species, traits = traits,
       bundle_dir = bundle_dir)
}

# ── match_tree_data ──────────────────────────────────────────────
match_tree_data <- function(tree, trait_df, species_col = "species") {
  tree_spp  <- tree$tip.label
  data_spp  <- trait_df[[species_col]]
  overlap   <- intersect(tree_spp, data_spp)
  dropped_tree <- setdiff(tree_spp, overlap)
  dropped_data <- setdiff(data_spp, overlap)
  if (length(dropped_tree) > 0)
    message(length(dropped_tree), " species in tree but not in data: ",
            paste(head(dropped_tree, 5), collapse = ", "),
            if (length(dropped_tree) > 5) ", ..." else "")
  if (length(dropped_data) > 0)
    message(length(dropped_data), " species in data but not in tree: ",
            paste(head(dropped_data, 5), collapse = ", "),
            if (length(dropped_data) > 5) ", ..." else "")
  message("Overlap: ", length(overlap), " species")
  pruned_tree <- ape::drop.tip(tree, setdiff(tree_spp, overlap))
  pruned_data <- trait_df[trait_df[[species_col]] %in% overlap, ]
  pruned_data <- pruned_data[match(pruned_tree$tip.label,
                                    pruned_data[[species_col]]), ]
  list(tree = pruned_tree, data = pruned_data)
}

# ── pull_fishbase ────────────────────────────────────────────────
pull_fishbase <- function(species_vec, fields, bundle_dir = NULL) {
  species_names <- gsub("_", " ", species_vec)
  tryCatch({
    suppressPackageStartupMessages(library(rfishbase))
    results <- data.frame(species = species_vec, stringsAsFactors = FALSE)
    if ("species" %in% names(fields)) {
      sp <- rfishbase::species(species_names,
                                fields = c("Species", fields[["species"]]))
      sp$species <- gsub(" ", "_", sp$Species)
      results <- merge(results, sp[, c("species", fields[["species"]])],
                       by = "species", all.x = TRUE)
    }
    if ("estimate" %in% names(fields)) {
      est <- rfishbase::estimate(species_names,
                                  fields = c("Species", fields[["estimate"]]))
      est$species <- gsub(" ", "_", est$Species)
      results <- merge(results, est[, c("species", fields[["estimate"]])],
                       by = "species", all.x = TRUE)
    }
    if ("morphometrics" %in% names(fields)) {
      mor <- rfishbase::morphometrics(species_names,
                                       fields = c("Species", fields[["morphometrics"]]))
      mor$species <- gsub(" ", "_", mor$Species)
      mor_agg <- aggregate(. ~ species,
                            data = mor[, c("species", fields[["morphometrics"]])],
                            FUN = function(x) x[1], na.action = na.pass)
      results <- merge(results, mor_agg, by = "species", all.x = TRUE)
    }
    if ("popgrowth" %in% names(fields)) {
      pg <- rfishbase::popgrowth(species_names,
                                  fields = c("Species", fields[["popgrowth"]]))
      pg$species <- gsub(" ", "_", pg$Species)
      pg_agg <- aggregate(. ~ species,
                           data = pg[, c("species", fields[["popgrowth"]])],
                           FUN = function(x) x[1], na.action = na.pass)
      results <- merge(results, pg_agg, by = "species", all.x = TRUE)
    }
    if ("ecology" %in% names(fields)) {
      eco <- rfishbase::ecology(species_names,
                                 fields = c("Species", fields[["ecology"]]))
      eco$species <- gsub(" ", "_", eco$Species)
      eco_agg <- aggregate(. ~ species, data = eco[, c("species", fields[["ecology"]])],
                           FUN = function(x) x[1], na.action = na.pass)
      results <- merge(results, eco_agg, by = "species", all.x = TRUE)
    }
    if ("reproduction" %in% names(fields)) {
      rep <- rfishbase::reproduction(species_names,
                                      fields = c("Species", fields[["reproduction"]]))
      rep$species <- gsub(" ", "_", rep$Species)
      rep_agg <- aggregate(. ~ species, data = rep[, c("species", fields[["reproduction"]])],
                           FUN = function(x) x[1], na.action = na.pass)
      results <- merge(results, rep_agg, by = "species", all.x = TRUE)
    }
    if ("stocks" %in% names(fields)) {
      stk <- rfishbase::stocks(species_names,
                                fields = c("Species", fields[["stocks"]]))
      stk$species <- gsub(" ", "_", stk$Species)
      stk_agg <- aggregate(. ~ species, data = stk[, c("species", fields[["stocks"]])],
                            FUN = function(x) x[1], na.action = na.pass)
      results <- merge(results, stk_agg, by = "species", all.x = TRUE)
    }
    results
  }, error = function(e) {
    message("rfishbase pull failed: ", e$message)
    if (!is.null(bundle_dir)) {
      csv <- file.path(bundle_dir, "species_traits.csv")
      if (file.exists(csv)) {
        message("Falling back to bundled species_traits.csv")
        return(read.csv(csv, stringsAsFactors = FALSE))
      }
    }
    message("No fallback CSV available. Returning empty data frame.")
    data.frame(species = species_vec, stringsAsFactors = FALSE)
  })
}

# ── safe_pgls ────────────────────────────────────────────────────
safe_pgls <- function(formula, df, tree) {
  tryCatch({
    rownames(df) <- df$species
    model <- nlme::gls(formula, data = df,
                       correlation = ape::corBrownian(1, tree,
                                                      form = ~species))
    s <- summary(model)
    cat("PGLS results:\n")
    cat("  Coefficient:", round(coef(model)[2], 4), "\n")
    cat("  t-value:    ", round(s$tTable[2, "t-value"], 3), "\n")
    cat("  p-value:    ", format.pval(s$tTable[2, "p-value"], digits = 3), "\n")
    model
  }, error = function(e) {
    message("PGLS failed: ", e$message)
    NULL
  })
}

# ── format_K_table ───────────────────────────────────────────────
format_K_table <- function(traits_df, tree,
                           vars = c("m", "A", "Jc", "Jt", "m_dS", "m_dL")) {
  results <- data.frame(
    Variable = character(0), K = numeric(0),
    p_value = numeric(0), stringsAsFactors = FALSE
  )
  for (v in vars) {
    if (!(v %in% colnames(traits_df))) {
      message("Skipping ", v, " -- not found in data")
      next
    }
    x <- setNames(traits_df[[v]], traits_df$species)
    x <- x[!is.na(x)]
    overlap <- intersect(names(x), tree$tip.label)
    if (length(overlap) < 5) {
      message("Skipping ", v, " -- fewer than 5 species with data")
      next
    }
    t_pruned <- ape::drop.tip(tree, setdiff(tree$tip.label, overlap))
    x <- x[t_pruned$tip.label]
    res <- phytools::phylosig(t_pruned, x, method = "K", test = TRUE)
    results <- rbind(results, data.frame(
      Variable = v, K = round(res$K, 3),
      p_value = round(res$P, 4), stringsAsFactors = FALSE
    ))
  }
  results
}
