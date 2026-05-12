#!/usr/bin/env Rscript
# part2-tree-with-images.R
#
# Demo 7 Part 2: read 3 exemplar images + the family's phylogeny, prune the
# tree to those 3 species, and render the tree with each species' image
# pinned at its tip.
#
# Inputs (all live in `data/my-fish/` after you unzip the FishView Part 2
# export bundle):
#   data/my-fish/images/<Genus_species>.png   (3 files — your exemplars)
#   data/my-fish/<family>-tree.tre            (the full family Newick tree)
#   data/my-fish/picked_tip_labels.txt        (the exact tree-tip labels to keep)
#   data/my-fish/picked_species.txt           (canonical "Genus species" names)
#
# Output:
#   results/part2-tree-with-images.pdf
#
# Run from the demo's top-level folder:
#   Rscript scripts/part2-tree-with-images.R
#
# Usage as a function (sourced from the Quarto chunk):
#   source("scripts/part2-tree-with-images.R")
#   plot_part2(family = "acanthuridae",
#              bundle_dir = "data/my-fish",
#              out_pdf = "results/part2-tree-with-images.pdf")

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(png)
})

# Robust image loader — handles PNG (and JPEG if jpeg::readJPEG is available)
.read_image <- function(path) {
  if (grepl("\\.png$", path, ignore.case = TRUE)) return(png::readPNG(path))
  if (grepl("\\.jpe?g$", path, ignore.case = TRUE)) {
    if (!requireNamespace("jpeg", quietly = TRUE)) {
      stop("jpeg package required for .jpg images. install.packages('jpeg').")
    }
    return(jpeg::readJPEG(path))
  }
  stop("Unsupported image extension: ", path)
}

# Match an image filename like "Acanthurus_lineatus.png" or "acanthurus-lineatus.jpg"
# back to a canonical "Genus_species" key.
.image_key <- function(filename) {
  base <- tools::file_path_sans_ext(basename(filename))
  base <- gsub("[-\\s]+", "_", base)
  # Capitalize first letter
  paste0(toupper(substring(base, 1, 1)), substring(base, 2))
}

plot_part2 <- function(family,
                       bundle_dir = "data/my-fish",
                       out_pdf    = "results/part2-tree-with-images.pdf",
                       image_height = 0.18,   # in plot-coordinate units (0..1 of tree depth)
                       image_aspect = 1.5,    # width/height ratio for the image box
                       label_cex = 0.95) {

  # --- 1. read inputs ------------------------------------------------------
  tree_path <- file.path(bundle_dir, paste0(family, "-tree.tre"))
  if (!file.exists(tree_path)) {
    # also try without -tree suffix
    tree_path <- list.files(bundle_dir, pattern = "\\.tre$", full.names = TRUE)[1]
  }
  if (is.na(tree_path) || !file.exists(tree_path)) {
    stop("Could not find tree file in ", bundle_dir)
  }

  # Read as newick first, fall back to nexus
  full_tree <- tryCatch(read.tree(tree_path),
                        error = function(e) NULL)
  if (is.null(full_tree)) {
    full_tree <- read.nexus(tree_path)
  }

  picked_path <- file.path(bundle_dir, "picked_tip_labels.txt")
  if (!file.exists(picked_path)) {
    stop("Missing ", picked_path,
         " — re-export the Part 2 bundle from FishView.")
  }
  picked_tips <- readLines(picked_path)
  picked_tips <- picked_tips[nchar(picked_tips) > 0]

  spp_path <- file.path(bundle_dir, "picked_species.txt")
  picked_species <- if (file.exists(spp_path)) {
    readLines(spp_path); gsub("\\s+", "_", readLines(spp_path))
  } else NULL

  # --- 2. prune the tree to those 3 tips -----------------------------------
  keep <- intersect(picked_tips, full_tree$tip.label)
  if (length(keep) < 2) {
    stop("Fewer than 2 of your picked tips were found in the tree.\n",
         "  picked_tips    : ", paste(picked_tips, collapse = ", "), "\n",
         "  first 5 tree tips: ", paste(head(full_tree$tip.label, 5), collapse = ", "))
  }
  sub_tree <- keep.tip(full_tree, keep)

  # --- 2b. require 3 different genera --------------------------------------
  canonical_kept <- regmatches(sub_tree$tip.label,
                               regexpr("^[A-Z][a-z]+_[a-z]+", sub_tree$tip.label))
  genera <- unique(sub("_.*", "", canonical_kept))
  if (length(genera) < 3) {
    stop(sprintf(
      "Your three picks come from only %d genus/genera (%s).\n",
      length(genera), paste(genera, collapse = ", ")),
      "Go back to FishView, pick one species each from THREE different genera, ",
      "then re-export the Part 2 bundle.")
  }

  # --- 3. find the image files --------------------------------------------
  img_dir <- file.path(bundle_dir, "images")
  if (!dir.exists(img_dir)) {
    # try root of bundle
    img_dir <- bundle_dir
  }
  img_files <- list.files(img_dir,
                          pattern = "\\.(png|jpe?g)$",
                          ignore.case = TRUE,
                          full.names = TRUE)
  if (length(img_files) == 0) {
    stop("No PNG/JPG images found under ", img_dir)
  }

  # Map each tip to its image. Try (a) "Genus_species" match against the
  # filename (canonical), (b) substring fallback.
  tip_to_image <- character(length(sub_tree$tip.label))
  names(tip_to_image) <- sub_tree$tip.label
  for (tip in sub_tree$tip.label) {
    # canonical = Genus_species prefix of the tip label
    m <- regmatches(tip, regexpr("^[A-Z][a-z]+_[a-z]+", tip))
    canonical <- if (length(m) > 0 && nchar(m) > 0) m else tip
    canonical_l <- tolower(canonical)

    keys <- vapply(img_files, .image_key, character(1))
    hit <- which(tolower(keys) == canonical_l)
    if (length(hit) == 0) {
      # substring fallback
      hit <- grep(canonical_l, tolower(img_files))
    }
    if (length(hit) > 0) tip_to_image[tip] <- img_files[hit[1]]
  }

  missing <- names(tip_to_image)[tip_to_image == ""]
  if (length(missing) > 0) {
    warning("No image found for these tips (skipped on the plot): ",
            paste(missing, collapse = ", "))
  }

  # --- 4. plot the tree and overlay the images ----------------------------
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  pdf(out_pdf, width = 11, height = 7)
  on.exit(dev.off(), add = TRUE)

  n_tips <- length(sub_tree$tip.label)
  tree_depth <- max(node.depth.edgelength(sub_tree))

  # Reserve the right ~30% of the plotting region for the images. No tip
  # labels on the tree itself — the italic species name lives under each image.
  x_max <- tree_depth * 1.5
  par(mar = c(2, 1, 3, 2), xpd = NA)

  plot(sub_tree,
       show.tip.label = FALSE,
       no.margin = FALSE,
       x.lim = c(0, x_max),
       y.lim = c(0.5, n_tips + 0.5),
       edge.width = 3,
       edge.color = "#333333")
  title(main = sprintf("Part 2: %s — your three exemplars on the tree", family))

  # Pull the tip coordinates set by the last plot.phylo
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  # Each tip gets its own image box. Width is capped so 3+ tips fit
  # vertically without overlapping (image_height is fraction of full y-range).
  img_h_units <- min(image_height * n_tips, 0.85)   # in tip-row units
  img_w_units <- 0.32 * tree_depth                  # in tree-depth units
  # If we have many tips, scale image height down so it fits in a row.
  row_height <- 1                                   # ape uses 1 unit per tip
  img_h_units <- min(img_h_units, 0.9 * row_height * max(1, 3 / n_tips))

  # x_left for all images: just past the tip end of each branch (small gap)
  x_left  <- tree_depth + 0.05 * tree_depth
  x_right <- x_left + img_w_units
  if (x_right > x_max) {
    img_w_units <- (x_max - x_left) * 0.95
    x_right <- x_left + img_w_units
  }

  for (i in seq_along(sub_tree$tip.label)) {
    tip   <- sub_tree$tip.label[i]
    img_p <- tip_to_image[tip]
    if (img_p == "") next

    img <- .read_image(img_p)

    tip_y <- pp$yy[i]
    y_bottom <- tip_y - img_h_units / 2
    y_top    <- tip_y + img_h_units / 2

    rasterImage(img, x_left, y_bottom, x_right, y_top, interpolate = TRUE)

    # Italic species label centered below the image
    canonical <- regmatches(tip, regexpr("^[A-Z][a-z]+_[a-z]+", tip))
    if (length(canonical) == 0 || nchar(canonical) == 0) canonical <- tip
    label_text <- gsub("_", " ", canonical)
    text(x = (x_left + x_right) / 2,
         y = y_bottom - 0.10 * img_h_units,
         labels = label_text,
         font = 3, cex = 0.95, adj = c(0.5, 1))
  }

  invisible(out_pdf)
}

# --- standalone CLI -----------------------------------------------------------
if (sys.nframe() == 0 && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  family <- if (length(args) >= 1) args[1] else "acanthuridae"
  bundle <- if (length(args) >= 2) args[2] else "data/my-fish"
  out    <- if (length(args) >= 3) args[3] else "results/part2-tree-with-images.pdf"
  plot_part2(family = family, bundle_dir = bundle, out_pdf = out)
  cat("Wrote", out, "\n")
}
