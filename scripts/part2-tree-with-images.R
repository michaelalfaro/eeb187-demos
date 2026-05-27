#!/usr/bin/env Rscript
# part2-tree-with-images.R
#
# Demo 7 Part 2: read the family's phylogeny + the per-tip exemplar images
# from a FishView bundle, prune the tree to the picked tips, and render the
# tree with each species' image pinned to its tip and the italic species
# name to the right of the image.
#
# Bundle layouts supported:
#   (a) HYBRID:       images/<Genus_species>.png            (flat, preferred)
#                   + images/<Genus_species>/exemplar.png
#                   + images/<Genus_species>/img-N.png      (replicates)
#   (b) flat only:    images/<Genus_species>.png
#   (c) nested only:  images/<Genus_species>/exemplar.png   (legacy)
#
# When BOTH a flat exemplar and a nested exemplar.png exist (the hybrid
# layout), the flat file wins — it is the simpler entry point and is
# guaranteed byte-identical to the nested exemplar.png.
#
# Required files (relative to bundle_dir):
#   <family>-tree.tre               OR  any single *.tre file in the dir
#   picked_tip_labels.txt           one tip label per line (must match tree)
#   images/...                      see layouts above
# Optional:
#   picked_species.txt              canonical "Genus species" names
#
# Output:
#   <out_pdf>  (default results/part2-tree-with-images.pdf)
#
# Usage as a function (sourced from the Quarto chunk):
#   source("scripts/part2-tree-with-images.R")
#   plot_part2(family = "chaetodontidae",
#              bundle_dir = "data/my-fish",
#              out_pdf = "results/part2-tree-with-images.pdf")
#
# Or from the command line:
#   Rscript scripts/part2-tree-with-images.R <family> <bundle_dir> <out_pdf>

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

# Map an image path back to a canonical "Genus_species" key.
# (a) flat:        images/Acanthurus_lineatus.png         -> "Acanthurus_lineatus"
# (b) per-species: images/Acanthurus_lineatus/exemplar.png -> "Acanthurus_lineatus"
.image_key <- function(path) {
  parent <- basename(dirname(path))
  if (grepl("^[A-Z][a-z]+[_ -][a-z]+", parent)) {
    base <- gsub("[- ]+", "_", parent)
  } else {
    base <- tools::file_path_sans_ext(basename(path))
    base <- gsub("[- ]+", "_", base)
  }
  # Capitalize first letter for stable case-insensitive matching downstream
  paste0(toupper(substring(base, 1, 1)), substring(base, 2))
}

plot_part2 <- function(family,
                       bundle_dir   = "data/my-fish",
                       out_pdf      = "results/part2-tree-with-images.pdf",
                       image_height = NULL,   # target image height in tip-row units (1 row = 1 unit)
                       image_aspect = 1.5,    # fallback width/height ratio if image dims fail
                       label_cex    = 0.95) {

  # --- 1. read inputs ------------------------------------------------------
  tree_path <- file.path(bundle_dir, paste0(family, "-tree.tre"))
  if (!file.exists(tree_path)) {
    tre_candidates <- list.files(bundle_dir, pattern = "\\.tre$", full.names = TRUE)
    if (length(tre_candidates) == 0) {
      stop("Could not find any *.tre file in ", bundle_dir)
    }
    tree_path <- tre_candidates[1]
  }

  full_tree <- tryCatch(read.tree(tree_path), error = function(e) NULL)
  if (is.null(full_tree)) full_tree <- read.nexus(tree_path)

  picked_path <- file.path(bundle_dir, "picked_tip_labels.txt")
  if (!file.exists(picked_path)) {
    stop("Missing ", picked_path,
         " - re-export the Part 2 bundle from FishView.")
  }
  picked_tips <- readLines(picked_path)
  picked_tips <- picked_tips[nchar(picked_tips) > 0]

  # --- 2. prune the tree ----------------------------------------------------
  keep <- intersect(picked_tips, full_tree$tip.label)
  if (length(keep) < 2) {
    stop("Fewer than 2 of your picked tips were found in the tree.\n",
         "  picked_tips      : ", paste(picked_tips, collapse = ", "), "\n",
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

  # --- 3. find image files --------------------------------------------------
  img_dir <- file.path(bundle_dir, "images")
  if (!dir.exists(img_dir)) img_dir <- bundle_dir

  img_files <- list.files(img_dir,
                          pattern = "\\.(png|jpe?g)$",
                          ignore.case = TRUE,
                          recursive = TRUE,
                          full.names = TRUE)
  if (length(img_files) == 0) {
    stop("No PNG/JPG images found under ", img_dir)
  }
  keys <- vapply(img_files, .image_key, character(1))

  # Map each tip -> single image path. Preference order:
  #   1. flat exemplar: images/<Genus_species>.<ext> (parent dir = "images")
  #   2. nested exemplar: images/<Genus_species>/exemplar.<ext>
  #   3. any other matching file in the per-species folder
  exemplar_re <- "(^|/)exemplar\\.(png|jpe?g)$"
  flat_re     <- "(^|/)[A-Z][a-z]+_[a-z]+\\.(png|jpe?g)$"
  tip_to_image <- character(length(sub_tree$tip.label))
  names(tip_to_image) <- sub_tree$tip.label
  for (tip in sub_tree$tip.label) {
    m <- regmatches(tip, regexpr("^[A-Z][a-z]+_[a-z]+", tip))
    canonical <- if (length(m) > 0 && nchar(m) > 0) m else tip
    canonical_l <- tolower(canonical)

    hit <- which(tolower(keys) == canonical_l)
    if (length(hit) == 0) hit <- grep(canonical_l, tolower(img_files))

    if (length(hit) > 0) {
      # 1) Flat exemplar wins (file lives directly under images/)
      flat_idx <- which(basename(dirname(img_files[hit])) == "images" &
                        grepl(flat_re, img_files[hit], ignore.case = TRUE))
      if (length(flat_idx) > 0) {
        pick <- hit[flat_idx[1]]
      } else {
        # 2) Nested exemplar.png next
        ex <- grep(exemplar_re, img_files[hit], ignore.case = TRUE)
        # 3) Otherwise first match
        pick <- if (length(ex) > 0) hit[ex[1]] else hit[1]
      }
      tip_to_image[tip] <- img_files[pick]
    }
  }

  missing <- names(tip_to_image)[tip_to_image == ""]
  if (length(missing) > 0) {
    warning("No image found for these tips (skipped on the plot): ",
            paste(missing, collapse = ", "))
  }

  # --- 4. plot --------------------------------------------------------------
  n_tips <- length(sub_tree$tip.label)
  tree_depth <- max(node.depth.edgelength(sub_tree))

  # Auto-scale PDF height so rows don't get squashed for large trees.
  pdf_height <- max(7, 0.45 * n_tips)
  pdf_width  <- 11

  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  pdf(out_pdf, width = pdf_width, height = pdf_height)
  on.exit(dev.off(), add = TRUE)

  # Preload images (so we know the longest label width AND each aspect)
  imgs <- vector("list", n_tips)
  aspects <- numeric(n_tips)     # native width/height
  for (i in seq_len(n_tips)) {
    p <- tip_to_image[i]
    if (is.na(p) || p == "") { aspects[i] <- NA_real_; next }
    im <- tryCatch(.read_image(p), error = function(e) NULL)
    if (is.null(im)) { aspects[i] <- NA_real_; next }
    imgs[[i]]   <- im
    # PNG/JPEG raster arrays are rows x cols (x channels) => height x width
    aspects[i]  <- ncol(im) / nrow(im)
  }
  aspects[!is.finite(aspects)] <- image_aspect

  # Label text per tip (canonical "Genus species", italicized at draw time)
  labels <- vapply(sub_tree$tip.label, function(tip) {
    m <- regmatches(tip, regexpr("^[A-Z][a-z]+_[a-z]+", tip))
    can <- if (length(m) > 0 && nchar(m) > 0) m else tip
    gsub("_", " ", can)
  }, character(1))

  # We want exactly one page. Strategy:
  #   1. Open a blank plot region with the same margins we'll use, so we can
  #      read par("pin") (plot region size in inches) and strwidth() in inches
  #      without committing to an x range.
  #   2. From that, compute target image height in inches, per-image widths
  #      in inches (from native aspect), and the longest label width.
  #   3. Choose x_max in data units so that (tree_depth + reserve) maps to
  #      the full plot-region width in inches.
  #   4. Draw the tree once with that x.lim, then overlay images + labels.
  par(mar = c(2, 1, 3, 2), xpd = NA)
  # Layout pass: open a coordinate system so par("pin") and strwidth() are
  # valid. This consumes the (blank) first page; we then suppress the page
  # advance for the real plot via par(new = TRUE).
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0.5, n_tips + 0.5))
  pin <- par("pin")    # c(width_in, height_in) of the plot region

  row_height     <- 1
  target_h_units <- if (is.null(image_height)) 0.9 * row_height else min(image_height, 0.9 * row_height)
  # y data-units per inch is fixed by the y range we'll use (0.5..n_tips+0.5 = n_tips)
  y_per_in       <- n_tips / pin[2]
  target_h_in    <- target_h_units / y_per_in

  # Per-image widths in inches (from native aspect ratio)
  widths_in  <- target_h_in * aspects
  heights_in <- rep(target_h_in, n_tips)

  # Label widths in inches (italic)
  lbl_widths_in <- vapply(labels, function(s)
    graphics::strwidth(s, units = "inches", cex = label_cex, font = 3),
    numeric(1))

  lbl_gap_in  <- 0.06
  edge_pad_in <- 0.10
  tip_gap_in  <- 0.06   # gap from tree tip to image left edge (inches)

  # Per-tip inches needed to the right of the tree tip
  per_tip_right_in <- widths_in + lbl_gap_in + lbl_widths_in + edge_pad_in
  required_right_in <- max(per_tip_right_in, na.rm = TRUE) + tip_gap_in

  # Tree-depth inches: whatever's left of the plot region
  tree_depth_in <- pin[1] - required_right_in
  if (tree_depth_in < 1.0) {
    # Plot region is too narrow; clamp so tree gets at least 1 inch.
    # In practice this means images/labels will need to overflow into edge_pad,
    # but at least the tree stays drawable.
    tree_depth_in <- 1.0
  }

  # Now set x_per_in so that tree_depth_units occupy tree_depth_in inches.
  x_per_in <- tree_depth / tree_depth_in
  # Total x range in data units = full plot width in inches * x_per_in
  x_max <- pin[1] * x_per_in

  # Draw the tree for real on the SAME page (suppress page advance).
  par(new = TRUE)
  plot(sub_tree,
       show.tip.label = FALSE,
       no.margin      = FALSE,
       x.lim          = c(0, x_max),
       y.lim          = c(0.5, n_tips + 0.5),
       edge.width     = 3,
       edge.color     = "#333333")

  # Title — use ASCII hyphen-minus (double-hyphen) to avoid em-dash mbcsToSbcs
  title(main = sprintf("Part 2: %s -- your exemplars on the tree", family))

  # Tip coords from the most recent plot.phylo
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  # x-left in data units (a small gap past tree depth)
  x_left <- tree_depth + tip_gap_in * x_per_in

  for (i in seq_len(n_tips)) {
    img <- imgs[[i]]
    if (is.null(img)) next

    # Budget for the image alone, in inches, after reserving label + pad
    avail_for_img_in <- pin[1] - tree_depth_in - tip_gap_in -
                        lbl_gap_in - lbl_widths_in[i] - edge_pad_in
    w_in <- widths_in[i]
    h_in <- heights_in[i]
    if (avail_for_img_in > 0 && w_in > avail_for_img_in) {
      shrink <- avail_for_img_in / w_in
      w_in <- avail_for_img_in
      h_in <- h_in * shrink   # preserve aspect
    }

    w_units <- w_in * x_per_in
    h_units <- h_in * y_per_in

    tip_y    <- pp$yy[i]
    y_bottom <- tip_y - h_units / 2
    y_top    <- tip_y + h_units / 2
    x_right  <- x_left + w_units

    rasterImage(img, x_left, y_bottom, x_right, y_top, interpolate = TRUE)

    # Italic label to the right of the image, vertically centered
    lbl_x <- x_right + lbl_gap_in * x_per_in
    text(x = lbl_x, y = tip_y,
         labels = labels[i],
         font = 3, cex = label_cex, adj = c(0, 0.5))
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
