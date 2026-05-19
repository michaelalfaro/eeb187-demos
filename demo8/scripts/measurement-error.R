# measurement-error.R
# -----------------------------------------------------------------
# Helper used by Demo 8 Parts 1 + 3.
#
# Runs the standard pavo pipeline on a single image and on a stack
# of N images of ONE species; returns the 6 pattern statistics that
# Demo 8 builds the rest of the analysis around:
#
#   m     — transition density (overall pattern complexity)
#   A     — directional aspect ratio = m_r / m_c
#   Jc    — chromatic Jaccard complexity (palette use)
#   Jt    — transition-Jaccard complexity (which color-pairs co-occur)
#   m_dS  — mean chromatic boundary strength (saturation step at edges)
#   m_dL  — mean achromatic boundary strength (luminance step at edges)
#
# m_dS and m_dL require a coldists table; we build one without a full
# visual model by feeding adjacent() opponent-color Euclidean distances
# computed directly from the classified image's RGB centroids. This
# matches the recipe from Alfaro et al. 2019 (ICB) and the Demo 6
# walkthrough — opponent transforms first, then dS / dL as Euclidean
# distances in that 2D + 1D space.
# -----------------------------------------------------------------

library(pavo)

# Six pavo variables Demo 8 keeps from adjacent() output.
.DEFAULT_PAVO_KEEP <- c("m", "A", "Jc", "Jt", "m_dS", "m_dL")


# Locate the tree file inside a bundle by globbing, so the worksheet
# doesn't need to hardcode <family>-mini-tree.tre per group.
bundle_tree_path <- function(bundle) {
  hits <- list.files(bundle, pattern = "-mini-tree\\.tre$", full.names = TRUE)
  if (length(hits) == 0)
    stop("No *-mini-tree.tre file in ", bundle)
  hits[1]
}


# Read a tree file, auto-detecting NEXUS vs Newick. The FishView server
# ships whichever format the source tree was — we've seen both in the
# wild (Chaetodontidae is Newick; Pomacanthidae is NEXUS).
#
# Strategy (defensive — file format isn't always what the header claims):
#   1. Sniff the first non-blank line for a NEXUS or Newick signature.
#   2. Try the preferred reader. If it errors, fall back to the other.
#   3. If both fail, stop with a useful message that includes the first
#      line of the file so the user can see what's actually in there.
#   4. NEXUS files can contain a multiPhylo (multiple trees). Return the
#      first one in that case.
read_tree_auto <- function(path) {
  if (!file.exists(path))
    stop("Tree file does not exist: ", path)
  raw <- readLines(path, n = 20, warn = FALSE)
  if (!length(raw))
    stop("Tree file is empty: ", path)
  first_nonblank <- raw[nzchar(trimws(raw))][1]

  looks_nexus  <- any(grepl("^\\s*#NEXUS",  raw, ignore.case = TRUE))
  looks_newick <- !looks_nexus &&
                  grepl("^\\s*\\(", first_nonblank %||% "")

  try_nexus  <- function() ape::read.nexus(path)
  try_newick <- function() ape::read.tree(path)

  readers <- if (looks_nexus) {
    list(nexus = try_nexus, newick = try_newick)
  } else if (looks_newick) {
    list(newick = try_newick, nexus = try_nexus)
  } else {
    # No clear signature — try both, Newick first.
    list(newick = try_newick, nexus = try_nexus)
  }

  errs <- character(0)
  for (nm in names(readers)) {
    tr <- tryCatch(readers[[nm]](),
                   error = function(e) {
                     errs[[nm]] <<- conditionMessage(e); NULL
                   })
    if (!is.null(tr)) {
      if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
      return(tr)
    }
  }

  stop(sprintf(
    "Could not parse tree at %s as either Newick or NEXUS.\n  first non-blank line: %s\n  errors: %s",
    path,
    if (is.null(first_nonblank)) "<all blank>" else first_nonblank,
    paste(sprintf("%s: %s", names(errs), errs), collapse = "; ")
  ))
}

# Internal: %||% fallback (defined here so we don't depend on rlang).
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a


# Strip specimen/voucher suffixes from tree tip labels so they reduce to
# bare Genus_species binomials.  Examples handled:
#   Pomacanthus_imperatorF490M998   -> Pomacanthus_imperator
#   Chaetodon_auriga2               -> Chaetodon_auriga
#   Chaetodon_auriga_PW648          -> Chaetodon_auriga
#   Acanthurus_olivaceus            -> Acanthurus_olivaceus   (unchanged)
#
# Rule: keep the first two underscore-separated tokens, and within each
# token keep the leading run of lowercase letters (so trailing digit/upper
# specimen codes are dropped from token 2). Genus is left untouched
# (token 1 already preserves capitalization).
to_binomial <- function(tips) {
  vapply(tips, function(t) {
    parts <- strsplit(t, "_", fixed = TRUE)[[1]]
    if (length(parts) < 2) return(parts[1])
    genus   <- parts[1]
    species <- sub("^([a-z]+).*", "\\1", parts[2])
    paste(genus, species, sep = "_")
  }, character(1), USE.NAMES = FALSE)
}


# Apply to_binomial() to a phylo object's tip labels and drop duplicates
# (keeping the first occurrence so the tree structure stays valid).
normalize_tree_tips <- function(tree) {
  new_labels <- to_binomial(tree$tip.label)
  dup <- duplicated(new_labels)
  if (any(dup)) {
    drop <- tree$tip.label[dup]
    tree <- ape::drop.tip(tree, drop)
    new_labels <- to_binomial(tree$tip.label)
  }
  tree$tip.label <- new_labels
  tree
}


# Build a coldists data frame from a classified image's class-mean RGBs.
# Uses opponent-color transforms:
#   col1 = (R - G) / (R + G)    red-green axis
#   col2 = (G - B) / (G + B)    green-blue axis
#   lum  = R + G + B            achromatic axis
# dS is Euclidean distance in (col1, col2). dL is |Δlum|.
#
# Returns a data frame with columns c1, c2, dS, dL (pavo's expected
# names for adjacent(coldists = ...)).
calc_euc_lum_dists <- function(classed) {
  rgb <- attr(classed, "classRGB")
  if (is.null(rgb))
    stop("classify() output has no classRGB attribute")

  R <- rgb[["R"]]; G <- rgb[["G"]]; B <- rgb[["B"]]
  col1 <- (R - G) / (R + G)
  col2 <- (G - B) / (G + B)
  lum  <- R + G + B

  pairs <- utils::combn(seq_len(nrow(rgb)), 2)
  data.frame(
    c1 = pairs[1, ],
    c2 = pairs[2, ],
    dS = sqrt((col1[pairs[1, ]] - col1[pairs[2, ]])^2 +
              (col2[pairs[1, ]] - col2[pairs[2, ]])^2),
    dL = abs(lum[pairs[1, ]] - lum[pairs[2, ]])
  )
}


# Identify the class with the whitest centroid (highest R+G+B sum). For
# FishView-exported bundle images, this is the painted background; pavo
# treats it as one of the k color classes if we don't mask it out.
pick_white_bg <- function(classed) {
  rgb <- attr(classed, "classRGB")
  which.max(rowSums(rgb))
}


# Run the per-image pavo pipeline on a single PNG and return a 1-row
# data.frame with the 6 columns named in `keep`.
#
# Bundle images are painted white outside the fish silhouette, so we
# identify the white class and pass it to adjacent() as bkgID with
# exclude = "background". Without this, m / A / m_dS / m_dL would all
# be dominated by the long uniform white field around each fish.
pavo_pipeline_one_image <- function(img_path,
                                    kcols  = 4,
                                    xpts   = 100,
                                    xscale = 100,
                                    keep   = .DEFAULT_PAVO_KEEP) {
  img     <- getimg(img_path)
  classed <- classify(img, kcols = kcols, plotnew = FALSE)
  cd      <- calc_euc_lum_dists(classed)
  bkg     <- pick_white_bg(classed)
  adj     <- adjacent(classed,
                      xpts     = xpts,
                      xscale   = xscale,
                      coldists = cd,
                      bkgID    = bkg,
                      exclude  = "background")

  missing <- setdiff(keep, colnames(adj))
  if (length(missing))
    stop("adjacent() did not return: ", paste(missing, collapse = ", "))

  adj[, keep, drop = FALSE]
}


# Run the pipeline on the N images of one species and return mean + SD.
#
# Args:
#   bundle_dir    — top-level bundle directory (e.g. "chaetodontidae-mini-bundle")
#   error_species — species name with the multi-image stack
#                   (subdir bundle_dir/images/<error_species>/img-{1..N}.png)
#   n_images      — how many images in the stack (default 5)
#
# Returns a list with:
#   $per_image — n_images × n_vars matrix of raw pavo stats
#   $mean      — named numeric vector of column means
#   $sd        — named numeric vector of column SDs
pavo_error_stats <- function(bundle_dir,
                             error_species,
                             n_images = 5,
                             kcols    = 4,
                             xpts     = 100,
                             xscale   = 100,
                             keep     = .DEFAULT_PAVO_KEEP) {

  dir   <- file.path(bundle_dir, "images", error_species)
  paths <- file.path(dir, sprintf("img-%d.png", seq_len(n_images)))

  missing <- !file.exists(paths)
  if (any(missing)) {
    stop("Missing image(s): ", paste(paths[missing], collapse = ", "))
  }

  per_image <- do.call(rbind, lapply(paths,
                                     pavo_pipeline_one_image,
                                     kcols  = kcols,
                                     xpts   = xpts,
                                     xscale = xscale,
                                     keep   = keep))
  rownames(per_image) <- sprintf("img-%d", seq_len(n_images))

  list(
    per_image = per_image,
    mean      = colMeans(per_image),
    sd        = apply(per_image, 2, sd)
  )
}
