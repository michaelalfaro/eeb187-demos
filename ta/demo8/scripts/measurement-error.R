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
