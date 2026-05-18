# measurement-error.R
# -----------------------------------------------------------------
# Helper used by Demo 8 Part 3.
# Runs the standard pavo pipeline on a stack of N images of ONE
# species and returns per-variable mean + SD across the stack.
#
# The SD vector is the demo's estimate of "how noisy is a single
# pavo measurement when you swap out the photograph?"
# -----------------------------------------------------------------

library(pavo)

# Default pavo variable set — matches what worksheet Part 1 extracts.
.DEFAULT_PAVO_KEEP <- c("p_1", "p_2", "p_3", "p_4",
                        "Sc", "m", "m_r", "m_c", "A")

# Run the per-image pavo pipeline on a single PNG and return a 1-row
# data.frame of pattern statistics.
pavo_pipeline_one_image <- function(img_path,
                                    kcols = 4,
                                    xpts  = 50,
                                    keep  = .DEFAULT_PAVO_KEEP) {
  img     <- getimg(img_path)
  classed <- classify(img, kcols = kcols, plotnew = FALSE)
  adj     <- adjacent(classed, xpts = xpts, xscale = 1)
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
                             xpts     = 50,
                             keep     = .DEFAULT_PAVO_KEEP) {

  dir   <- file.path(bundle_dir, "images", error_species)
  paths <- file.path(dir, sprintf("img-%d.png", seq_len(n_images)))

  missing <- !file.exists(paths)
  if (any(missing)) {
    stop("Missing image(s): ", paste(paths[missing], collapse = ", "))
  }

  per_image <- do.call(rbind, lapply(paths,
                                     pavo_pipeline_one_image,
                                     kcols = kcols, xpts = xpts, keep = keep))
  rownames(per_image) <- sprintf("img-%d", seq_len(n_images))

  list(
    per_image = per_image,
    mean      = colMeans(per_image),
    sd        = apply(per_image, 2, sd)
  )
}
