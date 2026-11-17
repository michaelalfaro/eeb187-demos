# =============================================================================
# EE BIOL 138 (Ecology and Evolution of Colour), Lab 7, Wednesday Nov 18, 2026.
# R self-check: run it once before lab.
#
# What it checks
#   1. R is installed on this computer itself (not Docker, not a cloud session).
#   2. The Lab 7 packages load:   ape, phytools, plotrix, png, jpeg
#   3. The project packages load: pavo, caper, ggplot2, tidyr
#   4. R can do what the lab asks: draw a tree, map a trait on it, test
#      phylogenetic signal, read a photograph, and measure a two-colour image.
#   5. The lab folder exists and holds the two data files.
#
# What it changes on the computer
#   It installs no packages. It makes one folder, EEB138-lab7, in the home
#   folder (Documents on Windows) and downloads the two Lab 7 data files into
#   it. Set make_lab_folder <- FALSE below to switch that off.
#
# How to run it
#   Open RStudio (or plain R). Then either
#     File > Open File... > this file > click "Source"
#   or paste one line into the Console:
#     source("https://michaelalfaro.github.io/eeb187-demos/scripts/lab07-r-self-check.R")
#   The last lines say PASS or NOT READY, and list exactly what to fix.
#
# Plain ASCII and base R only, so it runs on a brand-new R installation.
# =============================================================================

local({

  # ---- settings the instructor can change ----------------------------------
  lab_pkgs        <- c("ape", "phytools", "plotrix", "png", "jpeg")
  project_pkgs    <- c("pavo", "caper", "ggplot2", "tidyr")
  check_project   <- TRUE    # FALSE = check only the Lab 7 packages
  make_lab_folder <- TRUE    # FALSE = do not create the folder or download data
  min_r           <- "4.3.0" # oldest R that still has ready-built packages
  data_base <- "https://michaelalfaro.github.io/eeb187-demos/data/frederich"
  data_files <- c("acanthuridae-tree.tree", "acanthuridae-color-patterns.csv")
  lab_dir <- file.path(path.expand("~"), "EEB138-lab7")

  # ---- bookkeeping ----------------------------------------------------------
  lines <- character(0)   # everything printed in the report
  fixes <- character(0)   # blocking problems, each with its fix
  notes <- character(0)   # non-blocking remarks
  ok   <- function(msg) lines <<- c(lines, paste("[ OK ]", msg))
  note <- function(msg) { lines <<- c(lines, paste("[NOTE]", msg)); notes <<- c(notes, msg) }
  fix  <- function(msg, how) {
    lines <<- c(lines, paste("[ FIX]", msg))
    fixes <<- c(fixes, paste0(msg, "\n      Fix: ", how))
  }
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  loads <- function(p) isTRUE(tryCatch(quiet(requireNamespace(p, quietly = TRUE)),
                                      error = function(e) FALSE))
  installed <- function(p) nzchar(system.file(package = p))
  ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "?")

  # ---- 1. R itself ----------------------------------------------------------
  sysname <- Sys.info()[["sysname"]]
  os <- switch(sysname, Darwin = "macOS", Windows = "Windows", sysname)
  arch <- R.version$arch
  rver <- as.character(getRversion())

  if (getRversion() < min_r) {
    fix(sprintf("R %s is too old for this lab (need %s or newer).", rver, min_r),
        "install the current R from https://cran.r-project.org, then reinstall the packages.")
  } else {
    ok(sprintf("R %s on %s (%s)", rver, os, arch))
    if (getRversion() < "4.5.0")
      note(sprintf("R %s is older than the versions this lab was tested on (4.5, 4.6). Update R if a package fails.", rver))
  }

  in_docker <- file.exists("/.dockerenv")
  in_server <- identical(Sys.getenv("RSTUDIO_PROGRAM_MODE"), "server")
  if (in_docker || in_server) {
    fix("This R is running in a container or a browser/cloud session, not on your own computer.",
        "install R and RStudio on the laptop itself (prep page, steps 1 and 2) and run this check there.")
  } else if (!sysname %in% c("Darwin", "Windows")) {
    note(paste("This looks like a", sysname, "system. If that is Docker, Binder, Codespaces or Posit Cloud,",
               "install R on the laptop itself. If it is your own Linux laptop, tell the TA."))
  } else {
    ok("R is running on this computer (not Docker, not a cloud session)")
  }

  if (identical(sysname, "Darwin") && identical(arch, "x86_64")) {
    arm <- tryCatch(quiet(system2("sysctl", c("-n", "hw.optional.arm64"), stdout = TRUE, stderr = FALSE)),
                    error = function(e) "0")
    if (length(arm) && identical(arm[1], "1"))
      note("This Mac has an Apple chip but the Intel build of R. It works; it is only slower.")
  }

  lib <- .libPaths()[1]
  if (file.access(lib, 2) != 0) {
    fix(sprintf("R cannot write to its package folder (%s).", lib),
        "run install.packages(\"ape\") and answer Yes when R offers a personal library; then rerun this check.")
  }

  # ---- 2 and 3. packages ----------------------------------------------------
  check_set <- function(pkgs, label) {
    missing <- pkgs[!vapply(pkgs, installed, logical(1))]
    present <- setdiff(pkgs, missing)
    broken  <- present[!vapply(present, loads, logical(1))]
    good    <- setdiff(present, broken)
    if (length(good))
      ok(sprintf("%s: %s", label, paste(sprintf("%s %s", good, vapply(good, ver, "")), collapse = ", ")))
    if (length(missing))
      fix(sprintf("%s not installed: %s", label, paste(missing, collapse = ", ")),
          sprintf("install.packages(c(%s))   -- if R asks about installing \"from sources\", answer No.",
                  paste(sprintf("\"%s\"", missing), collapse = ", ")))
    if (length(broken))
      fix(sprintf("%s installed but will not load: %s", label, paste(broken, collapse = ", ")),
          sprintf("restart R (Session > Restart R), then install.packages(c(%s)) and answer No to \"from sources\".",
                  paste(sprintf("\"%s\"", broken), collapse = ", ")))
    good
  }
  good_lab <- check_set(lab_pkgs, "Lab 7 packages")
  good_proj <- if (check_project) check_set(project_pkgs, "Project packages") else character(0)

  # ---- 4. can R do what the lab asks? ---------------------------------------
  tmp <- tempfile("eeb138-check-"); dir.create(tmp)

  if (all(c("ape", "phytools", "plotrix") %in% good_lab)) {
    tree_ok <- tryCatch({
      set.seed(138)
      tr <- ape::rcoal(12)                       # a practice tree with 12 tips
      x  <- phytools::fastBM(tr)                 # a practice trait
      pdf_file <- file.path(tmp, "tree.pdf")
      grDevices::pdf(pdf_file, width = 6, height = 6)
      on.exit(if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE), add = TRUE)
      quiet(phytools::plotTree(tr, type = "fan", fsize = 0.7))   # needs plotrix
      cm <- quiet(phytools::contMap(tr, x, plot = FALSE))
      quiet(plot(cm, legend = FALSE))
      grDevices::dev.off()
      K <- quiet(phytools::phylosig(tr, x, method = "K", test = TRUE, nsim = 99))
      lam <- quiet(phytools::phylosig(tr, x, method = "lambda", test = TRUE))
      file.exists(pdf_file) && file.size(pdf_file) > 1000 &&
        is.finite(K$K) && is.finite(lam$lambda)
    }, error = function(e) conditionMessage(e))
    if (isTRUE(tree_ok)) {
      ok("Drew a fan tree, mapped a trait on it, and ran Blomberg's K and Pagel's lambda (practice data)")
    } else {
      fix(paste("The tree test failed:", if (is.character(tree_ok)) tree_ok else "no output was produced"),
          "restart R, run install.packages(c(\"ape\", \"phytools\", \"plotrix\")), answer No to \"from sources\", and rerun.")
    }
  }

  if (all(c("png", "jpeg") %in% good_lab)) {
    img_ok <- tryCatch({
      a <- array(0, dim = c(40, 80, 3)); a[, 1:40, 1] <- 1; a[, 41:80, 3] <- 1   # red half, blue half
      p <- file.path(tmp, "two-colour.png"); j <- file.path(tmp, "two-colour.jpg")
      png::writePNG(a, p); jpeg::writeJPEG(a, j, quality = 0.95)
      all(dim(png::readPNG(p))[1:2] == c(40, 80)) && all(dim(jpeg::readJPEG(j))[1:2] == c(40, 80))
    }, error = function(e) conditionMessage(e))
    if (isTRUE(img_ok)) ok("Wrote and read a PNG and a JPEG image")
    else fix(paste("The image test failed:", if (is.character(img_ok)) img_ok else "wrong image size"),
             "run install.packages(c(\"png\", \"jpeg\")), answer No to \"from sources\", and rerun.")
  }

  if (check_project && "pavo" %in% good_proj && "png" %in% good_lab) {
    pavo_ok <- tryCatch({
      a <- array(0, dim = c(60, 120, 3)); a[, 1:60, 1] <- 1; a[, 61:120, 3] <- 1
      p <- file.path(tmp, "pavo-practice.png"); png::writePNG(a, p)
      img <- quiet(pavo::getimg(p))
      cl  <- quiet(pavo::classify(img, kcols = 2, plotnew = FALSE))
      adj <- quiet(pavo::adjacent(cl, xpts = 30, xscale = 120))
      is.data.frame(adj) && "m" %in% names(adj) && is.finite(adj$m)
    }, error = function(e) conditionMessage(e))
    if (isTRUE(pavo_ok)) ok("pavo classified and measured a two-colour practice image")
    else fix(paste("The pavo test failed:", if (is.character(pavo_ok)) pavo_ok else "no pattern statistics came back"),
             "restart R, run install.packages(\"pavo\"), answer No to \"from sources\", and rerun.")
  }

  # ---- 5. the lab folder and the data ---------------------------------------
  if (make_lab_folder) {
    folder_ok <- tryCatch({
      dir.create(file.path(lab_dir, "data", "frederich"), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(lab_dir, "results"), showWarnings = FALSE)
      dir.exists(file.path(lab_dir, "data", "frederich")) && file.access(lab_dir, 2) == 0
    }, error = function(e) FALSE)
    if (!isTRUE(folder_ok)) {
      fix(sprintf("Could not create the lab folder (%s).", lab_dir),
          "make a folder called EEB138-lab7 yourself, anywhere you can find it, and tell the TA where it is.")
    } else {
      old_to <- getOption("timeout"); options(timeout = 30); on.exit(options(timeout = old_to), add = TRUE)
      got <- vapply(data_files, function(f) {
        dest <- file.path(lab_dir, "data", "frederich", f)
        if (file.exists(dest) && file.size(dest) > 1000) return(TRUE)
        isTRUE(tryCatch({
          quiet(utils::download.file(paste(data_base, f, sep = "/"), dest, quiet = TRUE, mode = "wb"))
          file.exists(dest) && file.size(dest) > 1000
        }, error = function(e) FALSE))
      }, logical(1))
      if (all(got) && "ape" %in% good_lab) {
        n_tips <- tryCatch(length(quiet(ape::read.nexus(file.path(lab_dir, "data", "frederich", data_files[1])))$tip.label),
                           error = function(e) NA_integer_)
        n_spp <- tryCatch(nrow(utils::read.csv(file.path(lab_dir, "data", "frederich", data_files[2]), check.names = FALSE)),
                          error = function(e) NA_integer_)
        if (isTRUE(n_tips == 65) && isTRUE(n_spp == 83))
          ok(sprintf("Lab folder ready: %s (tree with 65 tips, pattern table with 83 species)", lab_dir))
        else
          note(sprintf("The data files are in %s but did not read as expected (tips: %s, species: %s). Show this line to the TA.",
                       lab_dir, n_tips, n_spp))
      } else if (all(got)) {
        ok(sprintf("Lab folder ready: %s (data files downloaded; not read because ape is missing)", lab_dir))
      } else {
        note(sprintf("Lab folder made (%s) but the data did not download (no internet?). You can download it in lab.", lab_dir))
      }
    }
  }
  unlink(tmp, recursive = TRUE)

  # ---- report ---------------------------------------------------------------
  all_pkgs <- c(lab_pkgs, if (check_project) project_pkgs)
  n_good   <- length(c(good_lab, good_proj))
  bar <- strrep("-", 70)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M")
  cat("\n", bar, "\nEE BIOL 138  Lab 7 self-check\n", bar, "\n", sep = "")
  cat(paste(lines, collapse = "\n"), "\n", bar, "\n", sep = "")
  if (length(fixes) == 0) {
    cat("RESULT: PASS. This computer is ready for Lab 7.\n")
    if (make_lab_folder) {
      cat("Your lab folder is:  ", gsub("\\\\", "/", lab_dir), "\n", sep = "")
      cat("In lab you will start with:  setwd(\"~/EEB138-lab7\")\n")
    }
    summary_line <- sprintf("PASS | R %s | %s %s | %d of %d packages | %s", rver, os, arch, n_good, length(all_pkgs), stamp)
  } else {
    cat(sprintf("RESULT: NOT READY. %d thing%s to fix:\n\n", length(fixes), if (length(fixes) == 1) "" else "s"))
    for (i in seq_along(fixes)) cat(sprintf("  %d. %s\n\n", i, fixes[i]))
    cat("Do the fixes in order, then run this check again.\n")
    cat("If the same line fails twice, bring the laptop and this output to lab.\n")
    summary_line <- sprintf("NOT READY | R %s | %s %s | %d of %d packages | %s", rver, os, arch, n_good, length(all_pkgs), stamp)
  }
  cat("\nCopy this line into your lab notebook:\n  ", summary_line, "\n", bar, "\n", sep = "")
  invisible(list(pass = length(fixes) == 0, summary = summary_line, fixes = fixes, notes = notes))
})
