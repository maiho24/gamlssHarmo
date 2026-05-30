# Run once after cloning:
#   Rscript install_dependencies.R
#
# This script:
#   1. Installs required R packages.
#   2. Registers the gamlssHarmo command so it is available in your shell.
#
# For conda users (recommended):
#   conda env create -f environment.yml
#   conda activate gamlssHarmo
#   Rscript install_dependencies.R
#   -> gamlssHarmo is immediately available in the active environment.

options(repos = c(CRAN = Sys.getenv("CRAN_MIRROR",
                                     unset = "https://cloud.r-project.org")))

pkgs         <- c("gamlss", "gamlss.dist", "dplyr", "ggplot2",
                  "scales", "RColorBrewer", "yaml", "optparse", "logger")
missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) == 0) {
  message("All packages already installed.")
} else {
  message("Installing: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs)
}

# ---------------------------------------------------------------------------
# Register the gamlssHarmo command
# ---------------------------------------------------------------------------

register_cli <- function() {
  repo_root  <- normalizePath(".")
  launcher   <- file.path(repo_root, "gamlssHarmo")
  is_windows <- .Platform$OS.type == "windows"

  if (!file.exists(launcher))
    stop("Launcher script not found: ", launcher,
         "\nMake sure you run this script from the gamlssHarmo repo root.")

  # ---- Conda environment (highest priority) --------------------------------
  conda_prefix <- Sys.getenv("CONDA_PREFIX")
  if (nzchar(conda_prefix)) {
    if (is_windows) {
      dst <- file.path(conda_prefix, "Scripts", "gamlssHarmo.cmd")
      writeLines(c('@echo off',
                   paste0('Rscript "', launcher, '" %*')),
                 dst)
    } else {
      dst <- file.path(conda_prefix, "bin", "gamlssHarmo")
      unlink(dst)
      file.symlink(launcher, dst)
      Sys.chmod(launcher, "0755")
    }
    message("gamlssHarmo registered at: ", dst)
    return(invisible(dst))
  }

  # ---- Plain R, Linux / macOS: ~/.local/bin --------------------------------
  if (!is_windows) {
    local_bin <- file.path(path.expand("~"), ".local", "bin")
    dir.create(local_bin, recursive = TRUE, showWarnings = FALSE)
    dst <- file.path(local_bin, "gamlssHarmo")
    unlink(dst)
    file.symlink(launcher, dst)
    Sys.chmod(launcher, "0755")
    message("gamlssHarmo registered at: ", dst)
    in_path <- any(grepl(local_bin, strsplit(Sys.getenv("PATH"), ":")[[1]],
                         fixed = TRUE))
    return(invisible(dst))
  }

  # ---- Plain R, Windows: .cmd in repo root ---------------------------------
  dst <- file.path(repo_root, "gamlssHarmo.cmd")
  writeLines(c('@echo off',
               paste0('Rscript "', launcher, '" %*')),
             dst)
  message("gamlssHarmo.cmd created at: ", dst)
  message("\nTo use 'gamlssHarmo' from any directory, add the repo root to PATH:")
  message("  1. Open System Properties -> Advanced -> Environment Variables")
  message("  2. Under 'User variables', edit PATH and append: ", repo_root)
  message("  3. Open a new Command Prompt.")
  message("\nAlternatively, run the pipeline from the repo root using:")
  message("  gamlssHarmo.cmd diagnose --data data\\raw\\my_data.csv")
  return(invisible(dst))
}

register_cli()
