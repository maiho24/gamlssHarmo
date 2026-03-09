# Plain R fallback for non-conda users.
# For conda users, run: conda env create -f environment.yml

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
