#!/usr/bin/env Rscript
cran_repo <- "https://cloud.r-project.org"

install_archive_version <- function(pkg, version) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = cran_repo, dependencies = TRUE)
  }
  remotes::install_version(pkg, version = version, repos = cran_repo,
                           dependencies = TRUE, upgrade = "never")
}

install_required <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  if (getRversion() < "4.4.0" && identical(pkg, "rms")) {
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
      install.packages("Hmisc", repos = cran_repo, dependencies = TRUE)
    }
    install_archive_version("rms", "6.8-1")
  } else if (getRversion() < "4.4.0" && identical(pkg, "haplo.stats")) {
    if (!requireNamespace("arsenal", quietly = TRUE)) {
      install.packages("arsenal", repos = cran_repo, dependencies = TRUE)
    }
    install_required("rms")
    install_archive_version("haplo.stats", "1.9.7")
  } else {
    install.packages(pkg, repos = cran_repo, dependencies = TRUE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required R package could not be installed: %s", pkg))
  }
}

required_pkgs <- c(
  "haplo.stats", "ggplot2", "dplyr", "tidyr", "openxlsx", "readxl",
  "stringr", "scales", "forcats", "patchwork", "RColorBrewer",
  "gridExtra"
)
optional_pkgs <- c("logistf", "survival", "survminer")

cat("Installing required stage-15 R packages...\n")
for (pkg in required_pkgs) {
  install_required(pkg)
}

cat("Installing optional stage-15 R packages when possible...\n")
for (pkg in optional_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch(
      install.packages(pkg, repos = cran_repo, dependencies = TRUE),
      warning = function(w) invisible(NULL),
      error = function(e) invisible(NULL)
    )
  }
}

cat("Stage-15 R dependency bootstrap complete.\n")
