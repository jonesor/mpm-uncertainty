# Setup helpers for reproducible runs

setup_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "),
      "))",
      call. = FALSE
    )
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

setup_rstan <- function() {
  if (requireNamespace("rstan", quietly = TRUE)) {
    rstan::rstan_options(auto_write = TRUE)
    cores <- parallel::detectCores()
    if (is.na(cores) || cores < 1) cores <- 1
    options(mc.cores = cores)
  }
}

setup_prism <- function(prism_dir = "data/raw/prism") {
  if (!requireNamespace("prism", quietly = TRUE)) {
    stop(
      "Missing package: prism. Install with install.packages(\"prism\").",
      call. = FALSE
    )
  }
  if (!dir.exists(prism_dir)) {
    dir.create(prism_dir, recursive = TRUE)
  }
  prism::prism_set_dl_dir(prism_dir)
  invisible(prism_dir)
}
