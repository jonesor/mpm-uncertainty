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
    options(mc.cores = parallel::detectCores())
  }
}
