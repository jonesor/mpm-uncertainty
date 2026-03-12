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

get_compadre_version <- function(default = "6.26.3.0") {
  version <- getOption("mpm.compadre_version", default = Sys.getenv("MPM_COMPADRE_VERSION", unset = default))
  if (!nzchar(version)) {
    version <- default
  }
  return(version)
}

get_compadre_path <- function(corrected = FALSE, version = get_compadre_version()) {
  suffix <- if (corrected) "_Corrected" else ""
  version_tag <- if (grepl("^v", version)) {
    version
  } else {
    paste0("v", version)
  }
  path <- here::here("data", "raw", "compadre", paste0("COMPADRE_", version_tag, suffix, ".RData"))
  return(path)
}

load_compadre <- function(corrected = FALSE, version = get_compadre_version()) {
  path <- get_compadre_path(corrected = corrected, version = version)
  if (!file.exists(path)) {
    stop("COMPADRE file not found: ", path, call. = FALSE)
  }
  return(Rcompadre::cdb_fetch(path))
}

get_compadre_metadata <- function(corrected = FALSE, version = get_compadre_version()) {
  load_compadre(corrected = corrected, version = version)@version
}

get_compadre_version_label <- function(corrected = FALSE, version = get_compadre_version()) {
  meta <- get_compadre_metadata(corrected = corrected, version = version)
  paste0(meta$Database, " version ", meta$Version, " (created ", meta$DateCreated, ")")
}
