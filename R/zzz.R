.onAttach <- function(libname, pkgname) {
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
  }
  
  cat("NOT_CRAN:", Sys.getenv("NOT_CRAN"), "\n")
  cat("RCPP_PARALLEL_BACKEND:", Sys.getenv("RCPP_PARALLEL_BACKEND"), "\n")
  
  packageStartupMessage("Welcome to SpatMCA")
}