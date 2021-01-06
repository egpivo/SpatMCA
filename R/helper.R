set_cores <- function(ncores = NULL) {
  if (!is.null(ncores)) {
    if (!is.numeric(ncores))
      stop("Please enter valid type - but got ", class(ncores))
    if (ncores > 1) {
      tryCatch({
        RcppParallel::setThreadOptions(numThreads = ncores)
      }, error = print)
    } else {
      stop("The number of cores is not greater than 1 - but got ", ncores)
    }
  }
}