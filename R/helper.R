#' Internal function: Set the number of cores for parallel computing
#'
#' @keywords internal
#' @param ncores Number of number of cores for parallel computing. Default is NULL.
#' @return Logical
#'
set_cores <- function(ncores = NULL) {
  if (!is.null(ncores)) {
    if (!is.numeric(ncores))
      stop("Please enter valid type - but got ", class(ncores))
    
    defaultNumber <- RcppParallel::defaultNumThreads()
    if (ncores > defaultNumber)
      stop("The input number of cores are invalid - default is ",
           defaultNumber)
    if (ncores < 1)
      stop("The number of cores is not greater than 1 - but got ", ncores)
    tryCatch({
      RcppParallel::setThreadOptions(numThreads = ncores)
    }, error = print)
  }
}

