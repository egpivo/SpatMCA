#' Internal function: Set the number of cores for parallel computing
#'
#' @keywords internal
#' @param ncores Number of cores for parallel computing. Default is NULL.
#' @return Logical
#'
setCores <- function(ncores = NULL) {
  if (!is.null(ncores)) {
    if (!is.numeric(ncores)) {
      stop("Please enter valid type - but got ", class(ncores))
    }

    defaultNumber <- RcppParallel::defaultNumThreads()
    if (ncores > defaultNumber) {
      stop("The input number of cores is invalid - default is ", defaultNumber)
    }
    if (ncores < 1) {
      stop("The number of cores is not greater than 1 - but got ", ncores)
    }
    tryCatch(
      {
        RcppParallel::setThreadOptions(numThreads = ncores)
      },
      error = print
    )
  }
}

#'
#' Internal function: Validate input data for a spatpca object
#'
#' @keywords internal
#' @param x1 Location matrix (\eqn{p \times d}) correponding to Y1.
#' @param x2 Location matrix (\eqn{q \times d}) correponding to Y2.
#' @param Y1 Data matrix (\eqn{n \times p}) of the first variable stores the values at \eqn{p} locations with sample size \eqn{n}.
#' @param Y2 Data matrix (\eqn{n \times q}) of the second variable stores the values at \eqn{q} locations with sample size \eqn{n}.
#' @param M Number of folds for cross-validation
#' @return NULL
#'
checkInputData <- function(x1, x2, Y1, Y2, M) {
  if (nrow(x1) != ncol(Y1)) {
    stop("The number of rows of x1 should be equal to the number of columns of Y1.")
  }
  if (nrow(x2) != ncol(Y2)) {
    stop("The number of rows of x2 should be equal to the number of columns of Y2.")
  }
  if (nrow(x1) < 3 || nrow(x2) < 3) {
    stop("Number of locations must be larger than 2.")
  }
  if (ncol(x1) > 3 || ncol(x2) > 3) {
    stop("Dimension of locations must be less 4.")
  }
  if (nrow(Y1) != nrow(Y2)) {
    stop("The numbers of sample sizes of both data should be equal.")
  }
  if (M >= nrow(Y1)) {
    stop("Number of folds must be less than sample size, but got M = ", M)
  }
}

#'
#' Internal function: Detrend Y by column-wise centering
#'
#' @keywords internal
#' @param Y Data matrix
#' @return Detrended data matrix
#'
detrend <- function(Y, is_Y_detrended) {
  if (is_Y_detrended) {
    return(Y - rep(colMeans(Y), rep.int(nrow(Y), ncol(Y))))
  } else {
    return(Y)
  }
}
