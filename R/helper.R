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
#' Internal function: Validate input data for a spatmca object
#'
#' @keywords internal
#' @param x1 Location matrix (\eqn{p \times d}) corresponding to Y1.
#' @param x2 Location matrix (\eqn{q \times d}) corresponding to Y2.
#' @param Y1 Data matrix (\eqn{n \times p}) of the first variable stores the values at \eqn{p} locations with sample size \eqn{n}.
#' @param Y2 Data matrix (\eqn{n \times q}) of the second variable stores the values at \eqn{q} locations with sample size \eqn{n}.
#' @param M Number of folds for cross-validation
#' @return `NULL`
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

#' Internal function: Plot sequentially
#' @keywords internal
#' @param objs Valid ggplot2 objects
#' @return `NULL`
#' 
plot_sequentially <- function(objs) {
  originalPar <- par(no.readonly = TRUE)
  on.exit(par(par(originalPar)))
  par(ask = TRUE)
  for (obj in objs) {
    suppressWarnings(print(obj))
  }
  par(ask = FALSE)
}


#' Internal function: Plot 2D fields for cross validation results 
#' @keywords internal
#' @param cv_data A dataframe contains columns ``u``, ``v``, and ``cv`` 
#' @param variate A character represent the title
#' @return A ggplot object
plot_cv_field <- function(cv_data, variate) {
  default_theme <- theme_classic() +
    theme(
      text = element_text(size = 24),
      plot.title = element_text(hjust = 0.5)
    )

  result <- ggplot(cv_data, aes(x = u, y = v, z = cv, fill = cv)) +
    geom_tile() +
    scale_y_continuous(
      trans = log_trans(),
      breaks = trans_breaks("log", function(x) exp(x)),
      labels = trans_format("log", math_format(e^.x))
    ) +
    scale_x_continuous(
      trans = log_trans(),
      breaks = trans_breaks("log", function(x) exp(x)),
      labels = trans_format("log", math_format(e^.x))
    ) +
    ggtitle(variate) +
    default_theme
  return(result)
}