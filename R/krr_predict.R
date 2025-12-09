#' Predict method for Kernel Ridge Regression
#'
#' @description
#' Predict responses at new data points from a fitted \code{"krr"} object.
#'
#' @param object A fitted model of class \code{"krr"} returned by \code{krr_fit()}.
#' @param newdata New data matrix with the same number of columns as the
#'   training design matrix used to fit \code{object}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A numeric vector of predictions corresponding to the rows of \code{newdata}.
#'
#' @export
predict.krr <- function(object, newdata, ...) {
  if (missing(newdata)) {
    stop("newdata must be provided.")
  }

  newdata <- as.matrix(newdata)

  if (ncol(newdata) != ncol(object$X)) {
    stop("ncol(newdata) must match ncol of training X.")
  }

  K_new <- gauss_kernel(newdata, object$X, object$rho)
  drop(K_new %*% object$alpha)
}
