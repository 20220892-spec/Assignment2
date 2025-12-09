#' Fit Kernel Ridge Regression
#'
#' @description
#' Fits a kernel ridge regression model using the Gaussian kernel.
#'
#' @param X Numeric matrix of predictors of size \eqn{n \times p}.
#' @param y Numeric response vector of length \eqn{n}.
#' @param lambda Non-negative penalty parameter.
#' @param rho Positive bandwidth parameter of the Gaussian kernel.
#'
#' @details
#' The estimator solves
#' \deqn{\hat{f} = \arg\min_{f \in \mathcal{H}_K}
#' \frac{1}{n} \sum_{i=1}^n (y_i - f(x_i))^2 + \lambda \|f\|_{\mathcal{H}_K}^2.}
#' In the dual form, the fitted function has the representation
#' \eqn{\hat{f}(x) = \sum_{i=1}^n \alpha_i K(x_i, x)}, where
#' \eqn{\alpha = (K + \lambda I_n)^{-1} y}.
#'
#' @return
#' An object of class \code{"krr"} containing
#' \itemize{
#'   \item \code{X}: training design matrix
#'   \item \code{y}: training response
#'   \item \code{lambda}: penalty parameter
#'   \item \code{rho}: bandwidth parameter
#'   \item \code{alpha}: coefficient vector
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' x <- sort(runif(n, -1, 1))
#' X <- matrix(x, ncol = 1)
#' ftrue <- function(z) sin(2 * pi * z)
#' y <- ftrue(x) + rnorm(n, sd = 0.1)
#'
#' fit <- krr_fit(X, y, lambda = 0.01, rho = 5)
#'
#' x_new <- seq(-1, 1, length.out = 100)
#' X_new <- matrix(x_new, ncol = 1)
#' y_hat <- predict(fit, X_new)
#'
#' @export
krr_fit <- function(X, y, lambda, rho) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  n <- nrow(X)

  if (length(y) != n) {
    stop("length(y) must equal nrow(X).")
  }
  if (lambda < 0) {
    stop("lambda must be non-negative.")
  }
  if (rho <= 0) {
    stop("rho must be positive.")
  }

  # Gaussian kernel matrix on training data
  K <- gauss_kernel(X, X, rho)

  # Solve (K + lambda I_n) alpha = y
  alpha <- solve(K + lambda * diag(n), y)

  structure(
    list(
      X      = X,
      y      = y,
      lambda = lambda,
      rho    = rho,
      alpha  = alpha
    ),
    class = "krr"
  )
}
