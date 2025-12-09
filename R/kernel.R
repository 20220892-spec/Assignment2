#' Gaussian kernel matrix
#'
#' @description
#' Compute the Gaussian kernel matrix between two sets of points.
#'
#' @param X1 Numeric matrix of size \eqn{n_1 \times p}.
#' @param X2 Numeric matrix of size \eqn{n_2 \times p}.
#' @param rho Positive bandwidth parameter.
#'
#' @return
#' A numeric matrix of size \eqn{n_1 \times n_2} whose (i, j)-entry is
#' \eqn{\exp\{-\rho \|X1_i - X2_j\|^2\}}.
#'
#' @examples
#' X1 <- matrix(rnorm(10), ncol = 1)
#' X2 <- matrix(rnorm(15), ncol = 1)
#' K  <- gauss_kernel(X1, X2, rho = 1)
#' dim(K)
#'
#' @export
gauss_kernel <- function(X1, X2, rho) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  if (rho <= 0) {
    stop("rho must be positive.")
  }

  n1 <- nrow(X1)
  n2 <- nrow(X2)
  K  <- matrix(0, n1, n2)

  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      diff <- X1[i, ] - X2[j, ]
      K[i, j] <- exp(-rho * sum(diff * diff))
    }
  }
  K
}
