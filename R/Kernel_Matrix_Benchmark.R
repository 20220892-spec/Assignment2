if (!requireNamespace("Rcpp", quietly = TRUE)) {
  install.packages("Rcpp")
}
library(Rcpp)

## (a)

gauss_kernel_R <- function(X, rho) {
  X   <- as.matrix(X)
  rho <- as.numeric(rho)

  if (rho <= 0) stop("rho must be positive.")

  n <- nrow(X)
  p <- ncol(X)

  K <- matrix(0, n, n)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      diff  <- X[i, ] - X[j, ]
      dist2 <- sum(diff * diff)
      K[i, j] <- exp(-rho * dist2)
    }
  }

  K
}

## (b)

cppFunction('
Rcpp::NumericMatrix gauss_kernel_cpp(Rcpp::NumericMatrix X, double rho) {
  int n = X.nrow();
  int p = X.ncol();

  if (rho <= 0) {
    Rcpp::stop("rho must be positive.");
  }

  Rcpp::NumericMatrix K(n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double dist2 = 0.0;
      for (int k = 0; k < p; ++k) {
        double diff = X(i, k) - X(j, k);
        dist2 += diff * diff;
      }
      K(i, j) = std::exp(-rho * dist2);
    }
  }

  return K;
}
')

## (c)

n      <- 1000
p      <- 1
rho    <- 1
nrep   <- 50

time_R    <- numeric(nrep)
time_Rcpp <- numeric(nrep)

for (b in seq_len(nrep)) {
  X <- matrix(runif(n * p, 0, 1), ncol = p)

  time_R[b] <- system.time(
    gauss_kernel_R(X, rho)
  )[["elapsed"]]

  time_Rcpp[b] <- system.time(
    gauss_kernel_cpp(X, rho)
  )[["elapsed"]]
}

## (d)

bench_result <- data.frame(
  method = c("R", "Rcpp"),
  mean_time = c(mean(time_R), mean(time_Rcpp)),
  sd_time   = c(sd(time_R),   sd(time_Rcpp))
)

print(bench_result)
