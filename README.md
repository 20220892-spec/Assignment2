# Assignment2: Kernel Ridge Regression with Gaussian Kernel

This package implements **kernel ridge regression (KRR)** using a Gaussian kernel.

The main exported functions are:

- `gauss_kernel()` : Computes a Gaussian kernel matrix between two sets of points.
- `krr_fit()`      : Fits a kernel ridge regression model.
- `predict.krr()`  : Predicts responses at new data points from a fitted KRR model.

---

## Installation (for this assignment)

This package is developed inside an RStudio project.  
To load the package during development, use:

```r
# from the project root
devtools::load_all()


set.seed(1)

# Generate training data
n <- 100
x <- sort(runif(n, -1, 1))
X <- matrix(x, ncol = 1)

ftrue <- function(z) sin(2 * pi * z)
y <- ftrue(x) + rnorm(n, sd = 0.1)

# Fit KRR with Gaussian kernel
lambda <- 0.01
rho    <- 5

fit <- krr_fit(X, y, lambda = lambda, rho = rho)

# Prediction grid
x_new <- seq(-1, 1, length.out = 200)
X_new <- matrix(x_new, ncol = 1)
y_hat <- predict(fit, newdata = X_new)

# Plot: data, true function, and KRR fit
plot(x, y,
     pch = 16, col = "grey",
     main = sprintf("Kernel Ridge Regression (lambda = %.3f, rho = %.1f)",
                    lambda, rho),
     xlab = "x", ylab = "y")

lines(x_new, ftrue(x_new), lwd = 2)
lines(x_new, y_hat, col = "red", lwd = 2)

legend("topleft",
       legend = c("data", "true f", "KRR fit"),
       col    = c("grey", "black", "red"),
       pch    = c(16, NA, NA),
       lty    = c(NA, 1, 1),
       bty    = "n")

