
library(shiny)

set.seed(1)
n = 150
X = matrix(runif(n, -1, 1), ncol = 1)
ftrue = function(x) sin(2*pi*x) + 0.5*cos(8*pi*x)
y = ftrue(X[,1]) + rnorm(n, sd = 0.1)


gauss_kernel <- function(X1, X2, rho) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  if (rho <= 0) stop("rho must be positive.")

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

krr_fit <- function(X, y, lambda, rho) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)

  if (length(y) != n) stop("length(y) must equal nrow(X).")
  if (lambda < 0)     stop("lambda must be non-negative.")
  if (rho <= 0)       stop("rho must be positive.")

  K <- gauss_kernel(X, X, rho)
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

predict.krr <- function(object, newdata, ...) {
  newdata <- as.matrix(newdata)

  if (ncol(newdata) != ncol(object$X)) {
    stop("ncol(newdata) must match training X.")
  }

  K_new <- gauss_kernel(newdata, object$X, object$rho)
  drop(K_new %*% object$alpha)
}

ui <- fluidPage(
  titlePanel("Problem 3: KRR with Gaussian kernel"),

  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "rho",
        label = "rho (Gaussian kernel bandwidth):",
        min   = 0.1,
        max   = 50,
        value = 5,
        step  = 0.1
      ),
      sliderInput(
        "lambda",
        label = "lambda (penalty parameter):",
        min   = 1e-4,
        max   = 1,
        value = 0.01,
        step  = 1e-4
      )
    ),
    mainPanel(
      plotOutput("krrPlot", height = "450px")
    )
  )
)


server <- function(input, output, session) {

  output$krrPlot <- renderPlot({

    rho    <- input$rho
    lambda <- input$lambda

    fit <- krr_fit(X, y, lambda = lambda, rho = rho)

    x_grid <- seq(-1, 1, length.out = 400)
    X_grid <- matrix(x_grid, ncol = 1)
    y_hat  <- predict(fit, X_grid)

    plot(
      X[, 1], y,
      pch = 16, col = "grey",
      xlab = "x", ylab = "y",
      main = sprintf("rho = %.2f, lambda = %.4f", rho, lambda)
    )
    lines(x_grid, ftrue(x_grid), lwd = 2)
    lines(x_grid, y_hat, col = "red", lwd = 2)

    legend(
      "topleft",
      legend = c("data", "true f", "KRR fit"),
      col    = c("grey", "black", "red"),
      pch    = c(16, NA, NA),
      lty    = c(NA, 1, 1),
      bty    = "n"
    )
  }, res = 96)
}

shinyApp(ui = ui, server = server)

