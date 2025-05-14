age <- c(1.5, 2, 3.3, 4.3, 5.3, 6.3, 7.3, 8.3, 9.3, 10.3, 11.3, 12.3, 14)
length <- c(15.4, 28.03, 42.18, 46.2, 50.23, 54.26, 51.82, 57.27, 56.98, 58.93, 58.55, 59.91, 59.83)

nll <- function(par, age, length){
  Lmax <- par[1]
  K <- par[2]
  t0 <- par[3]
  sigma <- par[4]
  
  L <- Lmax * (1 - exp(-K * (age - t0)))
  L[L <= 0] <- 1e-6
  
  k <- length(length)
  nll_value <- k * log(sigma) + sum((log(length) - log(L)) ^ 2) / (2 * sigma ^ 2)
  
  return(nll_value)
}

init <- c(65, 0.2, -1, 0.1)
lower <- c(50, 0.01, -5, 0.001)
upper <- c(70, 1, 5, 1)

fit <- optim(par = init, fn = nll, age = age, length = length, method = "L-BFGS-B", lower = lower, upper = upper)

cat("NLL = ", fit$value, "\nL∞ = ", fit$par[1], "\nK = ", fit$par[2], "\nt0 = ", fit$par[3], "\nσ ^ 2 = ", fit$par[4] ^ 2)
