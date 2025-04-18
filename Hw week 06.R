# 1. Write a function to do root-finding based on Newton-Ralphson method and solve for “sin(x)/x-0.6=0”. 
# (Note: write your own codes, set tol=0.000001, try different initial values)

f <- function(x) sin(x) / x - 0.6
df <- function(x) (cos(x) * x - sin(x)) / x ^ 2

newton <- function(x0, tol = 0.000001, max = 1000) {
  x_old <- x0
  t <- 0
  
  while (TRUE) {
    x_new <- x_old - f(x_old) / df(x_old)
    t <- t + 1
    
    if (abs(x_new - x_old) < tol) {
      cat("stop at", t, "times", "\n")
      return(x_new)
    }
    
    if (t >= max) {
      cat("over maximum time of trying")
      return(NA)
    }
    
    x_old <- x_new
  }
}

newton(1)
newton(2)
newton(3)
newton(4)
newton(5)

# 2. Use data from Vidal (1980) and find the Belehradek’s equation for C2, C3, C4, C5 by minimizing the least square error, and set b=-2.05. 
# Plot the data and fitted curves.

setwd("C:/Users/dan91/Rstudio/statistic-in-ecology/stat_data")
data <- read.table("VidalTvsDuration.txt", header=T, sep='\t')

b <- -2.05

Belehradek_SSE <- function(params, T, y, b){
  a <- params[1]
  T0 <- params[2]
  
  y_hat <- a * (T - T0) ^ b
  sum((y- y_hat) ^ 2)
}

plot(NULL, xlim = c(8, 16), ylim = c(0, 30), xlab = "Temperature", ylab = "stage duration", main = "Calanus pacifics")

T_obs <- data$X.temperature

for (i in 2:ncol(data)) {
  y_obs <- data[, i]
  
  min <- optim(c(1, 5), Belehradek_SSE, T = T_obs, y = y_obs, b = b)
  a <- min$par[1]
  T0 <- min$par[2]
  
  points(T_obs, y_obs, col = "red", pch = 16)
  
  T_seq <- seq(min(T_obs), max(T_obs), length.out = 100)
  y_min <- a * (T_seq - T0)^b
  lines(T_seq, y_min, col = "blue", lwd = 3)
}
