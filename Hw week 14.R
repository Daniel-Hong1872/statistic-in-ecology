age <- c(0.5, 1, 1.5, 2, 3.3, 4.3, 5.3, 6.3, 7.3, 8.3, 9.3, 10.3, 11.3, 12.3, 14)
length <- c(10.67, 15.81, 26.63, 32.29, 39.85, 47.03, 43.65, 65.15, 49.68, 53.97, 
            52.69, 55.98, 53.51, 61.32, 67.56)
n <- length(age)

Lmax_grid <- seq(50, 100, length.out = 10)
K_grid <- seq(0.01, 0.6, length.out = 10)
t0_grid <- seq(-2, 1, length.out = 10)
sigma_grid <- seq(0.01, 0.5, length.out = 10)

grid_result <- expand.grid(Lmax = Lmax_grid, K = K_grid, t0 = t0_grid, sigma = sigma_grid)
log_posterior <- numeric(nrow(grid_result))
