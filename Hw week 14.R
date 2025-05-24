# ===== 第一題：Uniform Prior =====
library(dplyr)
library(ggplot2)

# 原始資料
age <- c(0.5, 1, 1.5, 2, 3.3, 4.3, 5.3, 6.3, 7.3, 8.3, 9.3, 10.3, 11.3, 12.3, 14)
length <- c(10.67, 15.81, 26.63, 32.29, 39.85, 47.03, 43.65, 65.15, 49.68, 53.97, 
            52.69, 55.98, 53.51, 61.32, 67.56)
n <- length(age)

# 定義參數範圍
Lmax_grid <- seq(50, 100, length.out = 25)
K_grid <- seq(0.01, 0.6, length.out = 25)
t0_grid <- seq(-2, 1, length.out = 25)
sigma_grid <- seq(0.01, 0.5, length.out = 25)

# 建立參數組合網格
grid_result_1 <- expand.grid(Lmax = Lmax_grid, K = K_grid, t0 = t0_grid, sigma = sigma_grid)
log_posterior_1 <- numeric(nrow(grid_result_1))

# 計算對數似然
for(i in seq_len(nrow(grid_result_1))) {
  Lmax <- grid_result_1$Lmax[i]
  K <- grid_result_1$K[i]
  t0 <- grid_result_1$t0[i]
  sigma <- grid_result_1$sigma[i]
  
  L_hat <- Lmax * (1 - exp(-K * (age - t0)))
  L_hat[L_hat <= 0] <- 1e-6  # 避免 log 為負
  
  log_l <- -sum((log(length) - log(L_hat))^2) / (2 * sigma^2) - n * log(sigma)
  log_posterior_1[i] <- log_l
}

# 標準化為 posterior
posterior_1 <- exp(log_posterior_1 - max(log_posterior_1))
posterior_1 <- posterior_1 / sum(posterior_1)
grid_result_1$posterior <- posterior_1

# Marginal posterior summaries
marginal_Lmax_1 <- grid_result_1 %>% 
  group_by(Lmax) %>% 
  summarise(posterior = sum(posterior))
mean_Lmax_1 <- sum(marginal_Lmax_1$Lmax * marginal_Lmax_1$posterior)
median_Lmax_1 <- marginal_Lmax_1$Lmax[which.min(abs(cumsum(marginal_Lmax_1$posterior) - 0.5))]

marginal_K_1 <- grid_result_1 %>% group_by(K) %>% summarise(posterior = sum(posterior))
mean_K_1 <- sum(marginal_K_1$K * marginal_K_1$posterior)
median_K_1 <- marginal_K_1$K[which.min(abs(cumsum(marginal_K_1$posterior) - 0.5))]

marginal_t0_1 <- grid_result_1 %>% group_by(t0) %>% summarise(posterior = sum(posterior))
mean_t0_1 <- sum(marginal_t0_1$t0 * marginal_t0_1$posterior)
median_t0_1 <- marginal_t0_1$t0[which.min(abs(cumsum(marginal_t0_1$posterior) - 0.5))]

marginal_sigma_1 <- grid_result_1 %>% group_by(sigma) %>% summarise(posterior = sum(posterior))
mean_sigma_1 <- sum(marginal_sigma_1$sigma * marginal_sigma_1$posterior)
median_sigma_1 <- marginal_sigma_1$sigma[which.min(abs(cumsum(marginal_sigma_1$posterior) - 0.5))]

# Lmax 圖（Uniform Prior + Posterior）
p_Lmax_1 <- ggplot(marginal_Lmax_1, aes(x = Lmax, y = posterior)) +
  geom_area(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = mean_Lmax_1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = median_Lmax_1, color = "blue", linetype = "dashed") +
  labs(title = "Posterior of L∞ (Uniform Prior)", x = expression(L[infinity]), y = "Density") +
  theme_minimal()

# K 圖
p_K_1 <- ggplot(marginal_K_1, aes(x = K, y = posterior)) +
  geom_area(fill = "forestgreen", alpha = 0.6) +
  geom_vline(xintercept = mean_K_1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = median_K_1, color = "blue", linetype = "dashed") +
  labs(title = "Posterior of K (Uniform Prior)", x = "K", y = "Density") +
  theme_minimal()

# t0 圖
p_t0_1 <- ggplot(marginal_t0_1, aes(x = t0, y = posterior)) +
  geom_area(fill = "orange", alpha = 0.6) +
  geom_vline(xintercept = mean_t0_1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = median_t0_1, color = "blue", linetype = "dashed") +
  labs(title = expression("Posterior of " * t[0]), x = expression(t[0]), y = "Density") +
  theme_minimal()

# sigma 圖
p_sigma_1 <- ggplot(marginal_sigma_1, aes(x = sigma, y = posterior)) +
  geom_area(fill = "purple", alpha = 0.6) +
  geom_vline(xintercept = mean_sigma_1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = median_sigma_1, color = "blue", linetype = "dashed") +
  labs(title = expression("Posterior of " * sigma), x = expression(sigma), y = "Density") +
  theme_minimal()

(p_Lmax_1 | p_K_1) / (p_t0_1 | p_sigma_1)

# ===== 第二題：Correlated Normal Prior =====
library(mvtnorm)
library(patchwork)

# 建立參數網格
grid_result_2 <- expand.grid(Lmax = Lmax_grid, K = K_grid, t0 = t0_grid, sigma = sigma_grid)

# prior 參數
mu <- c(86, 0.13)
Sigma <- matrix(c(10^2, -0.6 * 10 * 0.02, -0.6 * 10 * 0.02, 0.02^2), nrow = 2)

# posterior 計算
log_posterior_2 <- numeric(nrow(grid_result_2))

for(i in seq_len(nrow(grid_result_2))) {
  Lmax <- grid_result_2$Lmax[i]
  K <- grid_result_2$K[i]
  t0 <- grid_result_2$t0[i]
  sigma <- grid_result_2$sigma[i]
  
  L_hat <- Lmax * (1 - exp(-K * (age - t0)))
  L_hat[L_hat <= 0] <- 1e-6
  
  log_l <- -sum((log(length) - log(L_hat))^2) / (2 * sigma^2) - n * log(sigma)
  log_prior_Lmax_K <- dmvnorm(c(Lmax, K), mean = mu, sigma = Sigma, log = TRUE)
  log_prior_t0 <- ifelse(t0 >= -2 & t0 <= 1, 0, -Inf)
  log_prior_sigma <- ifelse(sigma >= 0.01 & sigma <= 0.5, 0, -Inf)
  
  log_posterior_2[i] <- log_l + log_prior_Lmax_K + log_prior_t0 + log_prior_sigma
}

posterior_2 <- exp(log_posterior_2 - max(log_posterior_2))
posterior_2 <- posterior_2 / sum(posterior_2)
grid_result_2$posterior <- posterior_2

# Marginal posteriors
marginal_Lmax_2 <- grid_result_2 %>% 
  group_by(Lmax) %>% 
  summarise(posterior = sum(posterior))
mean_Lmax_2 <- sum(marginal_Lmax_2$Lmax * marginal_Lmax_2$posterior)
median_Lmax_2 <- marginal_Lmax_2$Lmax[which.min(abs(cumsum(marginal_Lmax_2$posterior) - 0.5))]

marginal_K_2 <- grid_result_2 %>% group_by(K) %>% summarise(posterior = sum(posterior))
mean_K_2 <- sum(marginal_K_2$K * marginal_K_2$posterior)
median_K_2 <- marginal_K_2$K[which.min(abs(cumsum(marginal_K_2$posterior) - 0.5))]

marginal_t0_2 <- grid_result_2 %>% group_by(t0) %>% summarise(posterior = sum(posterior))
mean_t0_2 <- sum(marginal_t0_2$t0 * marginal_t0_2$posterior)
median_t0_2 <- marginal_t0_2$t0[which.min(abs(cumsum(marginal_t0_2$posterior) - 0.5))]

marginal_sigma_2 <- grid_result_2 %>% group_by(sigma) %>% summarise(posterior = sum(posterior))
mean_sigma_2 <- sum(marginal_sigma_2$sigma * marginal_sigma_2$posterior)
median_sigma_2 <- marginal_sigma_2$sigma[which.min(abs(cumsum(marginal_sigma_2$posterior) - 0.5))]
# ====== 畫圖 ======

# Lmax Posterior with Prior
p_Lmax_2 <- ggplot(marginal_Lmax_2, aes(x = Lmax, y = posterior)) +
  geom_area(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = mean_Lmax_2, color = "red", linetype = "dashed") +
  geom_vline(xintercept = median_Lmax_2, color = "blue", linetype = "dashed") +
  labs(title = "Posterior of L∞", x = expression(L[infinity]), y = "Density") +
  theme_minimal()

# K Posterior with Prior
p_K_2 <- ggplot(marginal_K_2, aes(x = K, y = posterior)) +
  geom_area(fill = "forestgreen", alpha = 0.6) +
  geom_vline(xintercept = mean_K_2, color = "red", linetype = "dashed") +
  geom_vline(xintercept = median_K_2, color = "blue", linetype = "dashed") +
  labs(title = "Posterior of K", x = "K", y = "Density") +
  theme_minimal()

# t0 Posterior
p_t0_2 <- ggplot(marginal_t0_2, aes(x = t0, y = posterior)) +
  geom_area(fill = "orange", alpha = 0.6) +
  geom_vline(xintercept = mean_t0_2, color = "red", linetype = "dashed") +
  geom_vline(xintercept = median_t0_2, color = "blue", linetype = "dashed") +
  labs(title = expression("Posterior of " * t[0]), x = expression(t[0]), y = "Density") +
  theme_minimal()

# sigma Posterior
p_sigma_2 <- ggplot(marginal_sigma_2, aes(x = sigma, y = posterior)) +
  geom_area(fill = "purple", alpha = 0.6) +
  geom_vline(xintercept = mean_sigma_2, color = "red", linetype = "dashed") +
  geom_vline(xintercept = median_sigma_2, color = "blue", linetype = "dashed") +
  labs(title = expression("Posterior of " * sigma), x = expression(sigma), y = "Density") +
  theme_minimal()

# 排版
(p_Lmax_2 | p_K_2) / (p_t0_2 | p_sigma_2)

joint_grid <- expand.grid(Lmax = Lmax_seq, K = K_seq)
joint_grid$density <- dmvnorm(cbind(joint_grid$Lmax, joint_grid$K), 
                              mean = mu, sigma = Sigma)

p_joint <- ggplot(joint_grid, aes(x = Lmax, y = K, fill = density)) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(title = "Joint Prior of Linf & K", x = expression(L[infinity]), y = "K") +
  theme_minimal()

p_Lmax_prior <- ggplot(marginal_Lmax_prior, aes(x = Lmax, y = density)) +
  geom_area(fill = "gray60") +
  labs(title = "Marginal Prior of L∞", x = expression(L[infinity]), y = "Density") +
  theme_minimal()

p_K_prior <- ggplot(marginal_K_prior, aes(x = K, y = density)) +
  geom_area(fill = "gray60") +
  coord_flip() +
  labs(title = "Marginal Prior of K", x = "K", y = "Density") +
  theme_minimal()

layout <- "
B#
AC
"

p_joint + p_Lmax_prior + p_K_prior + plot_layout(design = layout)
