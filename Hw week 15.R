# 1.
age <- c(0.5, 1, 1.5, 2, 3.3, 4.3, 5.3, 6.3, 7.3, 8.3, 9.3, 10.3, 11.3, 12.3, 14)
length <- c(10.67, 15.81, 26.63, 32.29, 39.85, 47.03, 43.65, 65.15, 49.68, 53.97, 
            52.69, 55.98, 53.51, 61.32, 67.56)

log_likelihood <- function(params){
  Linf <- params[1]
  K <- params[2]
  t0 <- params[3]
  sigma <- params[4]
  
  if (Linf <= 0 || K <= 0 || sigma <= 0) return(-Inf)
  
  L_hat <- Linf * (1 - exp(-K * (age - t0)))
  L_hat[L_hat <= 0] <- 1e-6
  
  ll <- -sum((log(length) - log(L_hat)) ^ 2) / (2 * sigma ^ 2) - length(length) * log(sigma)
  return(ll)
}

#prior function with uniform distribution
log_prior <- function(params){
  Linf <- params[1]
  K <- params[2]
  t0 <- params[3]
  sigma <- params[4]
  
  if(Linf < 50 || Linf > 100) return(-Inf)
  if (K < 0.01 || K > 0.6) return(-Inf)
  if (t0 < -2 || t0 > 1) return(-Inf)
  if (sigma < 0.01 || sigma > 0.5) return(-Inf)
  
  return(0)
}

#MCMC
set.seed(123)

n_iter <- 300000
params <- matrix(NA, nrow = n_iter, ncol = 4)
params[1, ] <- c(59.32, 0.34, -0.06, 0.10) #start by MLE value

for(i in 2:n_iter){
  current <- params[i - 1, ]
  
  proposal <- current + rnorm(4, mean = 0, sd = c(2, 0.03, 0.05, 0.01))
  
  log_post_current <- log_likelihood(current) + log_prior(current)
  log_post_proposal <- log_likelihood(proposal) + log_prior(proposal)
  
  r <- exp(log_post_proposal - log_post_current)
  
  if(runif(1) < r){
    params[i, ] <- proposal
  }else{
   params[i, ] <- current 
  }
}

burnin <- n_iter * 0.2
thinning <- 10
posterior <- as.data.frame(params[seq(burnin + 1, n_iter, by = thinning), ])
colnames(posterior) <- c("Linf", "K", "t0", "sigma")

#trace plot and posterior density
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))

for(i in 1:4){
  name <- colnames(posterior)[i]
  vec <- posterior[[i]]
  
  plot(vec, type = "l", col = i, main = paste("Trace:", name),
       xlab = "Iteration", ylab = name)
  
  dens <- density(vec)
  plot(dens, main = paste("Posterior:", name), col = i, lwd = 2,
       xlab = name)
  abline(v = mean(vec), col = "red", lty = 2)
  abline(v = median(vec), col = "black", lty = 2)
  legend("topright", legend = c(paste("mean =", round(mean(vec), 3)),
                                paste("median =", round(median(vec), 3))),
         text.col = c("red", "black"), bty = "n")
}

par(mfrow = c(1, 1))

ages <- seq(0, 30, by = 0.1)  # 要預測的年齡範圍

# 用 posterior 抽樣生成成長曲線
n_samples <- 1000
set.seed(123)
sample_indices <- sample(1:nrow(posterior), n_samples)

# 建立空矩陣儲存所有抽樣後的 L(t)
L_mat <- matrix(NA, nrow = n_samples, ncol = length(ages))

for (i in 1:n_samples) {
  p <- posterior[sample_indices[i], ]
  Linf <- p$Linf; K <- p$K; t0 <- p$t0; sigma <- p$sigma
  L_hat <- Linf * (1 - exp(-K * (ages - t0)))
  L_mat[i, ] <- L_hat
}

# 每個年齡點的 2.5%、50%、97.5% quantile → 成長可信區間
CI_lower <- apply(L_mat, 2, quantile, probs = 0.025)
CI_median <- apply(L_mat, 2, quantile, probs = 0.5)
CI_upper <- apply(L_mat, 2, quantile, probs = 0.975)

# 預測區間（含 measurement error）
# 每一個 L_hat 加上 lognormal 誤差的標準差
pred_mat <- L_mat * exp(rnorm(length(L_mat), mean = 0, sd = posterior[sample_indices, "sigma"][i]))

Pred_lower <- apply(pred_mat, 2, quantile, probs = 0.025)
Pred_upper <- apply(pred_mat, 2, quantile, probs = 0.975)

# 畫圖
plot(age, length, pch = 16, ylim = range(Pred_lower, Pred_upper), 
     xlim = c(0, 30), xlab = "Age (years)", ylab = "Length (cm)",
     main = "Observed Length-at-Age with 95% CI and 95% Predictive Interval")

# 加入預測區間
polygon(c(ages, rev(ages)), c(Pred_upper, rev(Pred_lower)),
        col = adjustcolor("mediumpurple2", alpha.f = 0.3), border = NA)

# 加入可信區間
polygon(c(ages, rev(ages)), c(CI_upper, rev(CI_lower)),
        col = adjustcolor("lightblue", alpha.f = 0.5), border = NA)

# 中位數成長線
lines(ages, CI_median, col = "blue", lwd = 2)

# 圖例
legend("bottomright", legend = c("Observed", "Median Growth", 
                                 "95% Growth CI", "95% Predictive CI"),
       col = c("black", "blue", "lightblue", "mediumpurple2"),
       pch = c(16, NA, 15, 15), lty = c(NA, 1, NA, NA),
       pt.cex = 2, bty = "n", cex = 0.9)

# 2.
library(MASS)  # for mvrnorm
library(mvtnorm)
library(ggplot2)
library(ggExtra)

# log-likelihood
log_likelihood_2 <- function(params){
  Linf <- params[1]; K <- params[2]; t0 <- params[3]; sigma <- params[4]
  if (any(is.na(c(Linf, K, t0, sigma)))) return(-Inf)
  if (Linf <= 0 || K <= 0 || sigma <= 0) return(-Inf)
  L_hat <- Linf * (1 - exp(-K * (age - t0)))
  L_hat[L_hat <= 0] <- 1e-6
  -sum((log(length) - log(L_hat))^2) / (2 * sigma^2) - length(length) * log(sigma)
}

# prior
mu_2 <- c(86, 0.13)
sigma_L_2 <- 10
sigma_K_2 <- 0.02
rho_2 <- -0.6
cov_matrix_2 <- matrix(c(
  sigma_L_2^2, rho_2 * sigma_L_2 * sigma_K_2,
  rho_2 * sigma_L_2 * sigma_K_2, sigma_K_2^2
), nrow = 2)

log_prior_2 <- function(params){
  Linf <- params[1]; K <- params[2]; t0 <- params[3]; sigma <- params[4]
  if (any(is.na(c(Linf, K, t0, sigma)))) return(-Inf)
  if (t0 < -2 || t0 > 1 || sigma < 0.01 || sigma > 0.5) return(-Inf)
  dmvnorm(c(Linf, K), mean = mu_2, sigma = cov_matrix_2, log = TRUE)
}

par(mfrow = c(1, 1))

#marginal prior of Linf
ggplot(data.frame(Linf = rnorm(10000, mean = 86, sd = 10)), aes(x = Linf)) +
  geom_density(fill = "gray", color = NA, adjust = 2) +
  xlim(60, 110) +
  labs(title = "Marginal Prior of Linf", x = expression(L[infinity]), y = "") +
  theme_minimal(base_size = 14)

#Marginal prior of K
ggplot(data.frame(K = rnorm(10000, mean = 0.13, sd = 0.02)), aes(x = K)) +
  geom_density(fill = "gray", color = NA, adjust = 2) +
  xlim(0.05, 0.2) +
  labs(title = "Marginal Prior of K", x = "K", y = "") +
  theme_minimal(base_size = 14)

#joint prior distribution of Linf and K

Linf_seq <- seq(60, 110, length.out = 200)
K_seq <- seq(0.05, 0.2, length.out = 200)
grid <- expand.grid(Linf = Linf_seq, K = K_seq)

grid$density <- dmvnorm(grid, mean = mu_2, sigma = cov_matrix_2)
# 生成 sample
set.seed(123)
prior_sample <- mvrnorm(100000, mu_2, cov_matrix_2)
prior_df <- as.data.frame(prior_sample)
colnames(prior_df) <- c("Linf", "K")

# 畫圖
ggplot(grid, aes(x = Linf, y = K, fill = density)) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(title = "Joint Prior Distribution of Linf and K",
       x = expression(L[infinity]), y = "K", fill = "Density") +
  theme_minimal(base_size = 14)


# MCMC
set.seed(111)
n_iter_2 <- 300000
params_2 <- matrix(NA, n_iter_2, 4)
params_2[1, ] <- c(86, 0.13, -0.5, 0.1)

for(i in 2:n_iter_2){
  current <- params_2[i - 1, ]
  proposal <- current + rnorm(4, mean = 0, sd = c(2, 0.01, 0.1, 0.01))
  log_post_current <- log_likelihood_2(current) + log_prior_2(current)
  log_post_proposal <- log_likelihood_2(proposal) + log_prior_2(proposal)
  r <- exp(log_post_proposal - log_post_current)
  
  if(runif(1) < r){
    params_2[i, ] <- proposal
  }else{
    params_2[i, ] <- current 
  }
}

burnin_2 <- 0.2 * n_iter_2
thinning_2 <- 10
posterior_2 <- as.data.frame(params_2[seq(burnin_2 + 1, n_iter_2, by = thinning_2), ])
colnames(posterior_2) <- c("Linf", "K", "t0", "sigma")

# 畫 trace + posterior density（不使用 coda）
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
for(i in 1:4){
  name <- colnames(posterior_2)[i]
  vec <- posterior_2[[i]]
  plot(vec, type = "l", col = i, main = paste("Trace:", name), xlab = "Iteration", ylab = name)
  dens <- density(vec)
  plot(dens, main = paste("Posterior:", name), col = i, lwd = 2, xlab = name)
  
  if (name == "Linf") {
    curve(dnorm(x, mean = 86, sd = 10), col = "gray", lty = 3, lwd = 2, add = TRUE)
    legend("topleft", legend = "prior", col = "gray", lty = 3, lwd = 2, bty = "n")
  } else if (name == "K") {
    curve(dnorm(x, mean = 0.13, sd = 0.02), col = "gray", lty = 3, lwd = 2, add = TRUE)
    legend("topleft", legend = "prior", col = "gray", lty = 3, lwd = 2, bty = "n")
  }
  abline(v = mean(vec), col = "red", lty = 2)
  abline(v = median(vec), col = "black", lty = 2)
  legend("topright", legend = c(paste("mean =", round(mean(vec), 3)), paste("median =", round(median(vec), 3))),
         text.col = c("red", "black"), bty = "n")
}

# Growth Curve Plot
ages_2 <- seq(0, 30, by = 0.1)
n_samples_2 <- 1000
sample_indices_2 <- sample(1:nrow(posterior_2), n_samples_2)
L_mat_2 <- matrix(NA, nrow = n_samples_2, ncol = length(ages_2))

for (i in 1:n_samples_2) {
  p <- posterior_2[sample_indices_2[i], ]
  L_mat_2[i, ] <- p$Linf * (1 - exp(-p$K * (ages_2 - p$t0)))
}
CI_lower_2 <- apply(L_mat_2, 2, quantile, probs = 0.025)
CI_median_2 <- apply(L_mat_2, 2, quantile, probs = 0.5)
CI_upper_2 <- apply(L_mat_2, 2, quantile, probs = 0.975)

# 預測區間
pred_mat_2 <- L_mat_2
for (i in 1:n_samples_2) {
  sigma_i <- posterior_2$"sigma"[sample_indices_2[i]]
  pred_mat_2[i, ] <- L_mat_2[i, ] * exp(rnorm(length(ages_2), 0, sigma_i))
}
PI_lower_2 <- apply(pred_mat_2, 2, quantile, probs = 0.025)
PI_upper_2 <- apply(pred_mat_2, 2, quantile, probs = 0.975)

# 畫圖
par(mfrow = c(1, 1))
plot(age, length, pch = 16, ylim = range(PI_lower_2, PI_upper_2),
     xlim = c(0, 30), xlab = "Age (years)", ylab = "Length (cm)",
     main = "Observed Length-at-Age with 95% CI and 95% Predictive Interval")
polygon(c(ages_2, rev(ages_2)), c(PI_upper_2, rev(PI_lower_2)),
        col = adjustcolor("mediumpurple2", alpha.f = 0.3), border = NA)
polygon(c(ages_2, rev(ages_2)), c(CI_upper_2, rev(CI_lower_2)),
        col = adjustcolor("lightblue", alpha.f = 0.5), border = NA)
lines(ages_2, CI_median_2, col = "blue", lwd = 2)
legend("bottomright", legend = c("Observed", "Median Growth", "95% Growth CI", "95% Predictive CI"),
       col = c("black", "blue", "lightblue", "mediumpurple2"),
       pch = c(16, NA, 15, 15), lty = c(NA, 1, NA, NA), pt.cex = 2, bty = "n")
