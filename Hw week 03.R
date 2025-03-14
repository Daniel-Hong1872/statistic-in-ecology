#bootstrap only need 999 times, first one is the original data
#always use the observed data mean, median, b0 and b1, not mean of bootstrap mean
#write together
#compare all together

setwd("C:/Users/dan91/Rstudio/statistic-in-ecology/stat_data")
library(readxl)
f_c_den <- read_excel("mock_dataTS/enviANDdensity.xls")
# 1. Compute the mean and SE(mean) for the fish and copepod densities
#Plot the histogram of bootstrapped means with bootstrap 1000 times.

# normal theory
# mean
Normal_f_mean <- mean(f_c_den$`FishDensity (ind./1000m3)`)
Normal_c_mean <- mean(f_c_den$`CopepodDensity ind./m3)`)
# SE:
Normal_f_SE <- sd(f_c_den$`FishDensity (ind./1000m3)`) / sqrt(length(f_c_den$`FishDensity (ind./1000m3)`))
Normal_c_SE <- sd(f_c_den$`CopepodDensity ind./m3)`) / sqrt(length(f_c_den$`CopepodDensity ind./m3)`))

#non-parametric bootstrap
# mean:
B <- 999
n <- length(f_c_den$`FishDensity (ind./1000m3)`)
bootstrap_fish_mean <- c(Normal_f_mean, numeric(B))
bootstrap_copepod_mean <- c(Normal_c_mean, numeric(B))

set.seed(123)
for(i in 1:B){
  r <- ceiling(runif(n, 1, n))
  resample_f_mean <- f_c_den$`FishDensity (ind./1000m3)`[r]
  resample_c_mean <- f_c_den$`CopepodDensity ind./m3)`[r]
  
  bootstrap_fish_mean[i + 1] <- mean(resample_f_mean)
  bootstrap_copepod_mean[i + 1] <- mean(resample_c_mean)
    }
par(mfrow = c(1, 2))
hist(bootstrap_fish_mean, breaks = 50)
hist(bootstrap_copepod_mean, breaks = 50)

#SE(mean)
bootstrap_fish_SE <- sd(bootstrap_fish_mean)
print(bootstrap_fish_SE)

bootstrap_copepod_SE <- sd(bootstrap_copepod_mean)
print(bootstrap_copepod_SE)

# 2. Compute the median and bootstrapped SE(median) for the fish and copepod densities. 
#Plot the histogram of bootstrapped medians with bootstrap 1000 times.

#fish median
f_median <- median(f_c_den$`FishDensity (ind./1000m3)`)

#copepod median
c_median <- median(f_c_den$`CopepodDensity ind./m3)`)

#bootstrap fish SE(median)
bootstrap_fish_median <- c(f_median, numeric(B))
bootstrap_copepod_median <- c(c_median, numeric(B))

#Bootstrap median
set.seed(123)
for(i in 1:B){
  r <- ceiling(runif(n, 1, n))
  resample_f_median <- f_c_den$`FishDensity (ind./1000m3)`[r]
  resample_c_median <- f_c_den$`CopepodDensity ind./m3)`[r]
  
  bootstrap_fish_median[i + 1] <- median(resample_f_median)
  bootstrap_copepod_median[i + 1] <- median(resample_c_median)
  }

#SE(median)
Bootstrap_f_SE_median <- sd(bootstrap_fish_median)
print(Bootstrap_f_SE_median)
hist(bootstrap_fish_median, breaks = 50)

Bootstrap_c_SE_median <- sd(bootstrap_copepod_median)
print(Bootstrap_c_SE_median)
hist(bootstrap_copepod_median, breaks = 50)

#Plot fish (dependent) v.s copepod (independent) and the regression line. 
#Compute the regression coefficients for fish = b0+b1*copepod and bootstrapped SE(b0) and SE(b1). 
#Plot the histogram of bootstrapped b0 and b1 with bootstrap 1000 times.

x <- cbind(1, f_c_den$`CopepodDensity ind./m3)`)
y <- f_c_den$`FishDensity (ind./1000m3)`
b <- solve(t(x) %*% x) %*% t(x) %*% y
b0 <- b[1]
b1 <- b[2]

y_hat <- x %*% b  #prediction
residuals <- y - y_hat  #residual
sigma_hat <- sqrt(sum(residuals^2) / (n - 2))  # standard deviation of residual

Normal_SE_b1 <- sigma_hat / sqrt(sum((f_c_den$`CopepodDensity ind./m3)` - mean(f_c_den$`CopepodDensity ind./m3)`))^2))
Normal_SE_b0 <- sigma_hat * sqrt(1/n + (mean(f_c_den$`CopepodDensity ind./m3)`)^2) / sum((f_c_den$`CopepodDensity ind./m3)` - mean(f_c_den$`CopepodDensity ind./m3)`))^2))

print(Normal_SE_b0)
print(Normal_SE_b1)

par(mfrow = c(1, 1))
plot(f_c_den$`CopepodDensity ind./m3)`, y, main = "Fish v.s. Copepod", xlab = "copepod density", ylab = "fish density")
abline(b0, b1, col = "lightblue", lwd = 1.5)

bootstrap_b0 <- c(b0, numeric(B))
bootstrap_b1 <- c(b1, numeric(B))

set.seed(123)
for(i in 1:B){
  r <- ceiling(runif(n, 1, n))
  
  resample_x <- x[r,]
  resample_y <- y[r]
  
  bootstrap_b <- solve(t(resample_x) %*% resample_x) %*% t(resample_x) %*% resample_y
  
  bootstrap_b0[i + 1] <- bootstrap_b[1]
  bootstrap_b1[i + 1] <- bootstrap_b[2]
}

bootstrap_SE_b0 <- sd(bootstrap_b0)
bootstrap_SE_b1 <- sd(bootstrap_b1)

par(mfrow = c(1, 2))
hist(bootstrap_b0, breaks = 50)
hist(bootstrap_b1, breaks = 50)

# 1. Compute the mean and stand error of the mean for the fish and copepod density (all data points) respectively using Jackknife. 
#Plot the histogram of Jackknife means.

jackknife_f_mean <- numeric(n)
jackknife_c_mean <- numeric(n)

for(i in 1:n){
  jackknife_f_mean[i] <- mean(f_c_den$`FishDensity (ind./1000m3)`[-i])
  jackknife_c_mean[i] <- mean(f_c_den$`CopepodDensity ind./m3)`[-i])
}
Jackknife_f_mean <- mean(jackknife_f_mean)
Jackknife_c_mean <- mean(jackknife_c_mean)

jackknife_f_SE_mean <- sqrt((n - 1) / n * sum((jackknife_f_mean - mean(jackknife_f_mean)) ^ 2))
jackknife_c_SE_mean <- sqrt((n - 1) / n * sum((jackknife_c_mean - mean(jackknife_c_mean)) ^ 2))

hist(jackknife_f_mean, breaks = 50)
hist(jackknife_c_mean, breaks = 50)

# 2. Compute the regression coefficients for fish = b0+b1*copepod and Jackknife SE of b0 and b1.
#Plot the histogram of Jackknife b0 and b1.

jackknife_b0 <- numeric(n)
jackknife_b1 <- numeric(n)

for(i in 1:n){
  jackknife_x <- cbind(1, f_c_den$`CopepodDensity ind./m3)`[-i])
  jackknife_y <- f_c_den$`FishDensity (ind./1000m3)`[-i]
  
  jackknife_b <- solve(t(jackknife_x) %*% jackknife_x) %*% t(jackknife_x) %*% jackknife_y
  
  jackknife_b0[i] <- jackknife_b[1]
  jackknife_b1[i] <- jackknife_b[2]
}

jackknife_SE_b0 <- sqrt((n - 1) / n * sum((jackknife_b0 - mean(jackknife_b0)) ^ 2))
jackknife_SE_b1 <- sqrt((n - 1) / n * sum((jackknife_b1 - mean(jackknife_b1)) ^ 2))

hist(jackknife_b0, breaks = 50)
hist(jackknife_b1, breaks = 50)

# 3. Compare the estimates for Q1 and Q2 obtained from normal theory, bootstrap, and jackknife.

mean <- c(Normal_f_mean, Normal_c_mean)
median <- c(f_median, c_median)
f_c_compare <- data.frame(mean, median, row.names = c("fish density", "copepod density"))
f_c_compare

b0_b1 <- data.frame(b0, b1)
b0_b1

fish_SE_mean <- c(Normal_f_SE, bootstrap_fish_SE, jackknife_f_SE_mean)
copepod_SE_mean <- c(Normal_c_SE, bootstrap_copepod_SE, jackknife_c_SE_mean)
SE_b0 <- c(Normal_SE_b0, bootstrap_SE_b0, jackknife_SE_b0)
SE_b1 <- c(Normal_SE_b1, bootstrap_SE_b1, jackknife_SE_b1)
method_compare <- data.frame(fish_SE_mean, copepod_SE_mean, SE_b0, SE_b1, row.names = c("Normal", "Bootstrap", "Jackknife"))
method_compare