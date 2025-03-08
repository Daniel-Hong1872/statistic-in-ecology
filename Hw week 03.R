setwd("C:/Users/dan91/Rstudio/statistic-in-ecology/stat_data")
library(readxl)
f_c_den <- read_excel("mock_dataTS/enviANDdensity.xls")
# 1. Compute the mean and SE(mean) for the fish and copepod densities
#Plot the histogram of bootstrapped means with bootstrap 1000 times.

# normal theory
# mean
mean(f_c_den$`FishDensity (ind./1000m3)`)
mean(f_c_den$`CopepodDensity ind./m3)`)
# SE:
sd(f_c_den$`FishDensity (ind./1000m3)`) / sqrt(length(f_c_den$`FishDensity (ind./1000m3)`))
sd(f_c_den$`CopepodDensity ind./m3)`) / sqrt(length(f_c_den$`CopepodDensity ind./m3)`))

#non-parametric bootstrap
set.seed(123)
# fish mean:
b <- 1000
n <- length(f_c_den$`FishDensity (ind./1000m3)`)
bootstrap_fish_mean <- numeric(b)

for(i in 1:b){
  x <- ceiling(runif(n, 1, n))
  resample_f_mean <- f_c_den$`FishDensity (ind./1000m3)`[x]
  
  bootstrap_fish_mean[i] <- mean(resample_f_mean)
  }

hist(bootstrap_fish_mean, breaks = 50)

#fish SE(mean)
bootstrap_fish_SE <- sd(bootstrap_fish_mean)
print(bootstrap_fish_SE)

#copepod mean
bootstrap_copepod_mean <- numeric(b)

for(i in 1:b){
  x <- ceiling(runif(n, 1, n))
  resample_c_mean <- f_c_den$`CopepodDensity ind./m3)`[x]
  
  bootstrap_copepod_mean[i] <- mean(resample_c_mean)
}

hist(bootstrap_copepod_mean, breaks = 50)

#copepod SE(mean)
bootstrap_copepod_SE <- sd(bootstrap_copepod_mean)
print(bootstrap_copepod_SE)

# 2. Compute the median and bootstrapped SE(median) for the fish and copepod densities. 
#Plot the histogram of bootstrapped medians with bootstrap 1000 times.

#fish median
median(f_c_den$`FishDensity (ind./1000m3)`)

#copepod median
median(f_c_den$`CopepodDensity ind./m3)`)

#bootstrap fish SE(median)
bootstrap_fish_median <- numeric(b)

for(i in 1:b){
  x <- ceiling(runif(n, 1, n))
  resample_f_median <- f_c_den$`FishDensity (ind./1000m3)`[x]
  
  bootstrap_fish_median[i] <- mean(resample_f_median)
}

hist(bootstrap_fish_median, breaks = 50)

#bootstrap copepod SE(median)
bootstrap_copepod_median <- numeric(b)

for(i in 1:b){
  x <- ceiling(runif(n, 1, n))
  resample_c_median <- f_c_den$`CopepodDensity ind./m3)`[x]
  
  bootstrap_copepod_median[i] <- mean(resample_c_median)
}

hist(bootstrap_copepod_median, breaks = 50)

