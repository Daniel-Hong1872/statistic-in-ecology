#Compute the regression coefficients for fish = b0+b1*copepod
#using bootstrap method to calculate the 95% confidence limits of b1 
#test whether b1 is significantly different from 0 (with bootstrap 1000 times). 
#Please calculate the CL using three methods: percentile, BC and BCa methods.

setwd("C:/Users/dan91/Rstudio/statistic-in-ecology/stat_data")
library(readxl)
f_c_den <- read_excel("mock_dataTS/enviANDdensity.xls")
fish_den <- f_c_den$`FishDensity (ind./1000m3)`
copepod_den <- f_c_den$`CopepodDensity ind./m3)`

x <- cbind(1, copepod_den)
y <- fish_den
b <- solve(t(x) %*% x) %*% t(x) %*% y
b0 <- b[1]
b1 <- b[2]

B <- 999
n <- length(fish_den)

bootstrap_b1 <- c(b1, numeric(B))

set.seed(123)
for(i in 1:B){
  r <- ceiling(runif(n, 1, n))
  
  resample_x <- x[r,]
  resample_y <- y[r]
  
  bootstrap_b <- solve(t(resample_x) %*% resample_x) %*% t(resample_x) %*% resample_y
  
  bootstrap_b1[i + 1] <- bootstrap_b[2]
}

#percentile
sort_b1 <- sort(bootstrap_b1)
lower_percentile <- sort_b1[floor(0.025 * (B + 1))]
upper_percentile <- sort_b1[floor(0.975 * (B + 1))]

lower_percentile
upper_percentile

#BC method
z0 <- qnorm(sum(bootstrap_b1 < b1)/(B + 1))
lower_BC <- sort_b1[floor(pnorm(2 * z0 - 1.96) * (B + 1))]
upper_BC <- sort_b1[floor(pnorm(2 * z0 + 1.96) * (B + 1))]

lower_BC
upper_BC

#BCa method
jackknife_b1 <- numeric(n)

set.seed(123)
for(i in 1:n){
  jackknife_x <- cbind(1, copepod_den[-i])
  jackknife_y <- fish_den[-i]
  
  jackknife_b <- solve(t(jackknife_x) %*% jackknife_x) %*% t(jackknife_x) %*% jackknife_y
  
  jackknife_b1[i] <- jackknife_b[2]
}

mean_jackknife_b1 <- mean(jackknife_b1)
numerator <- sum((mean_jackknife_b1 - jackknife_b1) ^ 3)
denominator <- 6 * (sum((mean_jackknife_b1 - jackknife_b1) ^ 2) ^ (3/2))
a <- numerator / denominator

lower_BCa <- sort_b1[floor(pnorm(z0 + (z0 - 1.96) / (1 - a * (z0 - 1.96))) * (B + 1))]
upper_BCa <- sort_b1[floor(pnorm(z0 + (z0 + 1.96) / (1 - a * (z0 + 1.96))) * (B + 1))]

lower_BCa
upper_BCa

#testing b1
p_value <- 2 * min(mean(bootstrap_b1 <= 0), mean(bootstrap_b1 >= 0))
if (p_value < 0.05){
  print("b1 is significant different from 0")
} else{
  print("b1 is not significant different from 0")
}



#Bootstrap test whether significant difference exists between the density of Oncaea Venusta and Canthocalanus pauper.

cop_compo <- read.table("mock_dataTS/Copepod_composition.txt", header = T)
cop_den <- read.table("mock_dataTS/Cop_density.txt")
cop_sp <- read.table("mock_dataTS/copepodSPlist.txt", fill = T, sep = "\n")

row.names(cop_compo) <- cop_sp$V1
Onve <- cop_compo["Oncaea venusta",] / 100
Capa <- cop_compo["Canthocalanus pauper",] / 100

Onve_den <- matrix(0, nrow = 1, ncol = length(cop_den$V1))
row.names(Onve_den) <- "Oncaea venusta"
for(i in 1:length(cop_den$V1)){
  Onve_den[, i] <- Onve[, i] * cop_den[i,]
}

Capa_den <- matrix(0, nrow = 1, ncol = length(cop_den$V1))
row.names(Capa_den) <- "Canthocalanus pauper"
for(i in 1:length(cop_den$V1)){
  Capa_den[, i] <- Capa[, i] * cop_den[i,]
}

Onve_den_mean <- mean(Onve_den)
Capa_den_mean <- mean(Capa_den)
original_diff <- Onve_den_mean - Capa_den_mean

Bootstrap_diff <- c(original_diff, numeric(B))

set.seed(111)
for(i in 1:B){
  resample_onve <- Onve_den[floor(runif(n, 1, n))]
  resample_capa <- Capa_den[floor(runif(n, 1, n))]
  
  Bootstrap_diff[i + 1] <- mean(resample_onve) - mean(resample_capa)
}

sort_diff <- sort(Bootstrap_diff)
lower_diff <- sort_diff[floor(0.025 * (1 + B))]
upper_diff <- sort_diff[floor(0.975 * (1 + B))]
lower_diff
upper_diff

if(lower_diff > 0 | upper_diff < 0) {
  print("the density of these two species have significant difference")
} else {
  print("the density of these two species have no significant difference")
}

#BCa

Z0 <- qnorm(sum(Bootstrap_diff < original_diff) / (B + 1))
jackknife_diff <- numeric(n)

for(i in 1:n){
  jackknife_onve <- Onve_den[, -i]
  jackknife_capa <- Capa_den[, -i]
  jackknife_diff[i] <- mean(jackknife_onve) - mean(jackknife_capa)
}

mean_jackknife_diff <- mean(jackknife_diff)
numerator_diff <- sum((mean_jackknife_diff - jackknife_diff) ^ 3)
denominator_diff <- 6 * (sum((mean_jackknife_diff - jackknife_diff) ^ 2) ^ (3/2))
A <- numerator_diff / denominator_diff

lower_BCa_diff <- sort(Bootstrap_diff)[floor(pnorm(Z0 + (Z0 - 1.96) / (1 - a * (Z0 - 1.96))) * (B + 1))]
upper_BCa_diff <- sort(Bootstrap_diff)[floor(pnorm(Z0 + (Z0 + 1.96) / (1 - a * (Z0 + 1.96))) * (B + 1))]

if(lower_BCa_diff > 0 | upper_BCa_diff < 0) {
  print("the density of these two species have significant difference")
} else {
  print("the density of these two species have no significant difference")
}
