# 1a. Generate 10000 random numbers from Gaussian distribution with mean=20 and variance =10, and plot the distribution.

g_dis <- rnorm(10000, 20, 10)

hist(g_dis)

# 1b. Generate 10000 random numbers from Binomial distribution with p=0.5 and n=40, and plot the distribution.
 
n <- 40
p <- 0.5
size <- 10000
 
b_dis <- rep(0, size)

for(i in 1:size){
  b_dis[i] <- sum(runif(n) < p)
}

hist(b_dis)

# Compare the distribution of 1a and 1b, what do you find?

par(mfrow = c(1, 2))
hist(g_dis)
hist(b_dis)

# 2. Make a program that can select our candidates for presentation next week. This program should select randomly but avoid selecting the numbers that had been selected before.

candi <- function (m, w){
  selected <- numeric(0)
  candidate <- matrix(0, nrow = w, ncol = 2)
  
  for(week in 1:w){
    repeat{
      num1 <- round(runif(1, 1, m))
      num2 <- round(runif(1, 1, m))
      
      if(num1 != num2 && !(num1 %in% selected) && !(num2 %in% selected)){
        candidate[week, ] <- c(num1, num2)
        selected <- c(selected, num1, num2)
        break
      }
    }
  }
  return(candidate)
}

candi(34, 16)
