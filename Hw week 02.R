# 1a. Generate 10000 random numbers from Gaussian distribution with mean=20 and variance =10, and plot the distribution.
#variance=10, sd=square root 10!!!
g_dis1 <- rnorm(10000)
g_dis <- g_dis1  * 10 + 20
mean(g_dis)
sd(g_dis)
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



x <- rnorm(10000) + 20
mean(x)
y <- x * 10
sd(y)


#generate random sequence every time and choose the first two
#remove the first two

candi2 <- function (m, w){
  selected <- numeric(0)
  
  for(week in 1:w){
    repeat{
      num1 <- round(runif(1, 1, m))
      num2 <- round(runif(1, 1, m))
      
      if(num1 != num2 && !(num1 %in% selected) && !(num2 %in% selected)){
        print(c(num1, num2))
        selected <- c(selected, num1, num2)
        break
      }
    }
  }
}
candi2(39, 15)

candi3 <- function(m){
  if (!exists("used", envir = .GlobalEnv)) {
    assign("used", numeric(0), envir = .GlobalEnv)
  }  
  available <- setdiff(1:m, used)
  
  x <- runif(length(available))
  candidate <- available[order(x)][1:2]
  
  assign("used", c(used, candidate), envir = .GlobalEnv)
  
  return(candidate)
}
set.seed(123)
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))
print(candi3(30))

