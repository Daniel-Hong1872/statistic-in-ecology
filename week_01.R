tag <- c(2, 3, 5, 7, 8, 9, 15, 21, 23, 26)
weight <- c(14.8, 21, 19.7, 23.2, 16, 16.1, 20, 29.3, 17.8, 21.2)
condition <- c("good", "fair", "fair", "poor", "fair", "good", "good", "fair", "fair", "poor")
fishData <- data.frame(tag, weight, condition) 
for(i in 1:length(fishData[,1])){
  if(fishData[i,2] >= 20 & fishData[i,3] !="poor"){
    fishData$survival[i] = 1
  }
  else{
    fishData$survival[i] = 0
  }
}
fishData

add2 <- function(x, y){
  output <- x + y
  return(output)
}
add2(3,9)


patients <- data.frame( 
  id = c(31, 62, 50, 99, 53, 75, 54, 58, 4, 74), 
  age = c(12, 18, 20, 17, 14, 8, 12, 24, 24, 21), 
  sex = c("M", "F", "F", "M", "F", "M", "M", "F", "F", "M") )
head(patients, n = 2)

#Use a logical operator to display ages that are larger than 20
patients[patients[,2] > 20, 2]

#Do the same as above but also display the corresponding id and sex
patients[patients[,2] > 20,]

#Display only female observations
patients[patients[,3] == "F",]

#Change the age of the 7th patient from 12 to 21
patients[7, 2] <- 21

#Calculate the proportion of subjects that are age 20 or greater
length(patients[patients[,2] >= 20, 2])/length(patients[,2])

#Calculate the proportion of males that are greater than 20
length(patients[patients$age > 20 & patients$sex == "M",])/length(patients[,2])

#Permanently delete the 10th subject 
patients[-10,]

#Permanently add two more subjects to this data frame (use rbind)
rbind(patients)

humidity <- c(63.33, NA, 64.63, 68.38, NA, 79.1, 77.46)
na.omit(humidity)
na.pass(humidity)

setwd("C:/Users/dan91/Rstudio/stat_data//mock_dataTS")
library(readxl)
read_excel("copepod_datasheet.xls")

x = seq(-pi, pi, by = 0.1)
y = sin(x)
#Line plot
plot(x, y, type = "l", col = 5, lwd = 9)
lines(x, cos(x), col = 1, lwd = 2, lty = 2)
#Adding label to the plot
plot(y~x, type = "l", main = "plot", xlab = "time", ylab = "abundance")
#Adjusting the axis
plot(x, y, xlim = c(-5, 5), ylim = c(-1, 1))
#Plot a histogram
hist(rnorm(30))


