

## Question 1

### Part1

# Setting up the kernel
SquaredExpKernel <- function(x1,x2,sigmaF=1,l=3){
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA,n1,n2)
    for (i in 1:n2){
        K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
    }
    return(K)
}


posteriorGP <- function(X, y, XStar, hyperParam, sigmaNoise) {
    k_star <- SquaredExpKernel(X, XStar, sigmaF = hyperParam$sigmaF, l = hyperParam$l)
    K <- SquaredExpKernel(X, X, sigmaF = hyperParam$sigmaF, l = hyperParam$l)
    noise <- diag(nrow(K)) * sigmaNoise
    L <- t(chol(K + noise))
    
    # Predictive mean
    alpha <- solve(L) %*% y
    posterior_mean <- t(k_star) %*% alpha
    # Predictive variance
    v <- solve(L) %*% k_star
    kk <- SquaredExpKernel(XStar, XStar, sigmaF = hyperParam$sigmaF, l = hyperParam$l)
    posterior_var <- kk - t(v) %*% v
    
    return(list(posterior_mean = posterior_mean, posterior_var = posterior_var))
    
}

### Part 2

library("mvtnorm")

SimGP <- function(m = 0,K,x,...){
    # Simulates nSim realizations (function) form a Gaussian process with mean m(x) and covariance K(x,x')
    # over a grid of inputs (x)
    n <- length(x)
    if (is.numeric(m)) meanVector <- rep(0,n) else meanVector <- m(x)
    if (is.numeric(K)) covMat <- K else covMat <- K(x,x,...)
    f <- rmvnorm(n, mean = meanVector, sigma = covMat)
    return(f)
}

# Information for updating the prior

X_observations <- c(0.4)
y_observations <- c(0.719)
hyperParam2 <- list(sigmaF = 1, l = 0.3)

# Use the data set and new given point to compute the posterior over the grid [-1,1]
xGrid <- seq(-1,1,length=50)
posteriorP2 <- posteriorGP(X = X_observations, y = y_observations, XStar = xGrid, hyperParam = hyperParam2, sigmaNoise = 0.1)

# find probability bounds
error <- qnorm(0.975)*sqrt(diag(posteriorP2$posterior_var))
upp_band <- as.vector(posteriorP2$posterior_mean) + as.vector(error) 
low_band <- as.vector(posteriorP2$posterior_mean) - as.vector(error) 
# Plot the posterior means over a grid
plot(x=xGrid,y=posteriorP2$posterior_mean, type="l",ylim=c(-2,2))
# Add 95% probability bounds (check bayesian course, we must have doen it there)
lines(x=xGrid,y=upp_band, type="l",col="blue")
lines(x=xGrid,y=low_band, type="l",col="blue")

### Part 3

# Update the information

X_observations[2] <- -0.6
y_observations[2] <- -0.044

posteriorP3 <- posteriorGP(X = X_observations, y = y_observations, XStar = xGrid, hyperParam = hyperParam2, sigmaNoise = 0.1)
# find probability bounds
error3 <- qnorm(0.975)*sqrt(diag(posteriorP3$posterior_var))
upp_band3 <- as.vector(posteriorP3$posterior_mean) + as.vector(error3) 
low_band3 <- as.vector(posteriorP3$posterior_mean) - as.vector(error3) 

# Plot the posterior means over a grid
plot(x=xGrid,y=posteriorP3$posterior_mean, type="l",ylim=c(-2,2))
# Add 95% probability bounds (check bayesian course, we must have doen it there)
lines(x=xGrid,y=upp_band3, type="l",col="blue")
lines(x=xGrid,y=low_band3, type="l",col="blue")

### Part 4
X_observations[1] <- -1
X_observations[3:5] <- c(-0.2,0.4,0.8)
y_observations[3:5] <- c(-0.940, 0.719, -0.664)

posteriorP4 <- posteriorGP(X = X_observations, y = y_observations, XStar = xGrid, hyperParam = hyperParam2, sigmaNoise = 0.1)
# find probability bounds
error4 <- qnorm(0.975)*sqrt(diag(posteriorP4$posterior_var))
upp_band4 <- as.vector(posteriorP4$posterior_mean) + as.vector(error4) 
low_band4 <- as.vector(posteriorP4$posterior_mean) - as.vector(error4) 

# Plot the posterior means over a grid
plot(x=xGrid,y=posteriorP4$posterior_mean, type="l",ylim=c(-2,2))
# Add 95% probability bounds (check bayesian course, we must have doen it there)
lines(x=xGrid,y=upp_band4, type="l",col="blue")
lines(x=xGrid,y=low_band4, type="l",col="blue")

### Part 5

hyperParam5 <- list(sigmaF = 1, l = 1)

posteriorP5 <- posteriorGP(X = X_observations, y = y_observations, XStar = xGrid, hyperParam = hyperParam5, sigmaNoise = 0.1)
# find probability bounds
error5 <- qnorm(0.975)*sqrt(diag(posteriorP5$posterior_var))
upp_band5 <- as.vector(posteriorP5$posterior_mean) + as.vector(error5) 
low_band5 <- as.vector(posteriorP5$posterior_mean) - as.vector(error5) 

# Plot the posterior means over a grid
plot(x=xGrid,y=posteriorP5$posterior_mean, type="l", ylim = c(-2,1))
# Add 95% probability bounds (check bayesian course, we must have doen it there)
lines(x=xGrid,y=upp_band5, type="l",col="blue")
lines(x=xGrid,y=low_band5, type="l",col="blue")

### Question2

### Part 1

#setwd( "C:/Users/Milda/Documents/Documents/Master degree/Adv Machine Learning/AdvML_Lab4")
#install.packages('kernlab')
#install.packages("AtmRay") # To make 2D grid like in Matlab's meshgrid
library(kernlab)
library(AtmRay)

dataQ2 <- read.csv("TempTullinge.csv", header=TRUE, sep=";")
time <- 1:2190
day <- 1:365
Stime <- seq(1,2190,by=5)
Sday <- seq(1,365,by=5)

# Function that returns kernel class object
Matern32 <- function(sigmaf = 1, ell = 1) 
{
  rval <- function(x1,x2){
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA,n1,n2)
    for (i in 1:n2){
      K[,i] <- sigmaf^2*exp(-0.5*( (x1-x2[i])/ell)^2 )
    }
    return(K)
  }
  class(rval) <- "kernel"
  return(rval)
} 

MaternFunc = Matern32(sigmaf = 1, ell = 2)
#evaluate it in the point x = 1; x' = 2
MaternFunc(1,2)

# use the kernelMatrix function to compute the covariance matrix K(X;X*)
#points <- matrix(c(1,3,4,2,3,4), ncol=2)
x1 <- c(1,3,4)
x2 <- c(2,3,4)

covM <- kernelMatrix(MaternFunc, x1,x2)




### Part 2

## With the subsetted data
temp <- dataQ2$temp[Stime]
# Find the noise value
polyFit <- lm(temp ~  Stime + Stime^2)
sigmaNoise = sd(polyFit$residuals)
MaternFunc2 = Matern32(sigmaf = 20, ell = 0.2)

# Fit the GP with built in Square expontial kernel (called rbfdot in kernlab)
GPfit <- gausspr(Stime, temp, kernel = MaternFunc2, var = sigmaNoise^2)
postMean2 <- predict(GPfit,Stime)

plot(Stime, temp, ylim=c(-30,40))
lines(Stime,postMean2, col="red", lwd = 2)


### Part 3

x<-scale(Stime)
xs<-scale(Stime) # XStar
n <- length(x)
SEkernel <- MaternFunc2
Kss <- kernelMatrix(kernel = SEkernel, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = SEkernel, x = x, y = x)
Kxs <- kernelMatrix(kernel = SEkernel, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), Kxs) # Covariance matrix of fStar

# Probability intervals for fStar
lines(Stime, postMean2 - 1.96*sqrt(diag(Covf)), col = "green")
lines(Stime, postMean2 + 1.96*sqrt(diag(Covf)), col = "green")

# Prediction intervals for yStar
lines(Stime, postMean2 - 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")
lines(Stime, postMean2 + 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")


### Part 4
# Find the noise value
polyFit4 <- lm(temp ~  rep.int(Sday,6) )
sigmaNoise4 = sd(polyFit4$residuals)

# Fit the GP with built in Square expontial kernel (called rbfdot in kernlab)
GPfit4 <- gausspr(rep.int(Sday,6), temp, kernel = MaternFunc2, var = sigmaNoise4^2)
postMean4 <- predict(GPfit4,rep.int(Sday,6))


plot(rep.int(Sday,6), temp)
lines(rep.int(Sday,6),postMean4, col="red", lwd = 2) #how and were to include the time variable?...

plot(Stime, temp)
lines(Stime,postMean4, col="red", lwd = 2) #how and were to include the time variable?...


### Part 5

sd <- sqrt(var(time))
d <- 365/sd
# Function that returns kernel class object
Matern5 <- function(sigmaf = 1, ell1 = 1, ell2 = 1) 
{
  rval <- function(x, y = NULL) {
    r = sqrt(crossprod(x-y));
    return(sigmaf^2*exp( -2*(sin(pi*abs(x-y)/d))^2/ell1)*exp(-0.5*(x-y)^2/ell2))
  }
  class(rval) <- "kernel"
  return(rval)
} 

MaternFunc5 = Matern5(sigmaf = 20, ell1 = 1, ell2=10)
GPfit5 <- gausspr(rep.int(Sday,6), temp, kernel = MaternFunc5, var = sigmaNoise4^2)
postMean5 <- predict(GPfit5,rep.int(Sday,6))

plot(rep.int(Sday,6), temp)
lines(rep.int(Sday,6),postMean4, col="red", lwd = 2) #how and were to include the time variable?...
lines(rep.int(Sday,6),postMean5, col="blue", lwd = 2) #how and were to include the time variable?...

## Question 3

# read data from github
data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv",
                 header=FALSE, sep=",")
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])
head(data)
set.seed(111); 
SelectTraining <- sample(1:dim(data)[1], size = 1000,
                                        replace = FALSE)
training <- data[SelectTraining,]
testing <- data[-SelectTraining,]
### Part 1
GPfitFraud <- gausspr(fraud ~ varWave + skewWave, data=training)

# class probabilities 
probPreds <- predict(GPfitFraud, training, type="probabilities")
x1 <- seq(min(training[,1]),max(training[,1]),length=100)
x2 <- seq(min(training[,2]),max(training[,2]),length=100)
gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))

gridPoints <- data.frame(gridPoints)
names(gridPoints) <- names(training)[1:2]
probPreds <- predict(GPfitFraud, gridPoints, type="probabilities")

# Plotting for Prob(setosa)
contour(x1,x2,matrix(probPreds[,1],100,byrow = TRUE), 20, xlab = "Petal.Length", ylab = "Petal.Width", main = 'Prob(Setosa) - Setosa is red')
points(training[training[,5]==0,1],training[training[,5]==0,2],col="blue")
points(training[training[,5]==1,1],training[training[,5]==1,2],col="red")

### Part 2

# predict on the testing set
predT <- predict(GPfitFraud,testing)
table(predT, testing$fraud) # confusion matrix

### Part 3
GPfitFraud3 <- gausspr(fraud ~ ., data=training)

# predict on the testing set
pred3 <- predict(GPfitFraud3,testing)
table(pred3, testing$fraud) # confusion matrix

