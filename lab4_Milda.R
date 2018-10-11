

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
plot(x=xGrid,y=posteriorP5$posterior_mean, type="l", ylim = c(-5,-1))
# Add 95% probability bounds (check bayesian course, we must have doen it there)
lines(x=xGrid,y=upp_band5, type="l",col="blue")
lines(x=xGrid,y=low_band5, type="l",col="blue")

### Question2

install.packages('kernlab')
install.packages("AtmRay") # To make 2D grid like in Matlab's meshgrid
library(kernlab)
library(AtmRay)

dataQ2 <- read.csv("/home/milpo192/Documents/AdvML_Lab4/TempTullinge.csv", header=TRUE, sep=";")
time <- 1:2190
day <- 1:365
Stime <- seq(1,2190,by=5)
Sday <- seq(1,365,by=5)
