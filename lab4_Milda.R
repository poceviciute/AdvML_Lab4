

### Question 1

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

## Question 2

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

# Generate the data set

XStar2 <- 0.4
yStar2 <- 0.719

# Use the data set and new given point to compute the posterior
hyperParam2 <- list(sigmaF = 1, l = 0.3)
# Which one is the right one??
#posteriorGP(X = xGrid, y = ySim, XStar = XStar2, hyperParam = hyperParam2, sigmaNoise = 0.1)
posteriorQ2 <- posteriorGP(X = XStar2, y = yStar2, XStar = XStar2, hyperParam = hyperParam2, sigmaNoise = 0.1)

## Somehow I need to use the posterior mean, posterior variance in generating values for plotting (potentially f ??)
# Plot the posterior mean of f over the interval [-1,1]???
xGrid <- seq(-1,1,length=20)
ySim2 <- rnorm(length(xGrid),mean=posteriorQ2$posterior_mean, sd = sqrt(posteriorQ2$posterior_var))

# Plot over a grid when I keep updating posterior with all X points in the interval [-1,1]
means <- c()
vars <- c()
temp <- c()
for (i in 1:length(xGrid)){
    temp <- posteriorGP(X = xGrid, y = ySim, XStar = xGrid[1:i], 
                             hyperParam = hyperParam2, sigmaNoise = 0.1)
means[i] <- mean(temp$posterior_mean)
vars[i] <- mean(temp$posterior_var)
}

#############
# The hint suggests that I should increase the number of points used as XStar, but then mt algorithm returns 
# multiple means. I did not find anywhere info how to deal with several XStar points.. I calculated the mean 
# of several means for different time points, but is it ok? 
# Example of what my function returns
posteriorGP(X = xGrid, y = ySim, XStar = xGrid[1:5], 
            hyperParam = hyperParam2, sigmaNoise = 0.1)
## Or maybe we don't need the loop - just have the posterior means calculated for all X observations?

temp <- posteriorGP(X = xGrid, y = ySim, XStar = xGrid, 
                    hyperParam = hyperParam2, sigmaNoise = 0.1)
means <- temp$posterior_mean
vars <- temp$posterior_var

### How to build confidence interval????

Upbound <- means + qnorm(0.95)*vars
Lobound <- means - qnorm(0.95)*vars

plot(x=xGrid,y=means,type="l")
lines(xGrid,Upbound, col = "blue", lwd = 2)
lines(xGrid,Lobound, col = "blue", lwd = 2)




