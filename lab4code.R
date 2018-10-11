#### Advanced machine learning

# Lab 4

# Question 1

#a)

conf_band <- function(pred_mean, pred_var){
  upp_band <- pred_mean+qnorm(0.975)*sqrt(diag(pred_var))
  low_band <- pred_mean-qnorm(0.975)*sqrt(diag(pred_var))
  return(list(upper=upp_band, lower=low_band))
}

SE_kernel <- function(x1,x2,sigmaF,l){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}

posteriorGP <- function(x, y, xStar, hyperParam, sigmaNoise){
  #### Input:
  #x = train data
  #y = targets
  #xStar = test data
  #hyperParam = list(sigmaF, l)
  #sigmaNoise
  #K covariance for indata, call se_kernel
  y <- as.matrix(y)
  sigmaF <- hyperParam$sigmaF
  l <- hyperParam$l
  k <- SE_kernel(x, x, sigmaF, l)
  L <- t(chol(k+sigmaNoise^2*diag(ncol(k))))
  kStar <- SE_kernel(x, xStar, sigmaF, l)
  alpha <- solve(t(L))%*%(solve(L)%*%y)
  fStar <- t(kStar)%*%alpha #predictive mean
  v <- solve(L)%*%kStar
  V <- SE_kernel(xStar, xStar, sigmaF, l)-t(v)%*%v #predictive variance
  #log_marg <- -0.5*t(y)%*%alpha-sum(L) - 0.5*length(x)*log(2*pi)
  #### Output:
  #posterior mean, posterior variance
  return(list("pred_mean"=fStar, "pred_var"=V))
}


#b)
noise <- 0.1
hyperParam <- list("sigmaF"=1, "l"=0.3)
xGrid <- seq(-1,1,length=50)

obs1 <- list(x=0.4, y=0.719)
post1<-posteriorGP(obs1$x, obs1$y, xGrid, hyperParam, noise)
confband1 <- conf_band(post1$pred_mean, post1$pred_var)

plot(xGrid, post1$pred_mean, type="l", ylim=c(-2,2), ylab="Posterior mean", xlab = "x", main="b) One observation")
lines(xGrid,confband1$upper, col="blue")
lines(xGrid,confband1$lower, col="blue")


#c)
obs2 <- list(x=c(obs1$x,-0.6), y=c(obs1$y,-0.044))
post2 <- posteriorGP(obs2$x, obs2$y, xGrid, hyperParam, noise)
confband2 <- conf_band(post2$pred_mean, post2$pred_var)

plot(xGrid, post2$pred_mean, type="l", ylim=c(-2,2), ylab="Posterior mean", xlab = "x", main="c) Two observations")
lines(xGrid,confband2$upper, col="blue")
lines(xGrid,confband2$lower, col="blue")

#d)
obs3 <- list(x=c(obs2$x, -1, -0.2, 0.8), y=c(obs2$y, 0.768, -0.94, -0.664))
post3 <- posteriorGP(obs3$x, obs3$y, xGrid, hyperParam, noise)
confband3 <- conf_band(post3$pred_mean, post3$pred_var)

plot(xGrid, post3$pred_mean, type="l", ylim=c(-2,2), ylab="Posterior mean", xlab = "x", main="d) Five observations")
lines(xGrid,confband3$upper, col="blue")
lines(xGrid,confband3$lower, col="blue")

#e)
hyperParam2 <- list(sigmaF=1, l=1)
post4 <- posteriorGP(obs3$x, obs3$y, xGrid, hyperParam2, noise)
confband4 <- conf_band(post4$pred_mean, post4$pred_var)

plot(xGrid, post4$pred_mean, type="l", ylim=c(-2,2), ylab="Posterior mean", xlab = "x", main="e) sigmaF=1, l=1")
lines(xGrid,confband4$upper, col="blue")
lines(xGrid,confband4$lower, col="blue")

# Question 2

#a)

#b)

#c)

#d)

#e)


# Question 3

#a)

#b)

#c)

