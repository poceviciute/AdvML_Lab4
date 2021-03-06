#### Advanced machine learning

# Lab 4

####################################### Question 1

#a)

conf_band <- function(pred_mean, pred_var){
  upp_band <- pred_mean+qnorm(0.975)*sqrt(diag(pred_var))
  low_band <- pred_mean-qnorm(0.975)*sqrt(diag(pred_var))
  return(list(upper=upp_band, lower=low_band))
}

# SE_kernel <- function(x1,x2,sigmaF,l){
#   n1 <- length(x1)
#   n2 <- length(x2)
#   K <- matrix(NA,n1,n2)
#   for (i in 1:n2){
#     K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
#   }
#   return(K)
# }

SE_kernel <- function(sigmaf, ell){
  val <- function(x1, x2){
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA,n1,n2)
    for (i in 1:n2){
      K[,i] <- sigmaf^2*exp(-0.5*( (x1-x2[i])/ell)^2 )
    }
    return(K)
  }
  class(val) <- "kernel"
  return(val)
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
  kernel <- SE_kernel(sigmaF, l)
  k <- kernel(x, x)
  L <- t(chol(k+sigmaNoise^2*diag(ncol(k))))
  kStar <- kernel(x, xStar)
  alpha <- solve(t(L))%*%(solve(L)%*%y)
  fStar <- t(kStar)%*%alpha #predictive mean
  v <- solve(L)%*%kStar
  V <- kernel(xStar, xStar)-t(v)%*%v #predictive variance
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

####################################### Question 2

#a)
library(readr)
library(kernlab)
library(AtmRay)

TempTullinge <- read.csv("D:/LiU/732A96/AdvML_Lab4/TempTullinge.csv", header=TRUE, sep=";")

time <- 1:length(TempTullinge$date)
day <- 1:365
fTime <- seq(from=time[1], to=time[length(time)], by = 5)
fDay <- seq(from=day[1], to=day[length(day)], by = 5)
temp <- TempTullinge$temp[fTime]


#What sigmaf and ell??
kernel_func <- SE_kernel(sigmaf=1, ell=0.3)
point1 <- c(1,2)
kernel_point1 <- kernel_func(point1[1], point1[2])
x_vec <- c(1,3,4)
xStar_vec <- c(2,3,4)
k_cov1<-kernelMatrix(kernel = kernel_func, x = x_vec, y = xStar_vec)

#b)
lm_fit <- lm(temp~time+time^2, data = TempTullinge)
sigma2n <- var(lm_fit$residuals)
kernel_func2 <- SE_kernel(sigmaf=20, ell=0.2)
gp_model <- gausspr(fTime,temp,data=TempTullinge, kernel=kernel_func2, var=sigma2n)
post_mean <- predict(gp_model, fTime)
plot(fTime, temp, ylim=c(-30,35))
lines(fTime,post_mean, col="blue", lwd=2)

#c)
x <- scale(fTime)
xs <- scale(fTime)
n <- length(fTime)
Kss <- kernelMatrix(kernel = kernel_func2, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = kernel_func2, x = x, y = x)
Kxs <- kernelMatrix(kernel = kernel_func2, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + sigma2n*diag(n), Kxs) # Covariance matrix of fStar

# Probability intervals for fStar
lines(fTime, post_mean - 1.96*sqrt(diag(Covf)), col = "darkgreen", lwd=2)
lines(fTime, post_mean + 1.96*sqrt(diag(Covf)), col = "darkgreen", lwd=2)

# Prediction intervals for yStar
lines(fTime, post_mean - 1.96*sqrt((diag(Covf) + sigma2n)), col = "purple", lwd=2)
lines(fTime, post_mean + 1.96*sqrt((diag(Covf) + sigma2n)), col = "purple", lwd=2)


#d)
fDay_rep <- rep(fDay,6)
lm_fit2 <- lm(temp~fDay_rep+fDay_rep^2)
sigma2n_day <- var(lm_fit2$residuals)
kernel_func3 <- SE_kernel(sigmaf=20, ell=1.2)
gp_model2 <- gausspr(fDay_rep,temp,data=TempTullinge, kernel=kernel_func3, var=sigma2n_day)
post_mean2 <- predict(gp_model2, fDay_rep)


lines(fTime,post_mean2, col="red", lwd=2)


#e)

periodic_kernel <- function(sigmaf, ell1, ell2, d){
  val <- function(x1, x2){
    n1 <- length(x1)
    n2 <- length(x2)
    K <- matrix(NA,n1,n2)
    for (i in 1:n2){
      K[,i] <- sigmaf^2*exp(-2*sin(pi*abs(x1-x2[i])/d)^2/ell1^2)*exp(-0.5*abs(x1-x2[i])^2/ell2^2)
    }
    return(K)
  }
  class(val) <- "kernel"
  return(val)  
}

periodic_k <- periodic_kernel(20, 1, 10, 365/sd(time))
gp_model3 <- gausspr(fTime,temp,data=TempTullinge, kernel=periodic_k, var=sigma2n_day)
post_mean3 <- predict(gp_model3, fTime)


lines(fTime,post_mean3, col="green", lwd=2)


####################################### Question 3

fraud_data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")
names(fraud_data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud") 
fraud_data[,5] <- as.factor(fraud_data[,5])
set.seed(111)
SelectTraining <- sample(1:dim(fraud_data)[1], size = 1000, replace = FALSE)
fraud_train <- fraud_data[SelectTraining,]
fraud_test <- fraud_data[-SelectTraining,]

#a)

gp_fraud1 <- gausspr(fraud~varWave+skewWave, data=fraud_train)
x1 <- seq(min(fraud_data$varWave),max(fraud_data$varWave),length=100)
x2 <- seq(min(fraud_data$skewWave),max(fraud_data$skewWave),length=100)
fraud_grid <- meshgrid(x1,x2)
grid_points <- data.frame(cbind(c(fraud_grid$x), c(fraud_grid$y)))
#grid_points <- data.frame(rbind("varWave"=fraud_grid$x, "skewWave"=fraud_grid$y))
names(grid_points) <- c("varWave","skewWave")
probPreds1 <- predict(gp_fraud1, grid_points, type="probabilities")

contour(x1,x2,matrix(probPreds1[,1], 100, byrow = TRUE), xlab = "varWave", ylab = "skewWave", main = 'Prob(fraud) - fraud is blue')
points(fraud_train[fraud_train$fraud==0,"varWave"],fraud_train[fraud_train$fraud==0,"skewWave"],col="red",pch=20)
points(fraud_train[fraud_train$fraud==1,"varWave"],fraud_train[fraud_train$fraud==1,"skewWave"],col="blue",pch=20)

conf_mat1<-table(predict(gp_fraud1, fraud_train[,1:2]), fraud_train$fraud)
accuracy1<-sum(diag(conf_mat1))/sum(conf_mat1)

#b)

conf_mat2<-table(predict(gp_fraud1, fraud_test[,1:2]), fraud_test$fraud)
accuracy2<-sum(diag(conf_mat2))/sum(conf_mat2)

#c)

gp_fraud2 <- gausspr(fraud~., data=fraud_train)
conf_mat3<-table(predict(gp_fraud2, fraud_test), fraud_test$fraud)
accuracy3<-sum(diag(conf_mat3))/sum(conf_mat3)
