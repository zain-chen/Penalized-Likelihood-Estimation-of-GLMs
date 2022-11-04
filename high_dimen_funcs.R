setwd("C:/Users/cq/OneDrive/文档")

# SCAD惩罚函数
SCAD <- function(lambda, theta){
  if(abs(theta) > lambda){
    if(3.7*lambda > abs(theta)){
      return((3.7*lambda - abs(theta))/2.7)
    }else{
      return(1)
    }
  }else{
    return(lambda)
  }
}


sigma.lqa <- function(beta, lambda){
  beta.copy = beta
  for(j in 1:p){    
    beta.copy[j] <- SCAD(lambda, beta.copy[j]) / abs(beta.copy[j]) 
  }
  sigma <- diag(as.vector(beta.copy))
  return(sigma)
}

sigma.mm <- function(beta, lambda){
  beta.copy = beta
  for(j in 1:p){    
    beta.copy[j] <- SCAD(lambda, beta.copy[j]) / (abs(beta.copy[j]) + tau)
  }
  sigma <- diag(as.vector(beta.copy))
  return(sigma)
}

Logi.LQA <- function(beta0, x, y, lambda, thd1, thd2){
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  beta <- beta0
  for(ii in 1:50){            # 为了不运行太久最多迭代50次
    ita <- x%*%beta
    sigma <- sigma.lqa(beta, lambda)
    mu <- exp(ita) / (exp(ita) + rep(1,n))
    V <- diag(as.vector(mu * (rep(1,n)-mu)))
    Z <- V %*% ita - (y - mu)
    # 为防不可逆停止程序用广义逆，同时求逆计算如果还是错误就跳出
    ggg <- try(ginv(t(x)%*%V%*%x - n*sigma) %*% t(x) %*% Z)
    if("try-error" %in% class(ggg)){
      break
    }else{
      newbeta <- ggg
    }
    if(NaN %in% newbeta) break      # 如果计算出错则跳出
    if(t(newbeta-beta)%*%(newbeta-beta)<thd2) break
    beta <- newbeta
  }
  for(q in 1:p){
    if(beta[q] < thd1) beta[q] <- 0
  }
  return(beta) 
}

Pois.LQA <- function(beta0, x, y, lambda, thd1, thd2){
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  beta <- beta0
  for(ii in 1:50){
    ita <- x%*%beta
    sigma <- sigma.lqa(beta, lambda)
    mu <- exp(ita) 
    V <- diag(as.vector(mu * mu))
    Z <- V %*% ita - (y - mu)
    ggg <- try(ginv(t(x)%*%V%*%x - n*sigma) %*% t(x) %*% Z)
    if("try-error" %in% class(ggg)){
      break
    }else{
      newbeta <- ggg
    }
    if(NaN %in% newbeta) break
    if(t(newbeta-beta)%*%(newbeta-beta)<thd2) break
    beta <- newbeta
  }
  for(q in 1:p){
    if(beta[q] < thd1) beta[q] <- 0
  }
  return(beta) 
}

Logi.MM <- function(beta0, x, y, lambda, thd1, thd2){
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  beta <- beta0
  for(ii in 1:50){
    ita <- x%*%beta
    sigma <- sigma.mm(beta, lambda)
    mu <- exp(ita) / (exp(ita) + rep(1,n))
    V <- diag(as.vector(mu * (rep(1,n)-mu)))
    Z <- V %*% ita - (y - mu)
    ggg <- try(ginv(t(x)%*%V%*%x - n*sigma) %*% t(x) %*% Z)
    if("try-error" %in% class(ggg)){
      break
    }else{
      newbeta <- ggg
    }
    if(NaN %in% newbeta) break
    if(t(newbeta-beta)%*%(newbeta-beta)<thd2) break
    beta <- newbeta
    for(q in 1:p){
      if(beta[q] < thd1) beta[q] <- 0
    }
  }
  return(beta)
}

Pois.MM <- function(beta0, x, y, lambda, thd1, thd2){
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  beta <- beta0
  for(ii in 1:50){
    ita <- x%*%beta
    sigma <- sigma.mm(beta, lambda)
    mu <- exp(ita) 
    V <- diag(as.vector(mu * mu))
    Z <- V %*% ita - (y - mu)
    ggg <- try(ginv(t(x)%*%V%*%x - n*sigma) %*% t(x) %*% Z)
    if("try-error" %in% class(ggg)){
      break
    }else{
      newbeta <- ggg
    }
    if(NaN %in% newbeta) break
    if(t(newbeta-beta)%*%(newbeta-beta)<thd2) break
    beta <- newbeta
    for(q in 1:p){
      if(beta[q] < thd1) beta[q] <- 0
    }
  }     
  return(beta) 
}


simu.logi <- function(sig, real.beta)
{
  # 按协方差矩阵的改变来模拟不同情况
  # 其他参数指定为全局变量
  correct.lqa <- c()
  incorrect.lqa <- c()
  correct.mm <- c()
  incorrect.mm <- c()
  correct.sis <- c()
  incorrect.sis <- c()
  correct.glmnet <- c()
  incorrect.glmnet <- c()
  ee.mm <- c()
  ee.lqa <- c()
  ee.sis <- c()
  ee.glmnet <- c()
  
  for(i in 1:200){
    # 更新
    x <- rmvnorm(n, rep(0,p), sig)
    colnames(x) <- c("x1","x2","x3","x4","x5","x6","x7","x8")
    u.logi <- exp(x %*% real.beta)/ (rep(1,n) + exp(x %*% real.beta))
    y.logi <- c()
    for(k in 1:n){
      y.logi[k] <- rbinom(1, 1, u.logi[k])
    }
    
    beta0.logi <- as.vector(glm(y.logi~x+0, family = binomial, maxit = 200)$coefficients)
   
    #tuning param：为方便对比，直接用glmnet交叉验证的最优值作为lambda
    cv.fit <- cv.glmnet(x, y.logi, family = binomial, lambda = seq(0,1,0.1))
    opt.lam <- cv.fit$lambda.min
    
    #  mm
    correct.mm[i] <- sum(abs(Logi.MM(beta0.logi, x, y.logi, opt.lam, thd1, thd2))[c(1,2,3,4)]>0)
    incorrect.mm[i] <- sum(abs(Logi.MM(beta0.logi, x, y.logi, opt.lam, thd1, thd2))[c(5,6,7,8)]>0)
    ee.mm[i] <- t(Logi.MM(beta0.logi, x, y.logi, opt.lam, thd1, thd2)-real.beta) %*% (Logi.MM(beta0.logi, x, y.logi, opt.lam, thd1, thd2)-real.beta)
    # lqa
    correct.lqa[i] <- sum(abs(Logi.LQA(beta0.logi, x, y.logi, opt.lam, thd1, thd2))[c(1,2,3,4)]>0)
    incorrect.lqa[i] <- sum(abs(Logi.LQA(beta0.logi, x, y.logi, opt.lam, thd1, thd2))[c(5,6,7,8)]>0)
    ee.lqa[i] <- t(Logi.LQA(beta0.logi, x, y.logi, opt.lam, thd1, thd2)-real.beta) %*% (Logi.LQA(beta0.logi, x, y.logi, opt.lam, thd1, thd2)-real.beta)
    # glmnet lasso
    lasso.fit <- glmnet(x, y.logi, family = binomial, lambda = opt.lam)
    correct.glmnet[i] <- sum(abs(as.vector(lasso.fit$beta))[c(1,2,3,4)]>thd1)
    incorrect.glmnet[i] <- sum(abs(as.vector(lasso.fit$beta))[c(5,6,7,8)]>thd1)
    ee.glmnet[i] <- t(as.vector(lasso.fit$beta)-real.beta) %*% (as.vector(lasso.fit$beta)-real.beta)
    # sis scad
    sis.fit <- SIS(x, y.logi, family = 'binomial', penalty = 'SCAD')
    beta.sis <- as.vector(sis.fit$coef.est)
    which.x <- as.vector(sis.fit$ix)
    beta.t <- rep(0,p)
    n.beta <- length(beta.sis)
    for(uu in 1:n.beta){
      beta.t[which.x[uu]] <- beta.sis[uu]
    }
    beta.t <- matrix(beta.t, p, 1)
    
    correct.sis[i] <- sum(real.beta[sis.fit$ix]>0)
    incorrect.sis[i] <- sum(real.beta[-sis.fit$ix]>0)
    ee.sis <- t(beta.t-real.beta) %*% (beta.t-real.beta)
  }
  
  # 把结果存放在dataframe里做成表格展示
  result.co <- colMeans(data.frame(correct.lqa, correct.mm, correct.sis, correct.glmnet))
  result.in <- colMeans(data.frame(incorrect.lqa, incorrect.mm, incorrect.sis, incorrect.glmnet))
  result.aee <- colMeans(data.frame(ee.lqa, ee.mm, ee.sis, ee.glmnet))
  result <- rbind.data.frame(result.co,result.in,result.aee)
  colnames(result) <- c("logi lqa", "logi mm", "sis scad", "glmnet lasso")
  rownames(result) <- c("correct", "incorrect", "AEE")
  return(result)
}

simu.pois <- function(sig, real.beta)
{
  # 按协方差矩阵的改变来模拟不同情况
  # 其他参数指定为全局变量
  correct.lqa <- c()
  incorrect.lqa <- c()
  correct.mm <- c()
  incorrect.mm <- c()
  correct.sis <- c()
  incorrect.sis <- c()
  correct.glmnet <- c()
  incorrect.glmnet <- c()
  ee.mm <- c()
  ee.lqa <- c()
  ee.sis <- c()
  ee.glmnet <- c()
  
  for(i in 1:200){
    # 更新
    x <- rmvnorm(n, rep(0,p), sig)
    colnames(x) <- c("x1","x2","x3","x4","x5","x6","x7","x8")
    u.pois <- exp(x %*% real.beta)
    y.pois <- c()
    for(k in 1:n){
      y.pois[k] <- rpois(1, u.pois[k])
    }
    
    beta0.pois <- as.vector(glm(y.pois~x+0, family = poisson, maxit = 200)$coefficients)
    
    #tuning param：为方便对比，直接用glmnet交叉验证的最优值作为lambda
    cv.fit <- cv.glmnet(x, y.pois, family = poisson, lambda = seq(0,1,0.1))
    opt.lam <- cv.fit$lambda.min
    
    # pois mm
    correct.mm[i] <- sum(abs(Pois.MM(beta0.pois, x, y.pois, opt.lam, thd1, thd2))[c(1,2,3,4)]>0)
    incorrect.mm[i] <- sum(abs(Pois.MM(beta0.pois, x, y.pois, opt.lam, thd1, thd2))[c(5,6,7,8)]>0)
    ee.mm[i] <- t(Pois.MM(beta0.pois, x, y.pois, opt.lam, thd1, thd2)-real.beta) %*% (Pois.MM(beta0.pois, x, y.pois, opt.lam, thd1, thd2)-real.beta)
    # pois lqa
    correct.lqa[i] <- sum(abs(Pois.LQA(beta0.pois, x, y.pois, opt.lam, thd1, thd2))[c(1,2,3,4)]>0)
    incorrect.lqa[i] <- sum(abs(Pois.LQA(beta0.pois, x, y.pois, opt.lam, thd1, thd2))[c(5,6,7,8)]>0)
    ee.lqa[i] <- t(Pois.LQA(beta0.pois, x, y.pois, opt.lam, thd1, thd2)-real.beta) %*% (Pois.LQA(beta0.pois, x, y.pois, opt.lam, thd1, thd2)-real.beta)
    # glmnet lasso
    lasso.fit <- glmnet(x, y.pois, family = poisson, lambda = opt.lam)
    correct.glmnet[i] <- sum(abs(as.vector(lasso.fit$beta))[c(1,2,3,4)]>thd1)
    incorrect.glmnet[i] <- sum(abs(as.vector(lasso.fit$beta))[c(5,6,7,8)]>thd1)
    ee.glmnet[i] <- t(as.vector(lasso.fit$beta)-real.beta) %*% (as.vector(lasso.fit$beta)-real.beta)
    # sis scad
    sis.fit <- SIS(x, y.pois, family = 'poisson', penalty = 'SCAD')
    beta.sis <- as.vector(sis.fit$coef.est)
    which.x <- as.vector(sis.fit$ix)
    beta.t <- rep(0,p)
    n.beta <- length(beta.sis)
    for(uu in 1:n.beta){
      beta.t[which.x[uu]] <- beta.sis[uu]
    }
    beta.t <- matrix(beta.t, p, 1)
    
    correct.sis[i] <- sum(real.beta[sis.fit$ix]>0)
    incorrect.sis[i] <- sum(real.beta[-sis.fit$ix]>0)
    ee.sis[i] <- t(beta.t-real.beta) %*% (beta.t-real.beta)
  }
  
  # 把结果存放在dataframe里做成表格展示
  result.co <- colMeans(data.frame(correct.lqa, correct.mm, correct.sis, correct.glmnet))
  result.in <- colMeans(data.frame(incorrect.lqa, incorrect.mm, incorrect.sis, incorrect.glmnet))
  result.aee <- colMeans(data.frame(ee.lqa, ee.mm, ee.sis, ee.glmnet))
  result <- rbind.data.frame(result.co,result.in,result.aee)
  colnames(result) <- c("pois lqa", "pois mm", "sis scad", "glmnet lasso")
  rownames(result) <- c("correct", "incorrect", "AEE")
  return(result)
}


library(mvtnorm)         # 用于生成多元正态分布的package
library(glmnet)          # Lasso估计
library(SIS)             # SACD估计
library(MASS)

# 参数初始化
n <- 100
p <- 8
real.beta.logi <- c(1, 1, 1, 1, 0, 0, 0, 0)
real.beta.pois <- c(0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)
rho1 <- 0.3
rho2 <- 0.5
rho3 <- 0.7
thd1 <- 1e-6
thd2 <- 1e-3
tau <- 1e-4

sig1 <- diag(p)
sig.2 <- function(rho){
  return(diag(p) - diag(rep(rho,p)) + matrix(rep(rho,p*p), p, p))
}
sig.3 <- function(rho){
  s <- diag(p)
  for(i in 1:p){
    for(j in 1:p){
      s[i,j] <- rho^abs(i-j)
    }
  }
  return(s)
}


set.seed(829)


sig1.logi <- simu.logi(sig1, real.beta.logi)
sig1.pois <- simu.pois(sig1, real.beta.pois)
sig2 <- sig.2(0.3)
sig2.logi <- simu.logi(sig2, real.beta.logi)
sig2.pois <- simu.pois(sig2, real.beta.pois)
sig3 <- sig.2(0.5)
sig3.logi <- simu.logi(sig3, real.beta.logi)
sig3.pois <- simu.pois(sig3, real.beta.pois)
sig4 <- sig.2(0.7)
sig4.logi <- simu.logi(sig4, real.beta.logi)
sig4.pois <- simu.pois(sig4, real.beta.pois)
sig5 <- sig.3(0.3)
sig5.logi <- simu.logi(sig5, real.beta.logi)
sig5.pois <- simu.pois(sig5, real.beta.pois)
sig6 <- sig.3(0.5)
sig6.logi <- simu.logi(sig6, real.beta.logi)
sig6.pois <- simu.pois(sig6, real.beta.pois)
sig7 <- sig.3(0.7)
sig7.logi <- simu.logi(sig7, real.beta.logi)
sig7.pois <- simu.pois(sig7, real.beta.pois)

