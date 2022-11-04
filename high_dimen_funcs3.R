

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

# 惩罚矩阵
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


yuce.logi <- function(sig, real.beta)
{
  # 模拟预测效果的函数
  a.lqa <- c()
  a.mm <- c()
  a.sis <- c()
  a.glmnet <- c()
  
  for(i in 1:200){
    # 训练，估计参数
    x <- rmvnorm(n, rep(0,p), sig)
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
    beta.mm <- Logi.MM(beta0.logi, x, y.logi, opt.lam, thd1, thd2)
    # lqa
    beta.lqa <- Logi.LQA(beta0.logi, x, y.logi, opt.lam, thd1, thd2)
    # glmnet lasso
    lasso.fit <- glmnet(x, y.logi, family = binomial, lambda = opt.lam)
    beta.lasso <- lasso.fit$beta
    # sis scad
    sis.fit <- SIS(x, y.logi, family = 'binomial', penalty = 'SCAD')
    beta.sis <- as.vector(sis.fit$coef.est)
    which.x <- as.vector(sis.fit$ix)
    beta.scad <- rep(0,p)
    n.beta <- length(beta.sis)
    for(uu in 1:n.beta){
      beta.scad[which.x[uu]] <- beta.sis[uu]
    }
    beta.scad <- matrix(beta.scad, p, 1)
    
    x.test <- rmvnorm(n, rep(0,p), sig)
    u.test <- exp(x.test %*% real.beta)/ (rep(1,n) + exp(x.test %*% real.beta))
    y.test <- c()
    u.lqa <- exp(x.test %*% beta.lqa)/ (rep(1,n) + exp(x.test %*% beta.lqa))
    y.lqa <- c()
    u.mm <- exp(x.test %*% beta.mm)/ (rep(1,n) + exp(x.test %*% beta.mm))
    y.mm <- c()
    u.lasso <- exp(x.test %*% beta.lasso)/ (rep(1,n) + exp(x.test %*% beta.lasso))
    y.lasso <- c()
    u.scad <- exp(x.test %*% beta.scad)/ (rep(1,n) + exp(x.test %*% beta.scad))
    y.scad <- c()
    for(k in 1:n){
      y.test[k] <- rbinom(1, 1, u.test[k])
      y.lqa[k] <- rbinom(1, 1, u.lqa[k])
      y.mm[k] <- rbinom(1, 1, u.mm[k])
      y.lasso[k] <- rbinom(1, 1, u.lasso[k])
      y.scad[k] <- rbinom(1, 1, u.scad[k])
    }
    
    a.lqa[i] <- sum(y.lqa==y.test)/n
    a.mm[i] <- sum(y.mm==y.test)/n
    a.sis[i] <- sum(y.scad==y.test)/n
    a.glmnet[i] <- sum(y.lasso==y.test)/n
  }
  
  # 把结果存放在dataframe里做成表格展示
  result <- c(mean(a.lqa), mean(a.mm), mean(a.glmnet), mean(a.sis))
  return(result)
}


yuce.pois <- function(sig, real.beta)
{
  # 按协方差矩阵的改变来模拟不同情况
  # 其他参数指定为全局变量
  mse.lqa <- c()
  mse.mm <- c()
  mse.sis <- c()
  mse.glmnet <- c()
  
  for(i in 1:200){
    # 更新
    x <- rmvnorm(n, rep(0,p), sig)
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
    beta.mm <- Pois.MM(beta0.pois, x, y.pois, opt.lam, thd1, thd2)
    # pois lqa
    beta.lqa <- Pois.LQA(beta0.pois, x, y.pois, opt.lam, thd1, thd2)
    # glmnet lasso
    lasso.fit <- glmnet(x, y.pois, family = poisson, lambda = opt.lam)
    beta.lasso <- lasso.fit$beta
    # sis scad
    sis.fit <- SIS(x, y.pois, family = 'poisson', penalty = 'SCAD')
    beta.sis <- as.vector(sis.fit$coef.est)
    which.x <- as.vector(sis.fit$ix)
    beta.scad <- rep(0,p)
    n.beta <- length(beta.sis)
    for(uu in 1:n.beta){
      beta.scad[which.x[uu]] <- beta.sis[uu]
    }
    beta.scad <- matrix(beta.scad, p, 1)
    
    x.test <- rmvnorm(n, rep(0,p), sig)
    u.test <- exp(x.test %*% real.beta)
    y.test <- c()
    u.lqa <- exp(x.test %*% beta.lqa)
    y.lqa <- c()
    u.mm <- exp(x.test %*% beta.mm)
    y.mm <- c()
    u.lasso <- exp(x.test %*% beta.lasso)
    y.lasso <- c()
    u.scad <- exp(x.test %*% beta.scad)
    y.scad <- c()
    for(k in 1:n){
      y.test[k] <- rpois(1, u.test[k])
      y.lqa[k] <- rpois(1, u.lqa[k])
      y.mm[k] <- rpois(1, u.mm[k])
      y.lasso[k] <- rpois(1, u.lasso[k])
      y.scad[k] <- rpois(1, u.scad[k])
    }
    
    mse.lqa <- t(y.lqa-y.test)%*%(y.lqa-y.test)
    mse.mm <- t(y.mm-y.test)%*%(y.mm-y.test)
    mse.sis <- t(y.scad-y.test)%*%(y.scad-y.test)
    mse.glmnet <- t(y.lasso-y.test)%*%(y.lasso-y.test)
  }
  # 把结果存放在dataframe里做成表格展示
  result <- c(mean(mse.lqa), mean(mse.mm), mean(mse.glmnet), mean(mse.sis))
  return(result)
}

library(mvtnorm)         # 用于生成多元正态分布的package
library(glmnet)          # Lasso估计
library(SIS)             # SACD估计
library(MASS)

# 参数初始化
n <- 500
p <- 40
real.beta.logi <- c(c(1, 1, 1, 1), rep(0,p-4))
real.beta.pois <- c(c(0.5, 0.5, 0.5, 0.5), rep(0,p-4))
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


a1.logi <- yuce.logi(sig1, real.beta.logi)
a1.pois <- yuce.pois(sig1, real.beta.pois)
sig2 <- sig.2(0.3)
a2.logi <- yuce.logi(sig2, real.beta.logi)
a2.pois <- yuce.pois(sig2, real.beta.pois)
sig3 <- sig.2(0.5)
a3.logi <- yuce.logi(sig3, real.beta.logi)
a3.pois <- yuce.pois(sig3, real.beta.pois)
sig4 <- sig.2(0.7)
a4.logi <- yuce.logi(sig4, real.beta.logi)
a4.pois <- yuce.pois(sig4, real.beta.pois)
sig5 <- sig.3(0.3)
a5.logi <- yuce.logi(sig5, real.beta.logi)
a5.pois <- yuce.pois(sig5, real.beta.pois)
sig6 <- sig.3(0.5)
a6.logi <- yuce.logi(sig6, real.beta.logi)
a6.pois <- yuce.pois(sig6, real.beta.pois)
sig7 <- sig.3(0.7)
a7.logi <- yuce.logi(sig7, real.beta.logi)
a7.pois <- yuce.pois(sig7, real.beta.pois)


simu.bic <- function(sig, real.beta)
{
  for(i in 1:200){
    # 更新
    x <- rmvnorm(n, rep(0,p), sig)
    u.pois <- exp(x %*% real.beta)
    y.pois <- c()
    for(k in 1:n){
      y.pois[k] <- rpois(1, u.pois[k])
    }
    beta0.pois <- as.vector(glm(y.pois~x+0, family = poisson, maxit = 200)$coefficients)
    lams <- seq(0,5,0.1)
    bic.mm <- c()
    for(lam in lams){
      beta.mm <- Pois.MM(beta0.pois, x, y.pois, lam, thd1, thd2)
      u.mm <- exp(x %*% beta.mm)
      y.mm <- c()
      for(k in 1:n){
        y.mm[k] <- rpois(1, u.mm[k])
      }
      bic.mm <- c(bic.mm, t(y.mm-y.pois)%*%(y.mm-y.pois)/n + log(n)/n * sum(beta.mm > 0))
    }
    opt.lam <- lams[which.min(bic.mm)]
    beta.mm <- Pois.MM(beta0.pois, x, y.pois, opt.lam, thd1, thd2)
    C <- sum(abs(beta.mm)[c(1,2,3,4)]>0)
    IC <- sum(abs(beta.mm)[-c(1,2,3,4)]>0)
    ee <- t(beta.mm-real.beta)%*%(beta.mm-real.beta)
  }
  # 把结果存放在dataframe里做成表格展示
  result <- c(C,IC, ee)
  return(result)
}

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

simu.bic(sig.2(0.5), real.beta.pois)




