library(mvtnorm)
library(glmnet)
library(ncvreg)
library(parcor)
library(EMVS)
library(foreach)
library(doMC)
registerDoMC(10)  


p=500
v0=seq(0.1,1,by=0.1)
v1=1000
beta_init=rep(0,p)
a=10
b=p-10
epsilon=10^{-5}

e1s1 <- foreach (i=1:100) %dopar% {
  n <- 100
  p <- 500
  gammatrue <- 1:8
  betatrue <- rep(0,p)
  set.seed(i)
  U <- rbinom(8,1,0.4)
  betatrue[gammatrue] <- (-1)^U*(log(n)/sqrt(n)+abs(rnorm(8)))
  sigma <- 1.5
  X <- matrix(rnorm(n*p),nrow=n)
  
  Y <- X%*%betatrue+rnorm(n,0,sigma)
  
  x <- scale(X)
  sdx <- attr(x,"scaled:scale")
  y <- Y-mean(Y)
  
  fit <- vector('list',7)
  fit[[1]] <- benchmark(X=X, Y=Y, s0=100, g = n, n.iter = 1e4, burnin = 2000, prior='bb')
  fit[[2]] <- pMTM.noc(X=X, Y=Y, s0=100, g = n, n.iter = 1e4, burnin = 2000, prior='bb')
  fit[[3]] <- pMTM(X=X, Y=Y, s0=100, g = n, n.iter = 1e4, burnin = 2000, prior='bb')
  
  fit[[4]] <- cv.glmnet(x,y,family='gaussian',alpha=1,intercept=F)
  fit[[5]] <- adalasso(x,y,use.Gram = F,intercept = F)
  fit[[6]] <- EMVS(y,x,v0=v0,v1=v1,beta_init = beta_init,type = 'betabinomial',sigma_init = 1,epsilon = epsilon,a=a,b=b)
  fit[[7]] <- cv.ncvreg(x,y)
  fit[[8]] <- cv.ncvreg(x,y,penalty='SCAD')
  
  beta <- matrix(NA,nrow=14,ncol=p)
  beta[7,] <- coef(fit[[4]],lambda='lambda.min')[-1]/sdx
  beta[8,] <- fit[[5]]$coefficients.adalasso/sdx
  beta[9,] <- emvs.best(fit[[6]],x=x,y=y)
  beta[10,] <- fit[[7]]$fit$beta[,fit[[7]]$min][-1]/sdx
  beta[11,] <- fit[[8]]$fit$beta[,fit[[8]]$min][-1]/sdx
  
  gamma <- vector('list',11)
  gamma[[1]] <- fit[[1]]$MPM
  gamma[[2]] <- fit[[1]]$HPM
  gamma[[3]] <- fit[[2]]$MPM
  gamma[[4]] <- fit[[2]]$HPM
  gamma[[5]] <- fit[[3]]$MPM
  gamma[[6]] <- fit[[3]]$HPM
  gamma[[7]] <- which(beta[7,]!=0)
  gamma[[8]] <- which(beta[8,]!=0)
  gamma[[9]] <- EMVSbest(fit[[6]])$indices
  gamma[[10]] <- which(beta[10,]!=0)
  gamma[[11]] <- which(beta[11,]!=0)
  
  for (j in 1:6){
    beta[j,] <- model.coef(gamma=gamma[[j]],g=n,x=x,Y=Y,sdx=sdx)
  }
  for (j in 12:14){
    beta[j,] <- bma.coef(x=x, Y=Y, sdx=sdx, g=n, gamma.model=fit[[j-11]]$gamma.model, post.prob=fit[[j-11]]$post.prob)
  }
  
  mse <- rep(NA,14)
  for (j in 1:14){
    mse[j] <- sqrt(sum((beta[j,]-betatrue)^2))
  }
  
  modelsize <- sapply(gamma,length)
  
  fp <- fn <- fdr <- rep(NA,11)
  for (j in 1:11){
    fp[j] <- sum(!gamma[[j]] %in% gammatrue)
    fn[j] <- sum(!gammatrue %in% gamma[[j]])
    fdr[j] <- ifelse(modelsize[j], fp[j]/modelsize[j], 0)
  }
  list(modelsize=modelsize,fp=fp,fn=fn,fdr=fdr,mse=mse)  
}