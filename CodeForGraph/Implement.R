library(parallel)
library(compiler)
library(foreach)
library(doMC)
enableJIT(3)
n.le <- 0
pMTMadpt <- cmpfun(pMTMadpt)
pMTM <- cmpfun(pMTM)
logMl <- cmpfun(logMl)
registerDoMC(5)
fit <- vector('list',5)
M <- c(100,500,1000,4000,8000)
n.iter <- c(8e6,1.8e6,8.5e5,3.5e5,2.2e5)

# WARNING: typically take 20 or more hours
fit <- foreach(ppp=1:5) %dopar% {
  fit <- pMTMadpt(X=X, Y=Y, s0=100, M=M[ppp], g = n, n.iter = n.iter[ppp], burnin = 0, prior='bb')
  fit
}


#################################################################
# run on another machine 


# n.le <- 0
# fit <- vector('list',5)
# M <- c(100,500,1000,4000,8000)
# n.iter <- c(7.5e6,1.8e6,8e5,3e5,1.2e5)
# fit <- foreach(ppp=1:5) %dopar% {
#   fit <- pMTM(X=X, Y=Y, s0=100, M=M[ppp], g = n, n.iter = n.iter[ppp], burnin = 0, prior='bb')
#   fit
# }
