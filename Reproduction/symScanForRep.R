n <- 500
p <- 1000

sigma <- 1
snr <- 3
beta <- snr*sigma*sqrt(log(p)/n)*c(2,-3,2,2,-3,3,-2,3,-2,3,rep(0,p-10))

X <- matrix(rnorm(n*p), nrow=n)
Y <- X%*%beta + rnorm(n,0,sigma)

K <- 5
g <- p^3

Y.norm <- t(Y)%*%Y


logPostc <- function(X, Y, g=p^3, s0=100, gamma){
  gamma <- which(gamma==1)
  gamma.abs <- length(gamma)
  n <- nrow(X)
  p <- ncol(X)
  if (gamma.abs==0) {return (-n/2*log1p(g))
  } else if (gamma.abs==1){
    X.gamma <- X[,gamma]
    rsq.gamma <- sum(Y*X.gamma)^2/sum(X.gamma^2)/Y.norm
    return(-3*gamma.abs*log(p)-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
  } else{
    X.gamma <- X[,gamma]
    cp <- crossprod(X.gamma,Y)
    rsq.gamma <- crossprod(cp,solve(crossprod(X.gamma),cp))/Y.norm
    return(-3*gamma.abs*log(p)-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
  }
}


symScan <- function(X, Y, K, init, s0=100, g=p^3, sym = TRUE){
  p <- ncol(X)
  gamma <- rep(0,p)
  if (init=='null'){
    gamma[sample.int(p-10,50,replace=F)+10] <- 1
  } else {
    gamma[c(1:10,sample.int(p-10,40,replace=F)+10)] <- 1
  }
  lp.cur <- logPostc(X,Y,g,s0,gamma)
  logPostc.tr <- rep(NA, K*p)
  for (k in 1:K){
    gamma.full <- sample(p,p, replace = !sym)
    for (j in 1:p){
      gamma.prime <- gamma
      gamma.prime[gamma.full[j]] <- 1 - gamma[gamma.full[j]]
      lp.prop <- logPostc(X,Y,g,s0,gamma.prime)
      if (log(runif(1)) < lp.prop - lp.cur) {
        gamma <- gamma.prime
        lp.cur <- lp.prop
      }
      logPostc.tr[(k-1)*p+j] <- lp.cur
    }
  }
  return(logPostc.tr)
}
  

lp <- matrix(NA,nrow=M,ncol=K*p)
for (m in 1:50) {
  lp[m,] <- symScan(X=X, Y=Y, K=K, init='null', sym=F)
  lp[m+50,] <- symScan(X=X, Y=Y, K=K, init='true', sym=F)
}

plot(lp[1,],type='l',col='grey')
for (i in 2:100){
  lines(lp[i,],type='l',col='grey')
}
abline(h=logPostc(X,Y,g,s0,rep(0,p)),col='blue',pch=2)
abline(h=logPostc(X,Y,g,s0,c(rep(1,10),rep(0,p-10))),col='red',pch=2)
