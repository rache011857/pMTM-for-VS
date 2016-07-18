n <- 500
p <- 1000

sigma <- 1
snr <- 3
beta <- snr*sigma*sqrt(log(p)/n)*c(2,-3,2,2,-3,3,-2,3,-2,3,rep(0,p-10))

X <- matrix(rnorm(n*p), nrow=n)
Y <- X%*%beta + rnorm(n,0,sigma)

x <- scale(X)
y <- Y-mean(Y)
y.norm <- sum(y^2)

library(doMC)
registerDoMC(8)

fit <- foreach(i=1:50) %dopar% {
fit1 <- pMTM(X = X, Y = Y, s0=100, init = 'null', prior = 'p', g=p^3)
fit2 <- pMTM(X = X, Y = Y, s0=100, init = 'true', prior = 'p', g=p^3)
list(fit1,fit2)
}


pdf("pMTM2.pdf")

plot(c(0,10000),c(-6800,-4800),ann=F,bty='n',type='n')

for (i in 1:50){
  lines(fit[[i]][[1]]$log.postc.tr, col='darkolivegreen3', type='l')
  lines(fit[[i]][[2]]$log.postc.tr, col='darkorchid3', type='l')
}

abline(h=logPostc(gamma=1:10, y=y, x=x, prior='p', y.norm=y.norm, g=p^3), col='red')
abline(h=logPostc(gamma=integer(0), y=y, x=x, prior='p', y.norm=y.norm, g=p^3), col='blue')
title(main='pMTM(dataset2)')
dev.off()


plot(c(0,400),c(-6800,-4800),ann=F,bty='n',type='n')

for (i in 1:50){
  lines(fit[[i]][[1]]$log.postc.tr[1:400], col='grey', type='l')
  lines(fit[[i]][[2]]$log.postc.tr[1:400], col='grey', type='l')
}




lp <- -4900
while (lp < -4850){
  X <- matrix(rnorm(n*p), nrow=n)
  Y <- X%*%beta + rnorm(n,0,sigma)
  
  x <- scale(X)
  y <- Y-mean(Y)
  y.norm <- sum(y^2)
  lp <- logPostc(gamma=1:10, y = y, x=x, y.norm = y.norm, prior='p', g = p^3)
  print(lp)
}
