library(mvtnorm)
n <- 50
p <- 8000
X1 <- matrix(NA,nrow=n,ncol=10)
set.seed(42)
Z <- matrix(rnorm(n*10),nrow=n)
Zstar <- rnorm(n,0,1)
vec <- c(1,3,5,8,9,10)
for (j in vec){
  X1[,j] <- 2*Z[,j]+0.1*Zstar
}
X1[,2] <- 0.6*X1[,1]+1.8*Z[,2]
X1[,4] <- 0.6*X1[,3]+1.8*Z[,4]
X1[,6] <- 0.6*X1[,5]+1.8*Z[,6]
X1[,7] <- 0.7*(X1[,8]+X1[,9]-X1[,10])+1.8*Z[,7]
Sigma <- matrix(c(1.000, -0.134, -0.283, -0.310, 0.221,
                  -0.134, 1.000, -0.125, -0.111, 0.234,
                  -0.283, -0.125, 1.000, 0.959, -0.092,
                  -0.310, -0.111, 0.959, 1.000, -0.182,
                  0.221, 0.234, -0.092, -0.182, 1.000), nrow=5)
X2 <- rmvnorm(n=n, sigma=Sigma)
X3 <- matrix(rnorm(n*(p-15)),nrow=n)
X <- cbind(X1,X2,X3)

gammatrue <- c(1,3,5,7,8,11,12,13)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(1.5,1,-1.5,-1,-0.8,0.9,0.8,1)
Y <- X%*%betatrue+rnorm(n,0,1)

x <- scale(X)
y <- Y-mean(Y)
y.norm <- sum(y^2)