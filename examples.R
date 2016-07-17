##e1s1
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

##e1s2
n <- 100
p <- 1000
gammatrue <- 1:8
betatrue <- rep(0,p)
set.seed(i+100)
U <- rbinom(8,1,0.4)
betatrue[gammatrue] <- (-1)^U*(log(n)/sqrt(n)+abs(rnorm(8)))
sigma <- 1.5
X <- matrix(rnorm(n*p),nrow=n)

Y <- X%*%betatrue+rnorm(n,0,sigma)

##e2s1
n <- 100
p <- 500
gammatrue <- 1:5
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2.5,3,-2,2.5,-2.5)
sigma <- 1.3
Sigma <- matrix(0.3,ncol = p,nrow = p)
diag(Sigma) <- 1
set.seed(i+200)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)

##e2s2
n <- 100
p <- 500
gammatrue <- 1:5
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2.5,3,-2,2.5,-2.5)
sigma <- 1.3
Sigma <- matrix(0.6,ncol = p,nrow = p)
diag(Sigma) <- 1
set.seed(i+300)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)

##e2s3
n <- 100
p <- 500
gammatrue <- 1:5
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2.5,3,-2,2.5,-2.5)
sigma <- 1.3
Sigma <- matrix(0.9,ncol = p,nrow = p)
diag(Sigma) <- 1
set.seed(i+400)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)

##e2s4
n <- 100
p <- 1000
gammatrue <- 1:5
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2.5,3,-2,2.5,-2.5)
sigma <- 1.3
Sigma <- matrix(0.3,ncol = p,nrow = p)
diag(Sigma) <- 1
set.seed(i+500)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)


##e2s5
n <- 100
p <- 1000
gammatrue <- 1:5
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2.5,3,-2,2.5,-2.5)
sigma <- 1.3
Sigma <- matrix(0.6,ncol = p,nrow = p)
diag(Sigma) <- 1
set.seed(i+600)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)


##e2s6
n <- 100
p <- 1000
gammatrue <- 1:5
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2.5,3,-2,2.5,-2.5)
sigma <- 1.3
Sigma <- matrix(0.9,ncol = p,nrow = p)
diag(Sigma) <- 1
set.seed(i+700)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)

##e3s1
n <- 100
p <- 500
gammatrue <- c(1,2,3,4,5,20,35,60,90,150,151,300)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2,-3,2,2,-3,3,-2,3,-2,3,2,-2)
sigma <- 2
Sigma <- matrix(0,ncol = p,nrow = p)
for(j in 1:p){
  for(k in 1:p){
    Sigma[j,k] <- 0.3^abs(j-k)
  }
}
set.seed(i+800)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)

##e3s2
n <- 100
p <- 500
gammatrue <- c(1,2,3,4,5,20,35,60,90,150,151,300)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2,-3,2,2,-3,3,-2,3,-2,3,2,-2)
sigma <- 2
Sigma <- matrix(0,ncol = p,nrow = p)
for(j in 1:p){
  for(k in 1:p){
    Sigma[j,k] <- 0.6^abs(j-k)
  }
}
set.seed(i+900)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)


##e3s3
n <- 100
p <- 500
gammatrue <- c(1,2,3,4,5,20,35,60,90,150,151,300)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2,-3,2,2,-3,3,-2,3,-2,3,2,-2)
sigma <- 2
Sigma <- matrix(0,ncol = p,nrow = p)
for(j in 1:p){
  for(k in 1:p){
    Sigma[j,k] <- 0.9^abs(j-k)
  }
}
set.seed(i+1000)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)

##e3s4
n <- 100
p <- 1000
gammatrue <- c(1,2,3,4,5,20,35,60,90,150,151,300)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2,-3,2,2,-3,3,-2,3,-2,3,2,-2)
sigma <- 2
Sigma <- matrix(0,ncol = p,nrow = p)
for(j in 1:p){
  for(k in 1:p){
    Sigma[j,k] <- 0.3^abs(j-k)
  }
}
set.seed(i+1100)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)

##e3s5
n <- 100
p <- 1000
gammatrue <- c(1,2,3,4,5,20,35,60,90,150,151,300)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2,-3,2,2,-3,3,-2,3,-2,3,2,-2)
sigma <- 2
Sigma <- matrix(0,ncol = p,nrow = p)
for(j in 1:p){
  for(k in 1:p){
    Sigma[j,k] <- 0.6^abs(j-k)
  }
}
set.seed(i+1200)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)


##e3s6
n <- 100
p <- 1000
gammatrue <- c(1,2,3,4,5,20,35,60,90,150,151,300)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(2,-3,2,2,-3,3,-2,3,-2,3,2,-2)
sigma <- 2
Sigma <- matrix(0,ncol = p,nrow = p)
for(j in 1:p){
  for(k in 1:p){
    Sigma[j,k] <- 0.9^abs(j-k)
  }
}
set.seed(i+1300)
X <- rmvnorm(n,mean=rep(0,p),sigma=Sigma)

Y <- X%*%betatrue+rnorm(n,0,sigma)


##e4s1
n <- 100
p <- 500
X <- matrix(NA,nrow=n,ncol=15)
set.seed(i+1400)
Z <- matrix(rnorm(n*15),nrow=n)
Zstar <- rnorm(n,0,1)
vec <- c(1,3,5,8,9,10,12,13,14,15)
for (j in vec){
  X[,j] <- 2*Z[,j]+0.1*Zstar
}
X[,2] <- 0.3*X[,1]+2*Z[,2]
X[,4] <- 0.3*X[,3]+2*Z[,4]
X[,6] <- 0.3*X[,5]+2*Z[,6]
X[,7] <- 0.4*(X[,8]+X[,9]-X[,10])+2.5*Z[,7]
X[,11] <- 0.4*(X[,14]+X[,15]-X[,12]-X[,13])+2.5*Z[,11]

X2 <- matrix(rnorm(n*(p-15)),nrow=n)

X <- cbind(X,X2)

gammatrue <- c(1,3,5,7,8,11,12,13)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(1.5,1.5,1.5,1.5,-1.5,1.5,1.5,1.5)
Y <- X%*%betatrue+rnorm(n,0,1)


##e4s2
n <- 100
p <- 500
X <- matrix(NA,nrow=n,ncol=15)
set.seed(i+1500)
Z <- matrix(rnorm(n*15),nrow=n)
Zstar <- rnorm(n,0,1)
vec <- c(1,3,5,8,9,10,12,13,14,15)
for (j in vec){
  X[,j] <- 2*Z[,j]+0.1*Zstar
}
X[,2] <- 0.6*X[,1]+1.8*Z[,2]
X[,4] <- 0.6*X[,3]+1.8*Z[,4]
X[,6] <- 0.6*X[,5]+1.8*Z[,6]
X[,7] <- 0.7*(X[,8]+X[,9]-X[,10])+1.8*Z[,7]
X[,11] <- 0.7*(X[,14]+X[,15]-X[,12]-X[,13])+1.8*Z[,11]

X2 <- matrix(rnorm(n*(p-15)),nrow=n)

X <- cbind(X,X2)

gammatrue <- c(1,3,5,7,8,11,12,13)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(1.5,1.5,1.5,1.5,-1.5,1.5,1.5,1.5)
Y <- X%*%betatrue+rnorm(n,0,1)


##e4s3
n <- 100
p <- 500
X <- matrix(NA,nrow=n,ncol=15)
set.seed(i+1600)
Z <- matrix(rnorm(n*15),nrow=n)
Zstar <- rnorm(n,0,1)
vec <- c(1,3,5,8,9,10,12,13,14,15)
for (j in vec){
  X[,j] <- 2*Z[,j]+0.5*Zstar
}
X[,2] <- X[,1]+Z[,2]
X[,4] <- X[,3]+Z[,4]
X[,6] <- X[,5]+Z[,6]
X[,7] <- X[,8]+X[,9]-X[,10]+0.1*Z[,7]
X[,11] <- X[,14]+X[,15]-X[,12]-X[,13]+0.1*Z[,11]

X2 <- matrix(rnorm(n*(p-15)),nrow=n)

X <- cbind(X,X2)

gammatrue <- c(1,3,5,7,8,11,12,13)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(1.5,1.5,1.5,1.5,-1.5,1.5,1.5,1.5)
Y <- X%*%betatrue+rnorm(n,0,1)


##e4s4
n <- 100
p <- 1000
X <- matrix(NA,nrow=n,ncol=15)
set.seed(i+1700)
Z <- matrix(rnorm(n*15),nrow=n)
Zstar <- rnorm(n,0,1)
vec <- c(1,3,5,8,9,10,12,13,14,15)
for (j in vec){
  X[,j] <- 2*Z[,j]+0.1*Zstar
}
X[,2] <- 0.3*X[,1]+2*Z[,2]
X[,4] <- 0.3*X[,3]+2*Z[,4]
X[,6] <- 0.3*X[,5]+2*Z[,6]
X[,7] <- 0.4*(X[,8]+X[,9]-X[,10])+2.5*Z[,7]
X[,11] <- 0.4*(X[,14]+X[,15]-X[,12]-X[,13])+2.5*Z[,11]

X2 <- matrix(rnorm(n*(p-15)),nrow=n)

X <- cbind(X,X2)

gammatrue <- c(1,3,5,7,8,11,12,13)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(1.5,1.5,1.5,1.5,-1.5,1.5,1.5,1.5)
Y <- X%*%betatrue+rnorm(n,0,1)


##e4s5
n <- 100
p <- 1000
X <- matrix(NA,nrow=n,ncol=15)
set.seed(i+1800)
Z <- matrix(rnorm(n*15),nrow=n)
Zstar <- rnorm(n,0,1)
vec <- c(1,3,5,8,9,10,12,13,14,15)
for (j in vec){
  X[,j] <- 2*Z[,j]+0.1*Zstar
}
X[,2] <- 0.6*X[,1]+1.8*Z[,2]
X[,4] <- 0.6*X[,3]+1.8*Z[,4]
X[,6] <- 0.6*X[,5]+1.8*Z[,6]
X[,7] <- 0.7*(X[,8]+X[,9]-X[,10])+1.8*Z[,7]
X[,11] <- 0.7*(X[,14]+X[,15]-X[,12]-X[,13])+1.8*Z[,11]

X2 <- matrix(rnorm(n*(p-15)),nrow=n)

X <- cbind(X,X2)

gammatrue <- c(1,3,5,7,8,11,12,13)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(1.5,1.5,1.5,1.5,-1.5,1.5,1.5,1.5)
Y <- X%*%betatrue+rnorm(n,0,1)


##e4s6
n <- 100
p <- 1000
X <- matrix(NA,nrow=n,ncol=15)
set.seed(i+1900)
Z <- matrix(rnorm(n*15),nrow=n)
Zstar <- rnorm(n,0,1)
vec <- c(1,3,5,8,9,10,12,13,14,15)
for (j in vec){
  X[,j] <- 2*Z[,j]+0.5*Zstar
}
X[,2] <- X[,1]+Z[,2]
X[,4] <- X[,3]+Z[,4]
X[,6] <- X[,5]+Z[,6]
X[,7] <- X[,8]+X[,9]-X[,10]+0.1*Z[,7]
X[,11] <- X[,14]+X[,15]-X[,12]-X[,13]+0.1*Z[,11]

X2 <- matrix(rnorm(n*(p-15)),nrow=n)

X <- cbind(X,X2)

gammatrue <- c(1,3,5,7,8,11,12,13)
betatrue <- rep(0,p)
betatrue[gammatrue] <- c(1.5,1.5,1.5,1.5,-1.5,1.5,1.5,1.5)
Y <- X%*%betatrue+rnorm(n,0,1)