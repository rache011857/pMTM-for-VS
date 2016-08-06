library(plyr)

pMTM <- function(X, Y, s0, g=nrow(X), M = 5, n.iter = 1e4, burnin = 2000, prior){
  n <- nrow(X)
  p <- ncol(X)
  
  x <- scale(X)
  y <- Y - mean(Y)
  y.norm <- sum(y^2)
  
  gamma <- integer(0)
  n.acpt <- rep(0,3) 
  n.prop <- rep(0,3) 
  
  impt.prob <- min(1,M/p)
  
  gamma.abs <- length(gamma)
  ml.cur <- logMl(gamma=gamma, y=y, x=x, y.norm=y.norm, g=g)
  
  move.prob <- matrix(1/3, s0, 3)
  move.prob[1,] <- c(1,0,0)
  move.prob[s0,] <- c(0,1,0)
  
  model.size <- rep(NA,n.iter)
  gamma.store <- vector('list', n.iter)
  
  tic <- proc.time()
  for (iter in 1:n.iter){
    move.type <-  sample(3,1,prob=move.prob[gamma.abs+1,])
    n.prop[move.type] <- n.prop[move.type] + 1
    
    if (move.type==1){ ## add
      weight <- rep(impt.prob,p)
      weight[gamma] <- 0
      
      eta <-   which(as.logical(rbinom(p,1,weight)))
      n.fwd <- length(eta)
      if (n.fwd>0){
        gamma.tilde <- lapply(eta, function(j) return(c(gamma, j)))
        fwd.nbhd.ls <- sapply(gamma.tilde, logMl, y=y, x=x, y.norm=y.norm, g=g)
        fwd.lp <- logSum(fwd.nbhd.ls)
        fwd.ix <- sample(n.fwd, 1, prob = exp(fwd.nbhd.ls - fwd.lp))
        gamma.prime <- gamma.tilde[[fwd.ix]]
        fwd.prob <- move.prob[gamma.abs+1,1] * impt.prob
        
        if (gamma.abs > 0){
          bwd.gamma.tilde <- lapply(1:gamma.abs, function(j) return(gamma.prime[-j]))
          bwd.nbhd.ls <- c(sapply(bwd.gamma.tilde, logMl, y=y, x=x, y.norm=y.norm, g=g), ml.cur)
          bwd.lp <- logSum(bwd.nbhd.ls)
        } else {
          bwd.lp <- ml.cur
        }
        bwd.prob <- move.prob[gamma.abs+2,2]
        
        
        acpt.rate <- fwd.lp - bwd.lp + log(bwd.prob) - log(fwd.prob) + logPrior(prior=prior, move.type=1, p=p, gamma.abs=gamma.abs)
        if (log(runif(1))<acpt.rate) {
          gamma <- sort(gamma.prime)
          n.acpt[1] <- n.acpt[1] + 1
          gamma.abs <- gamma.abs + 1
          ml.cur <- fwd.nbhd.ls[fwd.ix]
        }
      }
    } else if (move.type==2){ ## remove
      weight <- rep(impt.prob,p)
      
      gamma.tilde <- lapply(1:gamma.abs, function(j) return(gamma[-j]))
      fwd.nbhd.ls <- sapply(gamma.tilde, logMl, y=y, x=x, y.norm=y.norm, g=g)
      fwd.lp <- logSum(fwd.nbhd.ls)
      fwd.ix <- sample(gamma.abs, 1, prob = exp(fwd.nbhd.ls - fwd.lp))
      gamma.prime <- gamma.tilde[[fwd.ix]]
      fwd.prob <- move.prob[gamma.abs+1,2]
      
      bwd.prob <- move.prob[gamma.abs,1] * impt.prob
      weight[gamma] <- 0
      eta <- which(as.logical(rbinom(p,1,weight)))
      
      if (length(eta)>0){
        bwd.gamma.tilde <- lapply(eta, function(j) return(c(gamma.prime, j)))
        bwd.nbhd.ls <- c(sapply(bwd.gamma.tilde, logMl, y=y, x=x, y.norm=y.norm, g=g), ml.cur)
        bwd.lp <- logSum(bwd.nbhd.ls)
      } else {
        bwd.lp <- ml.cur
      }
      
      acpt.rate <- fwd.lp - bwd.lp + log(bwd.prob) - log(fwd.prob) + logPrior(prior=prior, move.type=-1, p=p, gamma.abs=gamma.abs)
      if (log(runif(1))<acpt.rate) {
        gamma <- gamma.prime
        n.acpt[2] <- n.acpt[2] + 1
        gamma.abs <- gamma.abs - 1
        ml.cur <- fwd.nbhd.ls[fwd.ix]
      }
    } else { ## swap
      weight <- rep(min(1, M/p/gamma.abs), p)
      weight[gamma] <- 0
      fwd.var.add <- integer(0)
      fwd.ix.rem <- integer(0)
      n.fwd <- 0
      for(ix in 1:gamma.abs){
        eta.add <- which(as.logical(rbinom(p,1,weight)))
        n.fwd.temp <- length(eta.add)
        fwd.var.add <- c(fwd.var.add, eta.add)
        fwd.ix.rem <- c(fwd.ix.rem, rep(ix, n.fwd.temp))
        n.fwd <- n.fwd + n.fwd.temp
      }
      
      if(n.fwd > 0){
        gamma.tilde <- lapply(1:n.fwd, function(j) return(c(gamma[-fwd.ix.rem[j]], fwd.var.add[j])))
        fwd.nbhd.ls <- sapply(gamma.tilde, logMl, y=y, x=x, y.norm=y.norm, g=g)
        fwd.lp <- logSum(fwd.nbhd.ls)
        fwd.ix <- sample(n.fwd, 1, prob = exp(fwd.nbhd.ls - fwd.lp))
        gamma.prime <- gamma.tilde[[fwd.ix]]
        var.add <- fwd.var.add[fwd.ix]
        var.rem <- gamma[fwd.ix.rem[fwd.ix]]
        
        weight[var.add] <- 0
        weight[var.rem] <- min(1, M/p/gamma.abs)
        bwd.var.add <- integer(0)
        bwd.ix.rem <- integer(0)
        n.bwd <- 0
        
        if (gamma.abs > 1){
          for(ix in 1:(gamma.abs-1)){
            eta.add <- which(as.logical(rbinom(p,1,weight)))
            n.bwd.temp <- length(eta.add)
            bwd.var.add <- c(bwd.var.add, eta.add)
            bwd.ix.rem <- c(bwd.ix.rem, rep(ix, n.bwd.temp))
            n.bwd <- n.bwd + n.bwd.temp
          }
        }
        weight[var.rem] <- 0
        eta.add <- which(as.logical(rbinom(p,1,weight)))
        n.bwd.temp <- length(eta.add)
        bwd.var.add <- c(bwd.var.add, eta.add)
        bwd.ix.rem <- c(bwd.ix.rem, rep(gamma.abs, n.bwd.temp))
        n.bwd <- n.bwd + n.bwd.temp
        
        if (n.bwd>0){
          bwd.gamma.tilde <- lapply(1:n.bwd, function(j) return(c(gamma.prime[-bwd.ix.rem[j]], bwd.var.add[j])))
          bwd.nbhd.ls <- c(sapply(bwd.gamma.tilde, logMl, y=y, x=x, y.norm=y.norm, g=g), ml.cur)
          bwd.lp <- logSum(bwd.nbhd.ls)
        } else {
          bwd.lp <- ml.cur
        }
        
        acpt.rate <- fwd.lp - bwd.lp
        if (log(runif(1))<acpt.rate) {
          gamma <- sort(gamma.prime)
          n.acpt[3] <- n.acpt[3] + 1
          ml.cur <- fwd.nbhd.ls[fwd.ix]
        }
      }
    }
    gamma.store[[iter]] <- gamma
    model.size[iter] <- gamma.abs
  }
  toc <- proc.time()
  inclusion <- as.matrix(count(unlist(gamma.store[-(1:burnin)])))
  MPM <- inclusion[which(inclusion[,2]>=(n.iter-burnin)/2),1]
  model.count <- count(sapply(gamma.store[-(1:burnin)], paste, collapse=" "))
  model.sort <- model.count[order(-model.count[2]),]
  post.prob <- model.sort[[2]]/(n.iter-burnin)
  model <- lapply(model.sort$x,modelSplit)
  HPM <- model[[1]]
  model.size.avg <- mean(model.size[-(1:burnin)])
  model.size.HPM <- length(MPM)
  model.size.MPM <- length(HPM)
  
  return(list(model=model, post.prob=post.prob, time.spend=toc-tic, n.prop=n.prop, n.acpt=n.acpt, n.iter=n.iter, MPM=MPM, HPM=HPM, burnin=burnin, model.size.avg=model.size.avg,model.size.MPM=model.size.MPM,model.size.HPM=model.size.HPM))
}