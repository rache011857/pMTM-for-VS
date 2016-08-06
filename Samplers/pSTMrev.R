library(plyr)

pSTMrev <- function(X, Y, s0, g = nrow(X), n.iter = 1e4, burnin = 2000, prior){
  
  n <- nrow(X)
  p <- ncol(X)
  
  x <- scale(X)
  y <- Y-mean(Y)
  y.norm <- sum(y^2)
  
  gamma.full <- 1:p
  gamma <- integer(0)
  
  move.prob <- matrix(1/3, s0, 3)
  move.prob[1,] <- c(1,0,0)
  move.prob[s0,] <- c(0,1,0)
  
  n.acpt <- rep(0,3)
  n.prop <- rep(0,3)
  gamma.abs <- length(gamma)
  gamma.store <- vector('list', n.iter)
  model.size <- rep(NA, n.iter)
  log.ml.cur <- logMl(gamma = gamma, y = y, x = x, y.norm = y.norm, g = g)
  
  tic <- proc.time()
  for (iter in 1:n.iter){
    move.type <-  sample(3,1,prob=move.prob[gamma.abs+1,])
    n.prop[move.type] <- n.prop[move.type] + 1
    
    if (move.type==1){ ## add
      nbhd <- gamma.full[!gamma.full %in% gamma]
      flip.ix <- 'if'(length(nbhd)==1, nbhd, sample(nbhd,1))
      gamma.prime <- c(gamma,flip.ix)
      log.ml.prop <- logMl(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, g = g)
      fwd.prob <- move.prob[gamma.abs+1,1]/(p-gamma.abs)
      bwd.prob <- move.prob[gamma.abs+2,2]/(gamma.abs+1)
      acpt.rate <- log(bwd.prob) - log(fwd.prob) + log.ml.prop - log.ml.cur + logPrior(prior=prior, move.type=1, p=p, gamma.abs=gamma.abs)
      
      if (log(runif(1)) < acpt.rate){
        log.ml.cur <- log.ml.prop
        gamma <- sort(gamma.prime)
        n.acpt[1] <- n.acpt[1] + 1
        gamma.abs <- gamma.abs + 1
      }
    } else if (move.type==2){ ## remove
      flip.ix <- 'if'(gamma.abs==1, gamma, sample(gamma,1))
      gamma.prime <- gamma[gamma!=flip.ix]
      log.ml.prop <- logMl(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, g = g)
      fwd.prob <- move.prob[gamma.abs+1,2]/gamma.abs
      bwd.prob <- move.prob[gamma.abs,1]/(p-gamma.abs+1)
      acpt.rate <- log(bwd.prob) - log(fwd.prob) + log.ml.prop - log.ml.cur + logPrior(prior=prior, move.type=-1, p=p, gamma.abs=gamma.abs)
      
      if (log(runif(1)) < acpt.rate){
        log.ml.cur <- log.ml.prop
        gamma <- gamma.prime
        n.acpt[2] <- n.acpt[2] + 1
        gamma.abs <- gamma.abs - 1
      }
    } else {
      flip.ix.rem <- 'if'(gamma.abs==1, gamma, sample(gamma,1))
      nbhd <- gamma.full[!gamma.full %in% gamma]
      flip.ix.add <- 'if'(length(nbhd)==1, nbhd, sample(nbhd,1))
      gamma.prime <- c(gamma[gamma!=flip.ix.rem],flip.ix.add)
      log.ml.prop <- logMl(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, g = g)
      
      acpt.rate <- log.ml.prop - log.ml.cur
      
      if (log(runif(1)) < acpt.rate){
        log.ml.cur <- log.ml.prop
        gamma <- sort(gamma.prime)
        n.acpt[3] <- n.acpt[3] + 1
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














