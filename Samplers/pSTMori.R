library(plyr)

pSTMori <- function(X, Y, s0, g = nrow(X), n.iter = 1e4, burnin = 2000, prior){
  
  n <- nrow(X)
  p <- ncol(X)
  
  x <- scale(X)
  y <- Y-mean(Y)
  y.norm <- sum(y^2)
  
  gamma.full <- 1:p
  gamma <- integer(0)
  
  n.acpt <- rep(0,3)
  n.prop <- rep(0,3)
  
  gamma.abs <- length(gamma)
  gamma.store <- vector('list', n.iter)
  model.size <- rep(NA, n.iter)
  log.ml.cur <- logMl(gamma = gamma, y = y, x = x, y.norm = y.norm, g = g)
  
  for (iter in 1:n.iter){
    move.type <-  rbinom(1,1,0.5)
    
    if (move.type==0){ ## single flip (add or remove)
      flip.ix <- sample.int(p,1)
      if (flip.ix %in% gamma){ # remove
        n.prop[1] <- n.prop[1] + 1
        gamma.prime <- gamma[gamma!=flip.ix]
        log.ml.prop <- logMl(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, g = g)
        acpt.rate <- log.ml.prop - log.ml.cur + logPrior(prior=prior, move.type=-1, p=p, gamma.abs=gamma.abs)
        
        if (log(runif(1)) < acpt.rate){
          gamma <- gamma.prime
          n.acpt[1] <- n.acpt[1] + 1
          log.ml.cur <- log.ml.prop
          gamma.abs <- gamma.abs - 1
        }
      } else {  # add
        n.prop[2] <- n.prop[2] + 1
        if (gamma.abs<=s0-1){
          gamma.prime <- c(gamma,flip.ix)
          log.ml.prop <- logMl(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, g = g)
          acpt.rate <- log.ml.prop - log.ml.cur + logPrior(prior=prior, move.type=1, p=p, gamma.abs=gamma.abs)
          
          if (log(runif(1)) < acpt.rate){
            gamma <- sort(gamma.prime)
            n.acpt[2] <- n.acpt[2] + 1
            log.ml.cur <- log.ml.prop
            gamma.abs <- gamma.abs + 1
          }
        }
      }
    } else { # swap
      n.prop[3] <- n.prop[3] + 1
      if(gamma.abs > 0){
        flip.ix.rem <- 'if'(gamma.abs==1, gamma, sample(gamma,1))
        flip.ix.add <- sample(gamma.full[!gamma.full %in% gamma],1)
        gamma.prime <- c(gamma[gamma!=flip.ix.rem],flip.ix.add)
        log.ml.prop <- logMl(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, g = g)
        acpt.rate <- log.ml.prop - log.ml.cur 
        
        if (log(runif(1)) < acpt.rate){
          gamma <- sort(gamma.prime)
          n.acpt[3] <- n.acpt[3] + 1
          log.ml.cur <- log.ml.prop
        }
      }
    }
    
    gamma.store[[iter]] <- gamma
    model.size[iter] <- gamma.abs
  }
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
