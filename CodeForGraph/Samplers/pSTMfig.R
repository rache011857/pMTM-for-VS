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
  
  gamma.abs <- length(gamma)
  log.ml.cur <- logMl(gamma = gamma, y = y, x = x, y.norm = y.norm, g = g)
  
  buildModelDf <- function(iter) {
    char <- paste(gamma,collapse = " ")
    index <- which(model.store$gamma==char)
    if (length(index)) {model.store[index,2] <<- model.store[index,2] + 1} else {
      len <- nrow(model.store)
      model.store[len+1,1] <<- char
      model.store[len+1,2] <<- 1
      model.store[len+1,3] <<- iter
    }
  }
  
  tic <- proc.time()
  for (iter in 1:n.iter){
    move.type <-  sample(3,1,prob=move.prob[gamma.abs+1,])
    
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
      }
    }
    buildModelDf(iter)
  }
  toc <- proc.time()  
  return(toc-tic)
}


