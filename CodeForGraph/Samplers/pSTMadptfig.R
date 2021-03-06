pSTMadpt <- function(X, Y, s0, zeta=2/3, g = nrow(X), cor.bound = 0.75, n.iter = 1e4, burnin = 2000, prior){
  
  n <- nrow(X)
  p <- ncol(X)
  
  x <- scale(X)
  y <- Y-mean(Y)
  y.norm <- sum(y^2)
  
  x.cor <- cor(X)
  cor.bound <- min(1, max(cor.bound, 1-1/p)) ## p>=2
  cor.cut <- as.numeric(quantile(abs(x.cor[upper.tri(x.cor)]), cor.bound))
  
  nhbr <- lapply(1:p, function(i) which(abs(x.cor[,i]) > cor.cut))
  nhbr.val <- lapply(1:p, function(i) return(fnB(abs(x.cor[nhbr[[i]],i]), cor.cut)))
  
  wt.update <- sapply(1:p, function(i) wtUpdate(p, nhbr[[i]], nhbr.val[[i]]))
  
  impt.prob <- rep(1,p)
  
  gamma.full <- 1:p
  gamma <- integer(0)
  
  move.prob <- matrix(1/3, s0, 3)
  move.prob[1,] <- c(1,0,0)
  move.prob[s0,] <- c(0,1,0)
  
  gamma.abs <- length(gamma)
  log.ml.cur <- logMl(gamma = gamma, y = y, x = x, y.norm = y.norm, g = g)
  indicator <- ifelse(1:n.iter<=burnin, 1:n.iter/burnin, (abs(1:n.iter-burnin))^(-zeta))
  
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
      weight <- impt.prob[nbhd]/sum(impt.prob[nbhd])
      flip.ix <- sample(nbhd,1,prob = weight)
      gamma.prime <- c(gamma,flip.ix)
      log.ml.prop <- logMl(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, g = g)
      fwd.prob <- move.prob[gamma.abs+1,1]*weight[which(nbhd==flip.ix)]
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
      nbhd <- c(gamma.full[!gamma.full %in% gamma],flip.ix)
      weight <- impt.prob[nbhd]/sum(impt.prob[nbhd])
      bwd.prob <- move.prob[gamma.abs,1]*weight[p-gamma.abs+1]
      acpt.rate <- log(bwd.prob) - log(fwd.prob) + log.ml.prop - log.ml.cur + logPrior(prior=prior, move.type=-1, p=p, gamma.abs=gamma.abs)
      
      if (log(runif(1)) < acpt.rate){
        log.ml.cur <- log.ml.prop
        gamma <- gamma.prime
        gamma.abs <- gamma.abs - 1
      }
    } else {
      flip.ix.rem <- 'if'(gamma.abs==1, gamma, sample(gamma,1))
      nbhd.fwd <- gamma.full[!gamma.full %in% gamma]
      weight.fwd <- impt.prob[nbhd.fwd]/sum(impt.prob[nbhd.fwd])
      flip.ix.add <- sample(nbhd.fwd,1,prob = weight.fwd)
      gamma.prime <- c(gamma[gamma!=flip.ix.rem],flip.ix.add)
      log.ml.prop <- logMl(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, g = g)
      fwd.prob <- weight.fwd[which(nbhd.fwd==flip.ix.add)]
      
      nbhd.bwd <- gamma.full[!gamma.full %in% gamma.prime]
      weight.bwd <- impt.prob[nbhd.bwd]/sum(impt.prob[nbhd.bwd])
      bwd.prob <- weight.bwd[which(nbhd.bwd==flip.ix.rem)]
      
      acpt.rate <- log(bwd.prob) - log(fwd.prob) + log.ml.prop - log.ml.cur
      
      if (log(runif(1)) < acpt.rate){
        log.ml.cur <- log.ml.prop
        gamma <- sort(gamma.prime)
      }
    }
    for(index in gamma) impt.prob <- impt.prob + indicator[iter]*wt.update[index,]
    buildModelDf(iter)
  }
  toc <- proc.time()
  
  return(toc-tic)
}
