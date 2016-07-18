pSTMori <- function(X, Y, s0, init, g = nrow(X), n.iter = 1e4, burnin = 2000, prior){
  
  n <- nrow(X)
  p <- ncol(X)
  
  x <- scale(X)
  y <- Y-mean(Y)
  y.norm <- sum(y^2)
  
  gamma.full <- 1:p
  gamma <- 'if'(init=='null', sample.int(p-10,50,replace=F)+10, c(1:10,sample.int(p-10,40,replace=F)+10))
  
  n.acpt <- rep(0,3)
  n.prop <- rep(0,3)
  gamma.abs <- length(gamma)
  gamma.store <- vector('list', n.iter)
  log.postc.tr <-  model.size <- rep(NA, n.iter)
  log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
  
  for (iter in 1:n.iter){
    move.type <-  rbinom(1,1,0.5)
    
    if (move.type==0){ ## single flip (add or remove)
      flip.ix <- sample.int(p,1)
      if (flip.ix %in% gamma){ # remove
        n.prop[1] <- n.prop[1] + 1
        gamma.prime <- gamma[gamma!=flip.ix]
        log.post.prop <- logPostc(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, prior = prior, g=g)
        
        if (log(runif(1)) < log.post.prop - log.post.cur){
          gamma <- gamma.prime
          n.acpt[1] <- n.acpt[1] + 1
          log.post.cur <- log.post.prop
          gamma.abs <- gamma.abs - 1
        }
      } else {  # add
        n.prop[2] <- n.prop[2] + 1
        if (gamma.abs<=s0-1){
          gamma.prime <- c(gamma,flip.ix)
          log.post.prop <- logPostc(gamma = gamma.prime, y = y, x = x, prior = prior, y.norm = y.norm, g=g)
          
          if (log(runif(1)) < log.post.prop - log.post.cur){
            gamma <- sort(gamma.prime)
            n.acpt[2] <- n.acpt[2] + 1
            log.post.cur <- log.post.prop
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
        log.post.prop <- logPostc(gamma = gamma.prime, y = y, x = x, prior = prior, y.norm = y.norm, g=g)
        
        if (log(runif(1)) < log.post.prop - log.post.cur){
          gamma <- sort(gamma.prime)
          n.acpt[3] <- n.acpt[3] + 1
          log.post.cur <- log.post.prop
        }
      }
    }
    
    gamma.store[[iter]] <- gamma
    log.postc.tr[iter] <- log.post.cur
    model.size[iter] <- gamma.abs
  }
  gamma.mat <- t(sapply(gamma.store[-(1:burnin)],inclusion, p=p))
  incl.prob <- apply(gamma.mat,2,sum)/(n.iter-burnin)
  MPM <- which(incl.prob>=0.5)
  gamma.count <- table(apply(gamma.mat, 1, paste, collapse=" "))
  gamma.sort <- sort(gamma.count,decreasing = T)
  post.prob <- as.numeric(gamma.sort)/(n.iter-burnin)
  gamma.model <- lapply(names(gamma.sort),gammaSplit)
  HPM <- gamma.model[[1]]
  model.size.avg <- mean(model.size[-(1:burnin)])
  model.size.HPM <- length(MPM)
  model.size.MPM <- length(HPM)
  return(list(gamma.model=gamma.model, post.prob=post.prob, log.postc.tr=log.postc.tr, n.prop=n.prop, n.acpt=n.acpt, n.iter=n.iter, MPM=MPM, HPM=HPM, burnin=burnin, gamma.store=gamma.store, model.size.avg=model.size.avg,model.size.MPM=model.size.MPM,model.size.HPM=model.size.HPM, incl.prob=incl.prob))
}


pSTMrev <- function(X, Y, s0, init, g = nrow(X), n.iter = 1e4, burnin = 2000, prior){
  
  n <- nrow(X)
  p <- ncol(X)
  
  x <- scale(X)
  y <- Y-mean(Y)
  y.norm <- sum(y^2)
  
  gamma.full <- 1:p
  gamma <- 'if'(init=='null', sample.int(p-10,50,replace=F)+10, c(1:10,sample.int(p-10,40,replace=F)+10))
  
  move.prob <- matrix(1/3, s0, 3)
  move.prob[1,] <- c(1,0,0)
  move.prob[s0,] <- c(0,1,0)
  
  n.acpt <- rep(0,3)
  n.prop <- rep(0,3)
  gamma.abs <- length(gamma)
  gamma.store <- vector('list', n.iter)
  log.postc.tr <-  model.size <- rep(NA, n.iter)
  log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
  
  for (iter in 1:n.iter){
    move.type <-  sample(3,1,prob=move.prob[gamma.abs+1,])
    n.prop[move.type] <- n.prop[move.type] + 1
    
    if (move.type==1){ ## add
      nbhd <- gamma.full[!gamma.full %in% gamma]
      flip.ix <- 'if'(length(nbhd)==1, nbhd, sample(nbhd,1))
      gamma.prime <- c(gamma,flip.ix)
      log.post.prop <- logPostc(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, prior = prior, g=g)
      fwd.prob <- move.prob[gamma.abs+1,1]/(p-gamma.abs)
      bwd.prob <- move.prob[gamma.abs+2,2]/(gamma.abs+1)
      acpt.rate <- log(bwd.prob) - log(fwd.prob) + log.post.prop - log.post.cur
      
      if (log(runif(1)) < acpt.rate){
        log.post.cur <- log.post.prop
        gamma <- sort(gamma.prime)
        n.acpt[1] <- n.acpt[1] + 1
        gamma.abs <- gamma.abs + 1
      }
    } else if (move.type==2){ ## remove
      flip.ix <- 'if'(gamma.abs==1, gamma, sample(gamma,1))
      gamma.prime <- gamma[gamma!=flip.ix]
      log.post.prop <- logPostc(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, prior = prior, g=g)
      fwd.prob <- move.prob[gamma.abs+1,2]/gamma.abs
      bwd.prob <- move.prob[gamma.abs,1]/(p-gamma.abs+1)
      acpt.rate <- log(bwd.prob) - log(fwd.prob) + log.post.prop - log.post.cur
      
      if (log(runif(1)) < acpt.rate){
        log.post.cur <- log.post.prop
        gamma <- gamma.prime
        n.acpt[2] <- n.acpt[2] + 1
        gamma.abs <- gamma.abs - 1
      }
    } else {
      flip.ix.rem <- 'if'(gamma.abs==1, gamma, sample(gamma,1))
      nbhd <- gamma.full[!gamma.full %in% gamma]
      flip.ix.add <- 'if'(length(nbhd)==1, nbhd, sample(nbhd,1))
      gamma.prime <- c(gamma[gamma!=flip.ix.rem],flip.ix.add)
      log.post.prop <- logPostc(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, prior = prior, g=g)
      
      acpt.rate <- log.post.prop - log.post.cur
      
      if (log(runif(1)) < acpt.rate){
        log.post.cur <- log.post.prop
        gamma <- sort(gamma.prime)
        n.acpt[3] <- n.acpt[3] + 1
      }
    }
    
    gamma.store[[iter]] <- gamma
    model.size[iter] <- gamma.abs
    log.postc.tr[iter] <- log.post.cur
  }
  gamma.mat <- t(sapply(gamma.store[-(1:burnin)],inclusion, p=p))
  incl.prob <- apply(gamma.mat,2,sum)/(n.iter-burnin)
  MPM <- which(incl.prob>=0.5)
  gamma.count <- table(apply(gamma.mat, 1, paste, collapse=" "))
  gamma.sort <- sort(gamma.count,decreasing = T)
  post.prob <- as.numeric(gamma.sort)/(n.iter-burnin)
  gamma.model <- lapply(names(gamma.sort),gammaSplit)
  HPM <- gamma.model[[1]]
  model.size.avg <- mean(model.size[-(1:burnin)])
  model.size.HPM <- length(MPM)
  model.size.MPM <- length(HPM)
  return(list(gamma.model=gamma.model, post.prob=post.prob, log.postc.tr=log.postc.tr, n.prop=n.prop, n.acpt=n.acpt, n.iter=n.iter, MPM=MPM, HPM=HPM, burnin=burnin, model.size.avg=model.size.avg,model.size.MPM=model.size.MPM,model.size.HPM=model.size.HPM))
}


pSTMadpt <- function(X, Y, s0, init, zeta=2/3, g = nrow(X), cor.bound = 0.75, n.iter = 1e4, burnin = 2000, prior){
  
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
  gamma <- 'if'(init=='null', sample.int(p-10,50,replace=F)+10, c(1:10,sample.int(p-10,40,replace=F)+10))
  
  move.prob <- matrix(1/3, s0, 3)
  move.prob[1,] <- c(1,0,0)
  move.prob[s0,] <- c(0,1,0)
  
  n.acpt <- rep(0,3) 
  n.prop <- rep(0,3)
  gamma.abs <- length(gamma)
  gamma.store <- vector('list', n.iter)
  log.postc.tr <-  model.size <- rep(NA, n.iter)
  log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
  
  for (iter in 1:n.iter){
    move.type <-  sample(3,1,prob=move.prob[gamma.abs+1,])
    n.prop[move.type] <- n.prop[move.type] + 1
    
    if (move.type==1){ ## add
      nbhd <- gamma.full[!gamma.full %in% gamma]
      weight <- impt.prob[nbhd]/sum(impt.prob[nbhd])
      flip.ix <- sample(nbhd,1,prob = weight)
      gamma.prime <- c(gamma,flip.ix)
      log.post.prop <- logPostc(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
      fwd.prob <- move.prob[gamma.abs+1,1]*weight[which(nbhd==flip.ix)]
      bwd.prob <- move.prob[gamma.abs+2,2]/(gamma.abs+1)
      acpt.rate <- log(bwd.prob) - log(fwd.prob) + log.post.prop - log.post.cur
      
      if (log(runif(1)) < acpt.rate){
        log.post.cur <- log.post.prop
        gamma <- sort(gamma.prime)
        n.acpt[1] <- n.acpt[1] + 1
        gamma.abs <- gamma.abs + 1
      }
    } else if (move.type==2){ ## remove
      flip.ix <- 'if'(gamma.abs==1, gamma, sample(gamma,1))
      gamma.prime <- gamma[gamma!=flip.ix]
      log.post.prop <- logPostc(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
      fwd.prob <- move.prob[gamma.abs+1,2]/gamma.abs
      nbhd <- c(gamma.full[!gamma.full %in% gamma],flip.ix)
      weight <- impt.prob[nbhd]/sum(impt.prob[nbhd])
      bwd.prob <- move.prob[gamma.abs,1]*weight[p-gamma.abs+1]
      acpt.rate <- log(bwd.prob) - log(fwd.prob) + log.post.prop - log.post.cur
      
      if (log(runif(1)) < acpt.rate){
        log.post.cur <- log.post.prop
        gamma <- gamma.prime
        n.acpt[2] <- n.acpt[2] + 1
        gamma.abs <- gamma.abs - 1
      }
    } else {
      flip.ix.rem <- 'if'(gamma.abs==1, gamma, sample(gamma,1))
      nbhd.fwd <- gamma.full[!gamma.full %in% gamma]
      weight.fwd <- impt.prob[nbhd.fwd]/sum(impt.prob[nbhd.fwd])
      flip.ix.add <- sample(nbhd.fwd,1,prob = weight.fwd)
      gamma.prime <- c(gamma[gamma!=flip.ix.rem],flip.ix.add)
      log.post.prop <- logPostc(gamma = gamma.prime, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
      fwd.prob <- weight.fwd[which(nbhd.fwd==flip.ix.add)]
      
      nbhd.bwd <- gamma.full[!gamma.full %in% gamma.prime]
      weight.bwd <- impt.prob[nbhd.bwd]/sum(impt.prob[nbhd.bwd])
      bwd.prob <- weight.bwd[which(nbhd.bwd==flip.ix.rem)]
      
      acpt.rate <- log(bwd.prob) - log(fwd.prob) + log.post.prop - log.post.cur
      
      if (log(runif(1)) < acpt.rate){
        log.post.cur <- log.post.prop
        gamma <- sort(gamma.prime)
        n.acpt[3] <- n.acpt[3] + 1
      }
    }
    indicator <- 'if'(iter<=burnin, iter/burnin, (iter-burnin)^(-zeta))
    for(index in gamma) impt.prob <- impt.prob + indicator*wt.update[index,]
    
    gamma.store[[iter]] <- gamma
    model.size[iter] <- gamma.abs
    log.postc.tr[iter] <- log.post.cur
  }
  gamma.mat <- t(sapply(gamma.store[-(1:burnin)],inclusion, p=p))
  incl.prob <- apply(gamma.mat,2,sum)/(n.iter-burnin)
  MPM <- which(incl.prob>=0.5)
  gamma.count <- table(apply(gamma.mat, 1, paste, collapse=" "))
  gamma.sort <- sort(gamma.count,decreasing = T)
  post.prob <- as.numeric(gamma.sort)/(n.iter-burnin)
  gamma.model <- lapply(names(gamma.sort),gammaSplit)
  HPM <- gamma.model[[1]]
  model.size.avg <- mean(model.size[-(1:burnin)])
  model.size.HPM <- length(MPM)
  model.size.MPM <- length(HPM)
  return(list(gamma.model=gamma.model, post.prob=post.prob, log.postc.tr=log.postc.tr, n.prop=n.prop, n.acpt=n.acpt, n.iter=n.iter, MPM=MPM, HPM=HPM, burnin=burnin, model.size.avg=model.size.avg,model.size.MPM=model.size.MPM,model.size.HPM=model.size.HPM))
}


pMTM <- function(X, Y, s0, init, g=nrow(X), M = 5, n.iter = 1e4, burnin = 2000, prior){
  n <- nrow(X)
  p <- ncol(X)
  
  x <- scale(X)
  y <- Y - mean(Y)
  y.norm <- sum(y^2)
  
  gamma <- 'if'(init=='null', sample.int(p-10,50,replace=F)+10, c(1:10,sample.int(p-10,40,replace=F)+10))
  n.acpt <- rep(0,3) 
  n.prop <- rep(0,3) 
  
  impt.prob <- rep(1,p)/p
  
  gamma.abs <- length(gamma)
  
  move.prob <- matrix(1/3, s0, 3)
  move.prob[1,] <- c(1,0,0)
  move.prob[s0,] <- c(0,1,0)
  
  log.postc.tr <-  model.size <- rep(NA, n.iter)
  gamma.store <- vector('list', n.iter)
  log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
  
  
  for (iter in 1:n.iter){
    move.type <-  sample(3,1,prob=move.prob[gamma.abs+1,])
    n.prop[move.type] <- n.prop[move.type] + 1
    
    if (move.type==1){ ## add
      weight <- pmin(1, impt.prob * M)
      weight[gamma] <- 0
      
      eta <-   which(as.logical(rbinom(p,1,weight)))
      n.fwd <- length(eta)
      if (n.fwd>0){
        gamma.tilde <- lapply(eta, function(j) return(c(gamma, j)))
        fwd.nbhd.ls <- sapply(gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, prior=prior, g=g)
        fwd.lp <- logSum(fwd.nbhd.ls)
        fwd.ix <- sample(n.fwd, 1, prob = exp(fwd.nbhd.ls - fwd.lp))
        gamma.prime <- gamma.tilde[[fwd.ix]]
        fwd.prob <- move.prob[gamma.abs+1,1] * weight[eta[fwd.ix]]
        
        if (gamma.abs > 0){
          bwd.gamma.tilde <- lapply(1:(gamma.abs+1), function(j) return(gamma.prime[-j]))
          bwd.nbhd.ls <- sapply(bwd.gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, prior=prior, g=g)
          bwd.lp <- logSum(bwd.nbhd.ls)
        } else {
          bwd.lp <- logPostc(gamma=gamma, y=y, x=x, y.norm=y.norm, prior=prior, g=g)
        }
        bwd.prob <- move.prob[gamma.abs+2,2]
        acpt.rate <- fwd.lp - bwd.lp + log(bwd.prob) - log(fwd.prob)
        
        if (log(runif(1))<acpt.rate) {
          gamma <- sort(gamma.prime)
          log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
          n.acpt[1] <- n.acpt[1] + 1
          gamma.abs <- gamma.abs + 1
        }
      }
    } else if (move.type==2){ ## remove
      gamma.tilde <- lapply(1:gamma.abs, function(j) return(gamma[-j]))
      fwd.nbhd.ls <- sapply(gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, prior=prior, g=g)
      fwd.lp <- logSum(fwd.nbhd.ls)
      fwd.ix <- sample(gamma.abs, 1, prob = exp(fwd.nbhd.ls - fwd.lp))
      gamma.prime <- gamma.tilde[[fwd.ix]]
      fwd.prob <- move.prob[gamma.abs+1,2]
      
      weight <- pmin(1, impt.prob * M)
      bwd.prob <- move.prob[gamma.abs,1] * weight[gamma[fwd.ix]]
      weight[gamma.prime] <- 0
      weight[gamma[fwd.ix]] <- 1
      
      eta <- which(as.logical(rbinom(p,1,weight)))
      n.bwd <- length(eta)
      bwd.gamma.tilde <- lapply(eta, function(j) return(c(gamma.prime, j)))
      bwd.nbhd.ls <- sapply(bwd.gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, prior=prior, g=g)
      bwd.lp <- logSum(bwd.nbhd.ls)
      
      
      acpt.rate <- fwd.lp - bwd.lp + log(bwd.prob) - log(fwd.prob)
      if (log(runif(1))<acpt.rate) {
        gamma <- sort(gamma.prime)
        log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
        n.acpt[2] <- n.acpt[2] + 1
        gamma.abs <- gamma.abs - 1
      }
    } else { ## swap
      weight <- pmin(1, impt.prob * M / gamma.abs)
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
        fwd.nbhd.ls <- sapply(gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, prior=prior, g=g)
        fwd.lp <- logSum(fwd.nbhd.ls)
        fwd.ix <- sample(n.fwd, 1, prob = exp(fwd.nbhd.ls - fwd.lp))
        gamma.prime <- gamma.tilde[[fwd.ix]]
        var.add <- fwd.var.add[fwd.ix]
        var.rem <- gamma[fwd.ix.rem[fwd.ix]]
        fwd.prob <- weight[var.add]
        
        weight <- pmin(1, impt.prob * M / gamma.abs)
        weight[gamma.prime] <- 0
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
        bwd.prob <- weight[var.rem]
        weight[var.rem] <- 1
        eta.add <- which(as.logical(rbinom(p,1,weight)))
        n.bwd.temp <- length(eta.add)
        bwd.var.add <- c(bwd.var.add, eta.add)
        bwd.ix.rem <- c(bwd.ix.rem, rep(gamma.abs, n.bwd.temp))
        n.bwd <- n.bwd + n.bwd.temp
        
        bwd.gamma.tilde <- lapply(1:n.bwd, function(j) return(c(gamma.prime[-bwd.ix.rem[j]], bwd.var.add[j])))
        bwd.nbhd.ls <- sapply(bwd.gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, prior=prior, g=g)
        bwd.lp <- logSum(bwd.nbhd.ls)
        
        acpt.rate <- fwd.lp - bwd.lp + log(bwd.prob) - log(fwd.prob)
        if (log(runif(1))<acpt.rate) {
          gamma <- sort(gamma.prime)
          log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
          n.acpt[3] <- n.acpt[3] + 1
        }
      }
    }
    gamma.store[[iter]] <- gamma
    model.size[iter] <- length(gamma)
    log.postc.tr[iter] <- log.post.cur
  }
  gamma.mat <- t(sapply(gamma.store[-(1:burnin)],inclusion, p=p))
  incl.prob <- apply(gamma.mat,2,sum)/(n.iter-burnin)
  MPM <- which(incl.prob>=0.5)
  gamma.count <- table(apply(gamma.mat, 1, paste, collapse=" "))
  gamma.sort <- sort(gamma.count,decreasing = T)
  post.prob <- as.numeric(gamma.sort)/(n.iter-burnin)
  gamma.model <- lapply(names(gamma.sort),gammaSplit)
  HPM <- gamma.model[[1]]
  model.size.avg <- mean(model.size[-(1:burnin)])
  model.size.HPM <- length(MPM)
  model.size.MPM <- length(HPM)
  return(list(gamma.model=gamma.model, post.prob=post.prob, log.postc.tr=log.postc.tr, n.prop=n.prop, n.acpt=n.acpt, n.iter=n.iter, MPM=MPM, HPM=HPM, burnin=burnin, model.size.avg=model.size.avg,model.size.MPM=model.size.MPM,model.size.HPM=model.size.HPM))
}


pMTMadpt <- function(X, Y, s0, init, zeta = 2/3, g=nrow(X), cor.bound = 0.75, M = 5, n.iter = 1e4, burnin=2000, prior){
  n <- nrow(X)
  p <- ncol(X)
  
  x <- scale(X)
  y <- Y - mean(Y)
  y.norm <- sum(y^2)
  
  
  x.cor <- cor(X)
  cor.bound <- min(1, max(cor.bound, 1-1/p)) ## p>=2
  cor.cut <- as.numeric(quantile(abs(x.cor[upper.tri(x.cor)]), cor.bound))
  
  nhbr <- lapply(1:p, function(i) which(abs(x.cor[,i]) > cor.cut))
  nhbr.val <- lapply(1:p, function(i) return(fnB(abs(x.cor[nhbr[[i]],i]), cor.cut)))
  
  wt.update <- sapply(1:p, function(i) wtUpdate(p, nhbr[[i]], nhbr.val[[i]]))
  
  gamma.full <- 1:p
  gamma <- 'if'(init=='null', sample.int(p-10,50,replace=F)+10, c(1:10,sample.int(p-10,40,replace=F)+10))
  
  n.acpt <- rep(0,3) 
  n.prop <- rep(0,3) 
  
  impt.probc <- rep(1,p)
  impt.prob <- impt.probc/sum(impt.probc)
  
  gamma.abs <- length(gamma)
  
  move.prob <- matrix(1/3, s0, 3)
  move.prob[1,] <- c(1,0,0)
  move.prob[s0,] <- c(0,1,0)
  
  log.postc.tr <-  model.size <- rep(NA, n.iter)
  gamma.store <- vector('list', n.iter)
  log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
  
  for (iter in 1:n.iter){
    move.type <-  sample(3,1,prob=move.prob[gamma.abs+1,])
    n.prop[move.type] <- n.prop[move.type] + 1
    
    if (move.type==1){ ## add
      weight <- pmin(1, impt.prob * M)
      weight[gamma] <- 0
      
      eta <-   which(as.logical(rbinom(p,1,weight)))
      n.fwd <- length(eta)
      if (n.fwd>0){
        gamma.tilde <- lapply(eta, function(j) return(c(gamma, j)))
        fwd.nbhd.ls <- sapply(gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, g=g, prior=prior)
        fwd.lp <- logSum(fwd.nbhd.ls)
        fwd.ix <- sample(n.fwd, 1, prob = exp(fwd.nbhd.ls - fwd.lp))
        gamma.prime <- gamma.tilde[[fwd.ix]]
        fwd.prob <- move.prob[gamma.abs+1,1] * weight[eta[fwd.ix]]
        
        if (gamma.abs > 0){
          bwd.gamma.tilde <- lapply(1:(gamma.abs+1), function(j) return(gamma.prime[-j]))
          bwd.nbhd.ls <- sapply(bwd.gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, g=g, prior=prior)
          bwd.lp <- logSum(bwd.nbhd.ls)
        } else {
          bwd.lp <- logPostc(gamma=gamma, y=y, x=x, y.norm=y.norm, g=g, prior=prior)
        }
        bwd.prob <- move.prob[gamma.abs+2,2]
        acpt.rate <- fwd.lp - bwd.lp + log(bwd.prob) - log(fwd.prob)
        
        if (log(runif(1))<acpt.rate) {
          gamma <- sort(gamma.prime)
          log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
          n.acpt[1] <- n.acpt[1] + 1
          gamma.abs <- gamma.abs + 1
        }
      }
    } else if (move.type==2){ ## remove
      gamma.tilde <- lapply(1:gamma.abs, function(j) return(gamma[-j]))
      fwd.nbhd.ls <- sapply(gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, g=g, prior=prior)
      fwd.lp <- logSum(fwd.nbhd.ls)
      fwd.ix <- sample(gamma.abs, 1, prob = exp(fwd.nbhd.ls - fwd.lp))
      gamma.prime <- gamma.tilde[[fwd.ix]]
      fwd.prob <- move.prob[gamma.abs+1,2]
      
      weight <- pmin(1, impt.prob * M)
      bwd.prob <- move.prob[gamma.abs,1] * weight[gamma[fwd.ix]]
      weight[gamma.prime] <- 0
      weight[gamma[fwd.ix]] <- 1
      
      eta <- which(as.logical(rbinom(p,1,weight)))
      n.bwd <- length(eta)
      bwd.gamma.tilde <- lapply(eta, function(j) return(c(gamma.prime, j)))
      bwd.nbhd.ls <- sapply(bwd.gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, g=g, prior=prior)
      bwd.lp <- logSum(bwd.nbhd.ls)
      
      
      acpt.rate <- fwd.lp - bwd.lp + log(bwd.prob) - log(fwd.prob)
      if (log(runif(1))<acpt.rate) {
        gamma <- sort(gamma.prime)
        log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
        n.acpt[2] <- n.acpt[2] + 1
        gamma.abs <- gamma.abs - 1
      }
    } else { ## swap
      weight <- pmin(1, impt.prob * M / gamma.abs)
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
        fwd.nbhd.ls <- sapply(gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, g=g, prior=prior)
        fwd.lp <- logSum(fwd.nbhd.ls)
        fwd.ix <- sample(n.fwd, 1, prob = exp(fwd.nbhd.ls - fwd.lp))
        gamma.prime <- gamma.tilde[[fwd.ix]]
        var.add <- fwd.var.add[fwd.ix]
        var.rem <- gamma[fwd.ix.rem[fwd.ix]]
        fwd.prob <- weight[var.add]
        
        weight <- pmin(1, impt.prob * M / gamma.abs)
        weight[gamma.prime] <- 0
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
        bwd.prob <- weight[var.rem]
        weight[var.rem] <- 1
        eta.add <- which(as.logical(rbinom(p,1,weight)))
        n.bwd.temp <- length(eta.add)
        bwd.var.add <- c(bwd.var.add, eta.add)
        bwd.ix.rem <- c(bwd.ix.rem, rep(gamma.abs, n.bwd.temp))
        n.bwd <- n.bwd + n.bwd.temp
        
        bwd.gamma.tilde <- lapply(1:n.bwd, function(j) return(c(gamma.prime[-bwd.ix.rem[j]], bwd.var.add[j])))
        bwd.nbhd.ls <- sapply(bwd.gamma.tilde, logPostc, y=y, x=x, y.norm=y.norm, g=g, prior=prior)
        bwd.lp <- logSum(bwd.nbhd.ls)
        
        acpt.rate <- fwd.lp - bwd.lp + log(bwd.prob) - log(fwd.prob)
        if (log(runif(1))<acpt.rate) {
          gamma <- sort(gamma.prime)
          log.post.cur <- logPostc(gamma = gamma, y = y, x = x, y.norm = y.norm, prior=prior, g=g)
          n.acpt[3] <- n.acpt[3] + 1
        }
      }
    }
    indicator <- 'if'(iter<=burnin, iter/burnin, (iter-burnin)^(-zeta))
    for(index in gamma) impt.probc <- impt.probc + indicator*wt.update[index,]
    impt.prob <- impt.probc/(M*impt.probc+p) 
    
    log.postc.tr[iter] <- log.post.cur
    gamma.store[[iter]] <- gamma
    model.size[iter] <- length(gamma)
  }
  gamma.mat <- t(sapply(gamma.store[-(1:burnin)],inclusion, p=p))
  incl.prob <- apply(gamma.mat,2,sum)/(n.iter-burnin)
  MPM <- which(incl.prob>=0.5)
  gamma.count <- table(apply(gamma.mat, 1, paste, collapse=" "))
  gamma.sort <- sort(gamma.count,decreasing = T)
  post.prob <- as.numeric(gamma.sort)/(n.iter-burnin)
  gamma.model <- lapply(names(gamma.sort),gammaSplit)
  HPM <- gamma.model[[1]]
  model.size.avg <- mean(model.size[-(1:burnin)])
  model.size.HPM <- length(MPM)
  model.size.MPM <- length(HPM)
  return(list(gamma.model=gamma.model, post.prob=post.prob, log.postc.tr=log.postc.tr, n.prop=n.prop, n.acpt=n.acpt, n.iter=n.iter, MPM=MPM, HPM=HPM, burnin=burnin, gamma.store=gamma.store, model.size.avg=model.size.avg,model.size.MPM=model.size.MPM,model.size.HPM=model.size.HPM, incl.prob=incl.prob))}