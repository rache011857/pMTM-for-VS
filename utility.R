############################################################################## 

### compute sum of numbers in log scale 

##  ------ARGUMENTS------
#    lx    : a vector containing values in log scale
##  ------OUTPUTS------
#    sum of numbers in lx

logSum <- function(lx) return(max(lx) + log(sum(exp(lx - max(lx)))))

############################################################################## 

### score function based on correlations between predictors

##  ------ARGUMENTS------
#    a    : correlation between two predictors
#    b    : a threshold of correlations
##  ------OUTPUTS------
#    score 

fnB <- function(a, b) {
  vv <- exp(-9 * (1-a)/(1-b))
  return(vv)
}

############################################################################## 

### add score to predictors 

##  ------ARGUMENTS------
#    p          : number of predictors
#    ix.vec     : a vector of indices which score need to be updated
#    wt.vector  : score vector  
##  ------OUTPUTS------
#    a score vector

wtUpdate <- function(p, ix.vec, wt.vector){
  score <- rep(0,p)
  score[ix.vec] <- wt.vector
  return(score)
}




############################################################################## 

### compute the logarithm of marginal likelihood (log L(Y | gamma)) up to a normalizing constant

##  ------ARGUMENTS------
#    x           : design matrix (after normalization)
#    y           : response (after centering)
#    g           : hyperparameter used in Zellner's g-prior, default is sample size (Unit information prior)
#    gamma       : a vector of indices of predictors
#    y.norm      : l_2 norm of response
##  ------OUTPUTS------
#    logarithm of marginal likelihood up to a normalizing constant
logMl <- function(gamma, y, x, y.norm, g=nrow(x)){
  gamma.abs <- length(gamma)
  n <- nrow(x)
  p <- ncol(x)
  if (gamma.abs==0) {return (-n/2*log1p(g))
  } else if (gamma.abs==1){
    x.gamma <- x[,gamma]
    rsq.gamma <- sum(y*x.gamma)^2/sum(x.gamma^2)/y.norm
    return(-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
  } else{
    x.gamma <- x[,gamma]
    cp <- crossprod(x.gamma,y)
    rsq.gamma <- crossprod(cp,solve(crossprod(x.gamma),cp))/y.norm
    return(-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
  }
}

############################################################################## 

### compute the difference of logarithm of prior during MH step for add/remove move

##  ------ARGUMENTS------
#    prior       : if 'p': model prior is (1/p)^modelsize; else: Beta-Binomial prior
#    move.type   : -1: remove a predictor; 1: add a predictor
#    p           : number of predictors
#    gamma.abs   : current model size
#    alpha       : hyperparameter of Beta-Binomial prior, default is 10
#    beta        : hyperparameter of Beta-Binomial prior, default is p-10
##  ------OUTPUTS------
#    the difference of logarithm of prior using in the MH updating 

logPrior <- function(prior, move.type, p, gamma.abs=0, alpha=10, beta=p-10){
  if (prior=='p'){
    return(-move.type*log(p))
  } else {
    if (move.type==1){
      return(log(gamma.abs+alpha)-log(p-gamma.abs-1+beta))
    }
    else {
      return(-log(gamma.abs+alpha-1)+log(p-gamma.abs+beta))
    }
  }
}

############################################################################## 

### used especially for reproducing the example 
### compute the logarithm of posterior probability of a given model (p(gamma | Y)) up to a normalizing constant

##  ------ARGUMENTS------
#    prior       : 'p': model prior is (1/p)^kappa*modelsize; 'uni': Uniform prior; else: Beta-Binomial prior
#    x           : design matrix (after normalization)
#    y           : response (after centering)
#    y.norm      : l_2 norm of response
#    gamma       : a vector of indices of predictors
#    g           : hyperparameter used in Zellner's g-prior, default is sample size (Unit information prior)
#    kappa       : hyperparameter used in sparisty prior, default is 3 (for reproducing)
#    alpha       : hyperparameter of Beta-Binomial prior, default is 10
#    beta        : hyperparameter of Beta-Binomial prior, default is p-10
##  ------OUTPUTS------
#    the logarithm of posterior probability up to a normalizing constant

logPostc <- function(gamma, y, x, y.norm, prior, g=nrow(x), kappa=3, alpha=10, beta=ncol(x)-10){
  gamma.abs <- length(gamma)
  n <- nrow(x)
  p <- ncol(x)
  if (prior=='p'){
    if (gamma.abs==0) {return (-n/2*log1p(g))
    } else if (gamma.abs==1){
      x.gamma <- x[,gamma]
      rsq.gamma <- sum(y*x.gamma)^2/sum(x.gamma^2)/y.norm
      return(-kappa*gamma.abs*log(p)-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
    } else{
      x.gamma <- x[,gamma]
      cp <- crossprod(x.gamma,y)
      rsq.gamma <- crossprod(cp,solve(crossprod(x.gamma),cp))/y.norm
      return(-kappa*gamma.abs*log(p)-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
    }
  } else if (prior=='uni') {
    if (gamma.abs==0) {return (-n/2*log1p(g))
    } else if (gamma.abs==1){
      x.gamma <- x[,gamma]
      rsq.gamma <- sum(y*x.gamma)^2/sum(x.gamma^2)/y.norm
      return(-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
    } else{
      x.gamma <- x[,gamma]
      cp <- crossprod(x.gamma,y)
      rsq.gamma <- crossprod(cp,solve(crossprod(x.gamma),cp))/y.norm
      return(-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
    }
  } else {
    if (gamma.abs==0) {return (log(beta(alpha,beta+p))-n/2*log1p(g))
    } else if (gamma.abs==1){
      x.gamma <- x[,gamma]
      rsq.gamma <- sum(y*x.gamma)^2/sum(x.gamma^2)/y.norm
      return(log(beta(alpha+1,beta+p-1))-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
    } else{
      x.gamma <- x[,gamma]
      cp <- crossprod(x.gamma,y)
      rsq.gamma <- crossprod(cp,solve(crossprod(x.gamma),cp))/y.norm
      return(log(beta(alpha+gamma.abs,beta+p-gamma.abs))-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
    }
  }
}

############################################################################## 

### convert gamma from indices form to 0/1 form
### e.g. input: gamma <- c(1,3,5,7), p <- 8 =>  output: c(1,0,1,0,1,0,1,0)

##  ------ARGUMENTS------
#    gamma  : a vector containing indices of predictors 
#    p      : number of predictors  
##  ------OUTPUTS------
#    incl   : a vector of indicators of indices   

inclusion <- function(gamma, p) {
  incl <- rep(0,p)
  incl[gamma] <- 1
  return(incl)
}



############################################################################## 

### split and convert gamma from character to 0/1 numerical values 

##  ------ARGUMENTS------
#    gamma.sort      : a vector of appearance times with model as vector name  
##  ------OUTPUTS------
#    an indicator vector

gammaSplit <- function(gamma.sort) {return(which(as.numeric(strsplit(gamma.sort, " ")[[1]])==1))}

############################################################################## 

### compute the beta given a single model 

##  ------ARGUMENTS------
#    x           : design matrix (after normalization)
#    Y           : response
#    sdx         : a vector containing standard deviations of each predictor
#    g           : hyperparameter used in Zellner's g-prior
#    gamma       : a vector of indices of predictors
##  ------OUTPUTS------
#    beta corresponding to the given model

modelCoef <- function(gamma, g, x, Y, sdx){
  gamma.abs <- length(gamma)
  p <- ncol(x)
  coef <- rep(0,p)
  if (gamma.abs==0) {return(coef)} else if (gamma.abs==1){
    coef[gamma] <- g/(g+1)*sum(x[,gamma]*Y)/sum(x[,gamma]^2)/sdx[gamma]
    return(coef)
  } else {
    coef[gamma] <- g/(g+1)*solve(x[,gamma],crossprod(x[,gamma],Y))/sdx[gamma]
    return(coef)
  }
}

############################################################################## 

### compute the beta using Bayesian model averaging

##  ------ARGUMENTS------
#    x           : design matrix (after normalization)
#    Y           : response
#    sdx         : a vector containing standard deviations of each predictor
#    g           : hyperparameter used in Zellner's g-prior
#    gamma.model : a list, every element is a vector of indices of predictors
#    post.prob   : posterior probabilities of models in gamma.model
##  ------OUTPUTS------
#    the beta using Bayesian model averaging

bmaCoef <- function(x, Y, sdx, g, gamma.model, post.prob){
  coef <- t(sapply(gamma.model, modelCoef, g=g, x=x, Y=Y, sdx=sdx))
  return(apply(coef*post.prob,2,sum))
}

############################################################################## 

### compute RSS of given design matrix, response and beta  

##  ------ARGUMENTS------
#    beta      : a coefficient vector 
#    x         : design matrix (after normalization)
#    y         : response (after centering)
##  ------OUTPUTS------
#    RSS

score <- function(beta,x,y){
  return(sum((y-x%*%beta)^2))
}

############################################################################## 

### get the best beta from EMVS output in terms of RSS (equivalent to BIC)  
### from the betas corresponding to highest probability model 

##  ------ARGUMENTS------
#    emvs.fit  : a list containing the result returned from EMVS()
#    x         : design matrix (after normalization)
#    y         : response (after centering)    
##  ------OUTPUTS------
#    the beta achieves minimal RSS    

emvsBest <- function(emvs.fit,x,y){
  log.post <- emvs.fit$log_g_function
  max.ix <- which(log.post==max(log.post))
  if (length(max.ix)==1){
    return(emvs.fit$betas[max.ix,])
  } else {
    best.ix <- which.min(apply(emvs.fit$betas[max.ix,], 1, score, x=x, y=y))
    return(emvs.fit$betas[max.ix[best.ix],])
  }
}