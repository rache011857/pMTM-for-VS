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
  n.le <<- n.le + 1
  gamma.abs <- length(gamma)
  n <- nrow(x)
  p <- ncol(x)
  if (gamma.abs>1) {
    x.gamma <- x[,gamma]
    cp <- crossprod(x.gamma,y)
    rsq.gamma <- crossprod(cp,solve(crossprod(x.gamma),cp))/y.norm
    return(-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
  } else if (gamma.abs==1){
    x.gamma <- x[,gamma]
    rsq.gamma <- sum(y*x.gamma)^2/sum(x.gamma^2)/y.norm
    return(-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
  } else{return (-n/2*log1p(g))}
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

# logPrior <- function(prior, move.type, p, gamma.abs=0, alpha=10, beta=p-10){
#   if (prior=='p'){
#     return(-move.type*log(p))
#   } else {
#     if (move.type==1){
#       return(log(gamma.abs+alpha)-log(p-gamma.abs-1+beta))
#     }
#     else {
#       return(-log(gamma.abs+alpha-1)+log(p-gamma.abs+beta))
#     }
#   }
# }


logPrior <- function(prior='bb', move.type, p, gamma.abs=0, alpha=10, beta=p-10){
  return('if'(move.type==1,log(gamma.abs+alpha)-log(p-gamma.abs-1+beta),-log(gamma.abs+alpha-1)+log(p-gamma.abs+beta)))
}

############################################################################## 

### used especially for reproducing the example 
### compute the logarithm of posterior probability of a given model (p(gamma | Y)) up to a normalizing constant

##  ------ARGUMENTS------
#    gamma       : a vector of indices of predictors
##  ------OUTPUTS------
#    the logarithm of posterior probability up to a normalizing constant

logPostc <- function(gamma){
  gamma.abs <- length(gamma)
  if (gamma.abs>1) {
    x.gamma <- x[,gamma]
    cp <- crossprod(x.gamma,y)
    rsq.gamma <- crossprod(cp,solve(crossprod(x.gamma),cp))/y.norm
    return(log(beta(alpha+gamma.abs,beta+p-gamma.abs))-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
  } else if (gamma.abs==1) {
    x.gamma <- x[,gamma]
    rsq.gamma <- sum(y*x.gamma)^2/sum(x.gamma^2)/y.norm
    return(log(beta(alpha+1,beta+p-1))-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
  } else {return (log(beta(alpha,beta+p))-n/2*log1p(g))}
}    

############################################################################## 

### convert gamma from indices form to 0/1 form
### e.g. input: gamma <- c(1,3,5,7), p <- 8 =>  output: c(1,0,1,0,1,0,1,0)

##  ------ARGUMENTS------
#    gamma  : a vector containing indices of predictors 
#    p      : number of predictors  
##  ------OUTPUTS------
#    incl   : a vector of indicators of indices   

# inclusion <- function(gamma, p) {
#   incl <- rep(0,p)
#   incl[gamma] <- 1
#   return(incl)
# }



############################################################################## 

### split and convert gamma from a factor list to index vector 

##  ------ARGUMENTS------
#    model.sort      : a index vector of factor type  
##  ------OUTPUTS------
#    a vector of index

modelSplit <- function(model.sort) {return(as.integer(strsplit(as.character(model.sort), " ")[[1]]))}

