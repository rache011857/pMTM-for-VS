gamma.store <-  fit$gamma.store
n.le <- fit$n.lhe

####################################################################
## posterior summary 
# method 1
model.count <- count(sapply(gamma.store, paste, collapse=" "))
model.store <- as.data.frame(matrix(NA, ncol = 4, nrow = nrow(model.count)))
model.store[,1] <- as.character(model.count$x)
model.store[,2] <- model.count$freq
model <- lapply(model.store[,1],modelSplit)
firstAppear <- function(gamma) {return(n.le[which.max(lapply(gamma.store,identical,gamma))])}

system.time(model.store[,3] <- simplify2array(mclapply(model, firstAppear, mc.cores=4)))



# method 2
model.store <- data.frame(gamma=character(), freq=integer(), first=integer(), post=double(), stringsAsFactors = F)
buildModelDf <- function(iter) {
  char <- paste(gamma.store[[iter]],collapse = " ")
  index <- which(model.store$gamma==char)
  if (length(index)) {model.store[index,2] <<- model.store[index,2] + 1} else {
    len <- nrow(model.store)
    model.store[len+1,1] <<- char
    model.store[len+1,2] <<- 1
    model.store[len+1,3] <<- n.le[iter]
  }
}

system.time(for (i in 1:1e5) {buildModelDf(i)})





# 
# alpha <- 10
# beta <- p-10
# x <- scale(X)
# y <- Y - mean(Y)
# y.norm <- sum(y^2)
# g <- n
# 
# logPostc <- function(gamma){
#   gamma <- modelSplit(gamma)
#   gamma.abs <- length(gamma)
#   if (gamma.abs==0) {return (log(beta(alpha,beta+p))-n/2*log1p(g))
#     } else if (gamma.abs==1){
#       x.gamma <- x[,gamma]
#       rsq.gamma <- sum(y*x.gamma)^2/sum(x.gamma^2)/y.norm
#       return(log(beta(alpha+1,beta+p-1))-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
#     } else{
#       x.gamma <- x[,gamma]
#       cp <- crossprod(x.gamma,y)
#       rsq.gamma <- crossprod(cp,solve(crossprod(x.gamma),cp))/y.norm
#       return(log(beta(alpha+gamma.abs,beta+p-gamma.abs))-gamma.abs/2*log1p(g)-n/2*log1p(g*(1-rsq.gamma)))
#     }
# }
# 
# model.sort <- model.store[order(-model.store$freq),][1:1000,]
# model.sort$post <- sapply(model.sort$gamma,logPostc)
# model.sort <- model.sort[order(model.sort$first),]
# 
# 
# xaxis <- model.sort$first
# yaxis <- cumsum(exp(model.sort$post))
# plot(xaxis,yaxis)
# 
# model.index <- function(gamma) {
#   incl <- rep(0,p)
#   incl[gamma] <- 1
#   return(strtoi(paste(incl,collapse=""),base=2))
# }
# 
# model.index(gamma.store[[1]])
