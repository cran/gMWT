# Version: 30-11-2012, Daniel Fischer

getSigTests <- function(X, alpha=NULL, crit="distance"){

 crit <- match.arg(crit,c("distance","ratio"))

  if(dim(X)[1]<2){
   opt.alpha <- calc.alpha(X,crit=crit)
   if(is.null(alpha))alpha <- opt.alpha
   exp <- floor(alpha * ncol(X))
   repCols <- which((X[1,]<=alpha)==TRUE)
   repMat <- data.frame(p.value=X[,repCols])
   output <- list(repMat=t(repMat),obs=dim(repMat)[1],exp=exp,alpha=alpha,opt.alpha=opt.alpha,X=X)
  } else {
   repMat <- list()
   alpha.out <- c()
   opt.alpha <- c()
   obs <- c()
   exp <- c()
   for(testRun in 1:dim(X)[1])
   {
     opt.alpha[testRun] <- calc.alpha(matrix(X[testRun,],nrow=1),crit=crit)
     ifelse(is.null(alpha), alpha.out[testRun] <- opt.alpha[testRun], alpha.out[testRun] <- alpha)
     exp[testRun] <- floor(alpha.out[testRun] * ncol(X))
     repCols <- which((X[testRun,]<=alpha.out[testRun])==TRUE)
     repMat[[testRun]] <- t(data.frame(p.value=X[testRun,repCols]))
     obs[testRun] <- dim(repMat[[testRun]])[2]
   }
   output <- list(repMat=repMat,obs=obs,exp=exp,alpha=alpha.out,opt.alpha=opt.alpha,X=X,crit=crit) 
   }
   
  class(output) <- "re"
  return(output)
}

# This function calculates the optimal alpha for a given set of p-Values
calc.alpha <- function(X,crit="distance"){

  sigTests <- matrix(NA,ncol=dim(X)[2],nrow=2)
  # Expected coordinates:
  sigTests[1,] <- 1:dim(X)[2]/dim(X)[2]
  # Now through possible sigTests we can provide, based on the expected.
  # Remember, that this resolution depends on the dimension of X!!!

  for(i in 1:dim(X)[2])
  {
    temp <- (X[1,] <= sigTests[1,i])
    sigTests[2,i] <- sum(temp)

  }
  sigTests[2,] <- sigTests[2,]/dim(X)[2]
 
  
  ratios <- sigTests[2,]/sigTests[1,]
  distances <- sigTests[2,] - sigTests[1,]

  if(crit=="ratio")
  {
    opt.alpha <- min(sigTests[1,ratios==max(ratios)])
  } else if (crit=="distance"){
    opt.alpha <- min(sigTests[1,distances==max(distances)])
  } else {
    stop("Something bad happened!!!\n")
  }
  opt.alpha 
}