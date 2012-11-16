rejectionPlot <- function(X,lCol="red",xlim=NULL,crit="distance"){

  crit <- match.arg(crit,c("distance","ratio"))

  if(is.null(xlim)) xlim <- c(0,1)
  wasMatrix <- is.matrix(X) && dim(X)[1]>1
  if(!is.matrix(X)) X <- t(as.matrix(X))

  if(length(lCol)!=dim(X)[1])
  {
    warning("Too less colors given, colors will be repeated!")
    lCol <- rep(lCol,dim(X)[1])[1:dim(X)[1]]
  }

  sigTests <- matrix(NA,ncol=dim(X)[2],nrow=dim(X)[1]+1)
  # Expected coordinates:
  sigTests[1,] <- 1:dim(X)[2]/dim(X)[2]
  # Now through possible sigTests we can provide, based on the expected.
  # Remeber, that this resolution depends on the dimension of X!!!

  for(i in 1:dim(X)[2])
  {
    temp <- (X <= sigTests[1,i])
    for(j in 2:dim(sigTests)[1])
    {
      sigTests[j,i] <- sum(temp[j-1,])
    }
  }
  for(i in 2:dim(sigTests)[1])
  {
    sigTests[i,] <- sigTests[i,]/dim(X)[2]
  }
  
  ratios <- matrix(NA,ncol=ncol(sigTests),nrow=nrow(sigTests)-1)
  distances <- matrix(NA,ncol=ncol(sigTests),nrow=nrow(sigTests)-1)
  for(i in 1:dim(ratios)[1])
  {
    ratios[i,] <- sigTests[i+1,]/sigTests[1,]
    distances[i,] <- sigTests[i+1,] - sigTests[1,]
  }
  
  ylim <- c(0,max(sigTests[,sum(sigTests[1,]<=xlim[2])]))
  
 #nf <- layout(matrix(c(1,2),ncol=1), widths=c(7), heights=c(7,2), TRUE)
 #nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
 nf <- layout(matrix(c(1,2),ncol=1),c(4,4), c(3,1), TRUE)
# layout.show(nf)

  
  par(oma=c(2,0,1,0)
 #     ,pty="s"
      ,mar=c(4,4,1,0)
     )
  plot(c(0,sigTests[1,]),c(0,sigTests[2,]),col=lCol[1],type="l",ylab="Observed Ratio",xlab=" ",xlim=xlim,ylim=ylim)
  lines(c(-10,10),c(-10,10),type="l")
  if(wasMatrix)
  {
    for(i in 3:dim(sigTests)[1])
    {
      lines(c(0,sigTests[1,]),c(0,sigTests[i,]),type="l",col=lCol[i-1])
    }
  }

 if(crit=="ratio")
 {
   plot(c(-1,2),c(1,1),type="l",xlim=xlim,ylim=c(min(ratios),max(ratios)),xlab="Expected Ratio",ylab="Ratio")
   for(i in 1:dim(ratios)[1])
   {
     lines(sigTests[1,],ratios[i,],col=lCol[i])
   }
  } else if(crit=="distance") {
     plot(c(-1,2),c(1,1),type="l",xlim=xlim,ylim=c(min(distances),max(distances)),xlab="Expected Ratio",ylab="Distance")
     for(i in 1:dim(ratios)[1])
     {
       lines(sigTests[1,],distances[i,],col=lCol[i])
     }
  } else {
    stop("We have a problem!!!\n")
  }
}