# Version: 30-11-2012, Daniel Fischer

kw.gmw <- function(X,g,cluster,goi,type,nper,mc,PARAMETERS,output){

 res <- list()
 diffTests <- t(as.matrix(sort(goi)))

 METHOD <- c("********* Kruskal-Wallis Test *********")
 DNAME <- PARAMETERS[[1]]
 TEST  <- PARAMETERS[[2]]
 TYPE  <- PARAMETERS[[3]]
 ALTERNATIVE <- PARAMETERS[[4]]
 STATISTIC   <- PARAMETERS[[5]]
 PVAL        <- PARAMETERS[[6]]

 dimX      <- PARAMETERS[[7]]
 XisVector <- PARAMETERS[[8]]

 #totalObs <- length(group)

## Case: X is vector
    if(XisVector){
##---------------------------------------------------------------------------------------------------------------------------------------
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, two sided, X is vector

	    for(testRun in 1:nrow(diffTests))
	    {
	      obsValue <- kwObs(X,g,goi)
	      # if a cluster object is given, then we permute cluster internal
	      if(is.null(cluster)==T)
	      {
		nullDist <- kwNullDist(X,g,goi,nper)
	      } else {
		nullDist <- kwNullDistCluster(X,g,cluster,goi,nper)
	      }
	      PVAL <- sum(nullDist>obsValue)/nper
	
	      names(PVAL) <- "p.value"
	      STATISTIC <- obsValue
	      names(STATISTIC) <- "obs.value"
	      ALTERNATIVE <- "two.sided"
	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
	      class(resTemp)<-"htest"
	      
              res[[testRun]] <- resTemp
	      names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=1,nrow=length(res))
	      colnames(resMin) <- "pValues"
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		resMin[i,1] <- res[[i]]$p.value
	      }
	      res <- resMin
	    }

	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, two sided, X is vector
	    res <- c()
            stop("We do not have a two-sided version for the Kruskal-Wallis test, sorry!!!")
	  } else if(type=="external"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: KW from base system, two sided, X is vector
	    for(testRun in 1:nrow(diffTests))
	    {
	      testResult <- kruskal.test(X[is.element(g,goi)],g[is.element(g,goi)])
	      PVAL <- testResult$p.value
	
	      names(PVAL) <- "p.value"
	      STATISTIC <- testResult$statistic
	      names(STATISTIC) <- "obs.value"
	      ALTERNATIVE <- "two.sided"
	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
	      class(resTemp)<-"htest"
	      
              res[[testRun]] <- resTemp
	      names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=1,nrow=length(res))
	      colnames(resMin) <- "pValues"
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		resMin[i,1] <- res[[i]]$p.value
	      }
	      res <- resMin
	    }

          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, two sided, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the Kruskal-Wallis test!")
	  }
## Case: X is a matrix
    } else{
##----------------------------------------------------------------------------------------------------------------------------------------
#Preparational things for the case that X is a matrix
    # First, restrict the cores to maximum of possible tests
    if(mc>detectCores()){
	mc <- detectCores()
	warning("You do not have so many cores on this machine! I automatically reduced it to your maximum number ",mc)
    }
    mc <- min(dimX[2],mc)

   if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, two sided, X is matrix
	   innerLoop2 <- function(i,testRun,nper){
             nullDist <- kwNullDist(X[,i],g,goi,nper)
	     obsValue <- kwObs(X[,i],g,goi)
             pValue <- sum(nullDist>obsValue)/nper
	     return(list(pValue=pValue,obsValue=obsValue))
            }

	    for(testRun in 1:nrow(diffTests))
	    { 
	      resTemp <- list()
	      resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoop2,testRun=testRun,nper=nper,mc.cores=mc))
	      for(i in 1:dimX[2])
	      {
		PVAL <- resInner[2*i-1]
		STATISTIC <- resInner[2*i]
		names(PVAL) <- "p.value"
		ALTERNATIVE <- "two.sided"
		names(STATISTIC) <- "obs.value"
		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
		class(resTemp[[i]])<-"htest"	    
	      }
	     res[[testRun]] <- resTemp
	     names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=dimX[2],nrow=length(res))
	      colnames(resMin) <- colnames(X)
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		for(j in 1:dimX[2])
		{
		  resMin[i,j] <- res[[i]][[j]]$p.value
		}
	      }
	      res <- resMin
	    }
	  } else if(type=="asymptotic"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: asymptotic, two sided, X is matrix
	    res <- c()
	    stop("We do not have a asymptotic version, sorry!!!")
          } else if(type=="external"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: KW from base system, two sided, X is matrix
	   innerLoop <- function(i,testRun){
	     testResult <- kruskal.test(X[is.element(g,goi),i],g[is.element(g,goi)])
             obsValue <- testResult$statistic
             pValue <- testResult$p.value
	     return(list(pValue=pValue,obsValue=obsValue))
            }

	    for(testRun in 1:nrow(diffTests))
	    { 
	      resTemp <- list()
	      resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoop,testRun=testRun,mc.cores=mc))
	      for(i in 1:dimX[2])
	      {
		PVAL <- resInner[2*i-1]
		STATISTIC <- resInner[2*i]
		names(PVAL) <- "p.value"
		ALTERNATIVE <- "two.sided"
		#DNAME <- paste("Data:",deparse(substitute(X)),", Groups:",deparse(substitute(g)),", Order: max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],")",sep="")
		names(STATISTIC) <- "obs.value"
		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
		class(resTemp[[i]])<-"htest"	    
	      }
	     res[[testRun]] <- resTemp
	     names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
	    }
	    if(output=="min")
	    {
	      resMin <- matrix(NA,ncol=dimX[2],nrow=length(res))
	      colnames(resMin) <- colnames(X)
	      rownames(resMin) <- names(res)
	      for(i in 1:length(res))
	      {
		for(j in 1:dimX[2])
		{
		  resMin[i,j] <- res[[i]][[j]]$p.value
		}
	      }
	      res <- resMin
	    }

          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, two sided, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the triple test!,O,T,M")
	  } 
    }
  res
}