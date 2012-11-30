# Version: 30-11-2012, Daniel Fischer

jt.star.gmw <- function(X,g,goi,type,nper,alternative,mc,PARAMETERS,output){

 res <- list()
 diffTests <- t(as.matrix(sort(goi)))

 METHOD <- c("********* Jonckheere-Terpstra* Test *********")
 DNAME <- PARAMETERS[[1]]
 TEST  <- PARAMETERS[[2]]
 TYPE  <- PARAMETERS[[3]]
 ALTERNATIVE <- PARAMETERS[[4]]
 STATISTIC   <- PARAMETERS[[5]]
 PVAL        <- PARAMETERS[[6]]

 dimX      <- PARAMETERS[[7]]
 XisVector <- PARAMETERS[[8]]

## Case: X is vector
    if(XisVector){
##---------------------------------------------------------------------------------------------------------------------------------------
       if(alternative=="two.sided"){
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, two sided, X is vector

	    for(testRun in 1:nrow(diffTests))
	    {
	      obsValue <- jtStar(X[is.element(g,diffTests[testRun,])],g[is.element(g,diffTests[testRun,])])
	      nullDist <- jtStarPTest(X[is.element(g,diffTests[testRun,])],g[is.element(g,diffTests[testRun,])],nper)
	      PVAL <- 2*min(sum(nullDist>=obsValue)/nper,sum(nullDist<obsValue)/nper)
	
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
            stop("We do not have a two-sided version for the triple test, sorry!!!A,2S,V")

	  } else  if(type=="external"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: JT from the clinfun package, two.sided, X is vector
# 	    for(testRun in 1:nrow(diffTests))
# 	    { # Our greater and base greater are different interpretations, remeber that!!!
# 	      testResult <- jonckheere.test(X[is.element(g,diffTests[testRun,])],g[is.element(g,diffTests[testRun,])] ,alternative="two.sided")
# 	      PVAL <- testResult$p.value
# 	
# 	      names(PVAL) <- "p.value"
# 	      STATISTIC <- testResult$statistic
# 	      names(STATISTIC) <- "obs.value"
# 	      ALTERNATIVE <- "two.sided"
# 	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
# 	      class(resTemp)<-"htest"
# 	      
#               res[[testRun]] <- resTemp
# 	      names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
# 	    }
# 	    if(output=="min")
# 	    {
# 	      resMin <- matrix(NA,ncol=1,nrow=length(res))
# 	      colnames(resMin) <- "pValues"
# 	      rownames(resMin) <- names(res)
# 	      for(i in 1:length(res))
# 	      {
# 		resMin[i,1] <- res[[i]]$p.value
# 	      }
# 	      res <- resMin
# 	    }
	    res <- c()
	    stop("We do not have this kind of type for the tjt* test!")

          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, two sided, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the triple test!,O,2S,V")
	  }
##---------------------------------------------------------------------------------------------------------------------------------------
       } else if(alternative=="greater"){
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, greater, X is vector

	    for(testRun in 1:nrow(diffTests))
	    {
	      obsValue <- jtStar(X[is.element(g,diffTests[testRun,])],g[is.element(g,diffTests[testRun,])])
	      nullDist <- jtStarPTest(X[is.element(g,diffTests[testRun,])],g[is.element(g,diffTests[testRun,])],nper)
	      PVAL <- sum(nullDist>=obsValue)/nper
	
	      names(PVAL) <- "p.value"
	      STATISTIC <- obsValue
	      names(STATISTIC) <- "obs.value"
	      ALTERNATIVE <- "greater"
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
# Case: asymptotic, greater, X is vector
	    res <- c()
            stop("We do not have a asymptotic version for the Jonckheere-Terpstra test, sorry!!!")
          
	  } else  if(type=="external"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: JT from the clinfun package, greater, X is vector
# 	    for(testRun in 1:nrow(diffTests))
# 	    { # Our greater and base greater are different interpretations, remeber that!!!
# 	      testResult <- jonckheere.test(X[is.element(g,diffTests[testRun,])],g[is.element(g,diffTests[testRun,])] ,alternative="increasing")
# 	      PVAL <- testResult$p.value
# 	
# 	      names(PVAL) <- "p.value"
# 	      STATISTIC <- testResult$statistic
# 	      names(STATISTIC) <- "obs.value"
# 	      ALTERNATIVE <- "increasing"
# 	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
# 	      class(resTemp)<-"htest"
# 	      
#               res[[testRun]] <- resTemp
# 	      names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
# 	    }
# 	    if(output=="min")
# 	    {
# 	      resMin <- matrix(NA,ncol=1,nrow=length(res))
# 	      colnames(resMin) <- "pValues"
# 	      rownames(resMin) <- names(res)
# 	      for(i in 1:length(res))
# 	      {
# 		resMin[i,1] <- res[[i]]$p.value
# 	      }
# 	      res <- resMin
# 	    }
	    res <- c()
	    stop("We do not have this kind of type for the jt* test!")
	  } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, greater, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the triple test!,O,G,V")
	  }
       } else if(alternative=="smaller"){
##---------------------------------------------------------------------------------------------------------------------------------------
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, smaller, X is vector

	    for(testRun in 1:nrow(diffTests))
	    {
	      obsValue <- jtStar(X[is.element(g,diffTests[testRun,])],g[is.element(g,diffTests[testRun,])])
	      nullDist <- jtStarPTest(X[is.element(g,diffTests[testRun,])],g[is.element(g,diffTests[testRun,])],nper)
	      PVAL <- sum(nullDist<obsValue)/nper
	
	      names(PVAL) <- "p.value"
	      STATISTIC <- obsValue
	      names(STATISTIC) <- "obs.value"
	      ALTERNATIVE <- "smaller"
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
# Case: asymptotic, smaller, X is vector
	    res <- c()
            stop("We do not have a two-sided version for the triple test, sorry!!!,A,S,V")
	  
	  } else  if(type=="external"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: JT from the clinfun package, smaller, X is vector
# 	    for(testRun in 1:nrow(diffTests))
# 	    { # Our greater and base greater are different interpretations, remeber that!!!
# 	      testResult <- jonckheere.test(X[is.element(g,diffTests[testRun,])],g[is.element(g,diffTests[testRun,])] ,alternative="decreasing")
# 	      PVAL <- testResult$p.value
# 	
# 	      names(PVAL) <- "p.value"
# 	      STATISTIC <- testResult$statistic
# 	      names(STATISTIC) <- "obs.value"
# 	      ALTERNATIVE <- "increasing"
# 	      resTemp<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
# 	      class(resTemp)<-"htest"
# 	      
#               res[[testRun]] <- resTemp
# 	      names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
# 	    }
# 	    if(output=="min")
# 	    {
# 	      resMin <- matrix(NA,ncol=1,nrow=length(res))
# 	      colnames(resMin) <- "pValues"
# 	      rownames(resMin) <- names(res)
# 	      for(i in 1:length(res))
# 	      {
# 		resMin[i,1] <- res[[i]]$p.value
# 	      }
# 	      res <- resMin
# 	    }
	    res <- c()
	    stop("We do not have this kind of type for the jt* test!")
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, one sided, X is vector
	    res <- c()
	    stop("We do not have this kind of type for the jt* test!")
	  }
       } else {
	    res <- c()
	    stop("There is no other option than small, greater or two-sided...All other")
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

    if(alternative=="two.sided"){
       	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, two sided, X is matrix
	      # Define the function, that is performed for column i (important for parallelization)
	   innerLoop <- function(i,testRun){
             nullDist <- jtStarPTest(X[is.element(g,diffTests[testRun,]),i],g[is.element(g,diffTests[testRun,])],nper) 
             obsValue <- jtStar(X[is.element(g,diffTests[testRun,]),i],g[is.element(g,diffTests[testRun,])])
             pValue <- 2*min(sum(nullDist<obsValue)/nper,sum(nullDist>=obsValue)/nper)
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
		ALTERNATIVE <- "smaller"
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
	    warning("We do not have a two-sided version for the triple test, sorry!!!,A,T,M")

	} else if(type=="external"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: base, two sided, X is matrix
# 	   innerLoop <- function(i,testRun){
# 	     testResult <- jonckheere.test(X[is.element(g,diffTests[testRun,]),i],g[is.element(g,diffTests[testRun,])] ,alternative="two.sided")
#              obsValue <- testResult$statistic
#              pValue <- testResult$p.value
# 	     return(list(pValue=pValue,obsValue=obsValue))
#             }
# 
# 	    for(testRun in 1:nrow(diffTests))
# 	    { 
# 	      resTemp <- list()
# 	      resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoop,testRun=testRun,mc.cores=mc))
# 	      for(i in 1:dimX[2])
# 	      {
# 		PVAL <- resInner[2*i-1]
# 		STATISTIC <- resInner[2*i]
# 		names(PVAL) <- "p.value"
# 		ALTERNATIVE <- "two.sided"
# 		#DNAME <- paste("Data:",deparse(substitute(X)),", Groups:",deparse(substitute(g)),", Order: max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],")",sep="")
# 		names(STATISTIC) <- "obs.value"
# 		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
# 		class(resTemp[[i]])<-"htest"	    
# 	      }
# 	     res[[testRun]] <- resTemp
# 	     names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
# 	    }
# 	    if(output=="min")
# 	    {
# 	      resMin <- matrix(NA,ncol=dimX[2],nrow=length(res))
# 	      colnames(resMin) <- colnames(X)
# 	      rownames(resMin) <- names(res)
# 	      for(i in 1:length(res))
# 	      {
# 		for(j in 1:dimX[2])
# 		{
# 		  resMin[i,j] <- res[[i]][[j]]$p.value
# 		}
# 	      }
# 	      res <- resMin
# 	    }
	    res <- c()
	    stop("We do not have this kind of type for the jt* test!")
          } else {
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other options, two sided, X is matrix
	    res <- c()
	    stop("We do not have this kind of type for the triple test!,O,T,M")
	  }
    } else if(alternative=="greater"){
     #  res <- do.call(rbind,mclapply(c(1:dimX[2]),innerLoop,mc.cores=mc))
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, greater, X is matrix
          # Define the function, that is performed for column i (important for parallelization)
	   innerLoop <- function(i,testRun){
             nullDist <- jtStarPTest(X[is.element(g,diffTests[testRun,]),i],g[is.element(g,diffTests[testRun,])],nper) 
             obsValue <- jtStar(X[is.element(g,diffTests[testRun,]),i],g[is.element(g,diffTests[testRun,])])
             pValue <- sum(nullDist>=obsValue)/nper
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
		ALTERNATIVE <- "smaller"
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
# Case: asymptotic, greater, X is matrix
	    res <- c()
            stop("We do not have a two-sided version for the triple test, sorry!!!A,2S,V")

	} else if(type=="external"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: base, greater, X is matrix
# 	   innerLoop <- function(i,testRun){
# 	     testResult <- jonckheere.test(X[is.element(g,diffTests[testRun,]),i],g[is.element(g,diffTests[testRun,])] ,alternative="increasing")
#              obsValue <- testResult$statistic
#              pValue <- testResult$p.value
# 	     return(list(pValue=pValue,obsValue=obsValue))
#             }
# 
# 	    for(testRun in 1:nrow(diffTests))
# 	    { 
# 	      resTemp <- list()
# 	      resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoop,testRun=testRun,mc.cores=mc))
# 	      for(i in 1:dimX[2])
# 	      {
# 		PVAL <- resInner[2*i-1]
# 		STATISTIC <- resInner[2*i]
# 		names(PVAL) <- "p.value"
# 		ALTERNATIVE <- "increasing"
# 		#DNAME <- paste("Data:",deparse(substitute(X)),", Groups:",deparse(substitute(g)),", Order: max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],")",sep="")
# 		names(STATISTIC) <- "obs.value"
# 		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
# 		class(resTemp[[i]])<-"htest"	    
# 	      }
# 	     res[[testRun]] <- resTemp
# 	     names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
# 	    }
# 	    if(output=="min")
# 	    {
# 	      resMin <- matrix(NA,ncol=dimX[2],nrow=length(res))
# 	      colnames(resMin) <- colnames(X)
# 	      rownames(resMin) <- names(res)
# 	      for(i in 1:length(res))
# 	      {
# 		for(j in 1:dimX[2])
# 		{
# 		  resMin[i,j] <- res[[i]][[j]]$p.value
# 		}
# 	      }
# 	      res <- resMin
# 	    }
	    res <- c()
	    stop("We do not have this kind of type for the jt* test!")
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other, greater, X is matrix

	  } else {
	    res <- c()
	    warning("We do not have this kind of type for the UIT!,O,G,M")
	  }
    } else if(alternative=="smaller"){
     #  res <- do.call(rbind,mclapply(c(1:dimX[2]),innerLoop,mc.cores=mc))
	  if(type=="permutation"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: permutation, smaller, X is matrix
	      # Define the function, that is performed for column i (important for parallelization)
	   innerLoop <- function(i,testRun){
             nullDist <- jtStarPTest(X[is.element(g,diffTests[testRun,]),i],g[is.element(g,diffTests[testRun,])],nper) 
             obsValue <- jtStar(X[is.element(g,diffTests[testRun,]),i],g[is.element(g,diffTests[testRun,])])
             pValue <- sum(nullDist<obsValue)/nper
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
		ALTERNATIVE <- "smaller"
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
# Case: asymptotic, smaller, X is matrix
	    res <- c()
	    warning("We do not have this kind of type for the triple test!,A,S,M")

	} else if(type=="external"){
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: base, smaller, X is matrix
# 	   innerLoop <- function(i,testRun){
# 	     testResult <- jonckheere.test(X[is.element(g,diffTests[testRun,]),i],g[is.element(g,diffTests[testRun,])] ,alternative="decreasing")
#              obsValue <- testResult$statistic
#              pValue <- testResult$p.value
# 	     return(list(pValue=pValue,obsValue=obsValue))
#             }
# 
# 	    for(testRun in 1:nrow(diffTests))
# 	    { 
# 	      resTemp <- list()
# 	      resInner <-  unlist(mclapply(c(1:dimX[2]),innerLoop,testRun=testRun,mc.cores=mc))
# 	      for(i in 1:dimX[2])
# 	      {
# 		PVAL <- resInner[2*i-1]
# 		STATISTIC <- resInner[2*i]
# 		names(PVAL) <- "p.value"
# 		ALTERNATIVE <- "decreasing"
# 		#DNAME <- paste("Data:",deparse(substitute(X)),", Groups:",deparse(substitute(g)),", Order: max(P",diffTests[testRun,1],diffTests[testRun,3],",P",diffTests[testRun,2],diffTests[testRun,3],")",sep="")
# 		names(STATISTIC) <- "obs.value"
# 		resTemp[[i]]<-c(list(method=METHOD,data.name=DNAME,alternative=ALTERNATIVE,statistic=STATISTIC,test=TEST,p.value=PVAL,type=TYPE))
# 		class(resTemp[[i]])<-"htest"	    
# 	      }
# 	     res[[testRun]] <- resTemp
# 	     names(res)[testRun] <- paste(diffTests[testRun,],collapse="")
# 	    }
# 	    if(output=="min")
# 	    {
# 	      resMin <- matrix(NA,ncol=dimX[2],nrow=length(res))
# 	      colnames(resMin) <- colnames(X)
# 	      rownames(resMin) <- names(res)
# 	      for(i in 1:length(res))
# 	      {
# 		for(j in 1:dimX[2])
# 		{
# 		  resMin[i,j] <- res[[i]][[j]]$p.value
# 		}
# 	      }
# 	      res <- resMin
# 	    }
	    res <- c()
	    stop("We do not have this kind of type for the jt* test!")
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other, smaller, X is matrix

          } else {
	    res <- c()
	    stop("We do not have this kind of type for the triple test!,O,S,M")
	  }
#----------------------------------------------------------------------------------------------------------------------------------------
# Case: other, other, X is matrix

    } else {
	    res <- c()
	    stop("There are no other alternatives possible, sorry! All other....")
     }
  }
  res
}