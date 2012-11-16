# This function calculates the different probabilistic indices for P_t , P_tt' and P_tt't''
# It supports (under Linux) parallel computing and offers several different algorithms to
# calculate the probability estimators, which should in general provide the same results, but
# have of course different calculation times. Different algorithms are mainly here to validate
# the different approaches as well as compare different implementation with respect to calculation
# time.

estPI <- function(X,g,type="pair",goi=NULL,mc=1,order=TRUE,alg="Cnaive"){

 type <- match.arg(type,c("single","pair","triple"))
 alg <- match.arg(alg,c("Cnaive","Rsub","Rnaive","Rgrid"))
 g <- relabelGroups(g)

 # Define the output object
 result <- list()

 # Get flags, if the input X will be analyzed as vector or matrix (matrix input can be handled parallel)
 dimX <- dim(X)
 XisVector <- is.null(dimX)
 if(is.null(goi)) goi <- g
 goi <- unique(goi)
 Ngoi <- length(goi)

 # First remove those entries, which are not requested via goi:
 pickThose <- is.element(g,goi)
 ifelse(XisVector , X <- X[pickThose] , X <- X[pickThose,])
 g <- g[pickThose]

 # Input check for amount of cores
 if(mc>detectCores())
 {
   mc <- detectCores()
   warning(paste("You do not have so many cores on this machine! I automatically reduced it to your maschines maximum number:",mc,"\n"))
 }
 mc <- min(dimX[2],mc)

 # Check first, how many probabilistic indices are requested (depending on input type, goi and order)
 NlistItems <- c()
 if(type=="single")
 {
  # In case of P_t possibe PE are for each single group
  NlistItems <- Ngoi
  comb <- goi
 } else if (type=="pair"){
  NlistItems <- getNListItems(Ngoi,"pair",order)
  comb <- getComb(goi,type="pairs",order)
 } else if (type=="triple"){
  NlistItems <- getNListItems(Ngoi,"triple",order)
  comb <- getComb(goi,type="triple",order)
 }

# Loops for parallel calcuations in the 3 cases (only available when X is a matrix):
# Remember, those functions are ONLY to apply within the following FOR-loops and hence
# they contain variables, which are only visible within this environment!!!
inner1 <- function(i,X,g,t,goi,alg){
  PE1(X[,i],g,t,goi,alg) 
}
inner2 <- function(i,X,g,comb,alg){
  PE2(X[,i],g,comb,alg) 
}
inner3 <- function(i,X,g,comb,alg){
  PE3(X[,i],g,comb,alg) 
}

#inner3 <- function(i,alg="C"){
#   if(alg=="C")
#   {
#     getP.C(X[g==comb[oRun,1],i],X[g==comb[oRun,2],i],X[g==comb[oRun,3],i])
#   } else if(alg=="Rsub") {
#     triple.subMatEstimate(X[g==comb[oRun,1],i],X[g==comb[oRun,2],i],X[g==comb[oRun,3],i],alt=c(1,0,0,0,0,0))
#   } else {
#     stop("We do not have this kind of PI for triples, sorry!")
#   }
#}


 # Now go through all the possible options for the PI
 for(oRun in 1:NlistItems)
 {  
    if(XisVector)
    {
      # Calculate here, if just one variable is given
      if(type=="single")
      {
	result[[oRun]] <- PE1(X,g,goi[oRun],goi,alg) 
	names(result)[oRun] <- paste("P(",goi[oRun],")",sep="")

      } else if(type=="pair")
	{
	result[[oRun]] <- PE2(X,g,comb[oRun,],alg)
	names(result)[oRun] <- paste("P(",paste(comb[oRun,],collapse="<"),")",sep="")

      } else if(type=="triple"){
	  

#           if(alg=="C")
# 	  {
# 	    result[[oRun]] <- getP.C(X[g==comb[oRun,1]],X[g==comb[oRun,2]],X[g==comb[oRun,3]])
# 	  } else if(alg=="sub") {
#  	    result[[oRun]] <-  triple.subMatEstimate(X[g==comb[oRun,1]],X[g==comb[oRun,2]],X[g==comb[oRun,3]])
# 	  }
	result[[oRun]] <- PE3(X,g,comb[oRun,],alg)
	names(result)[oRun] <- paste("P(",paste(comb[oRun,],collapse="<"),")",sep="")
      }
  
    } else {
      # Calculate PE, if more than one variable is given
      if(type=="single")
      {
	result[[oRun]] <- unlist(mclapply(c(1:dimX[2]),inner1,X,g,goi[oRun],goi,alg,mc.cores=mc))
	names(result)[oRun] <- paste("P(",goi[oRun],")",sep="")

      } else if(type=="pair"){
	result[[oRun]] <- unlist(mclapply(c(1:dimX[2]),inner2,X,g,comb[oRun,],alg,mc.cores=mc))
	names(result)[oRun] <- paste("P(",paste(comb[oRun,],collapse="<"),")",sep="")

      } else if(type=="triple"){
	result[[oRun]] <- unlist(mclapply(c(1:dimX[2]),inner3,X,g,comb[oRun,],alg,mc.cores=mc))
	#result[[oRun]] <- unlist(mclapply(c(1:dimX[2]),inner3,alg,mc.cores=mc))
	names(result)[oRun] <- paste("P(",paste(comb[oRun,],collapse="<"),")",sep="")
      }
    }
 }

# Rearrange the output, depending on Vector/Matrix input
 if(XisVector)
 {
  result <- unlist(result) 
 } else {
  temp <- matrix(NA,nrow=NlistItems,ncol=dimX[2])
  for(i in 1:NlistItems)
  {
    temp[i,] <- result[[i]]
  }
  rownames(temp) <- names(result)
  colnames(temp) <- colnames(X)
  result <- temp
 }

 res <- list(probs=result,type=type,goi=goi,order=order,obs=X,g=g,alg=alg)

 class(res) <- "estPI"
 return(res)
} # End of function: estPI
