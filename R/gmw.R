gmw <- function(X,g,goi=NULL,test="mw",type="permutation",prob="pair",nper=2000,alternative="greater",mc=1,output="min",cluster=NULL,order=TRUE){

 type <- match.arg(type,c("permutation","base","asymptotic","submat"))
 test <- match.arg(test,c("uit","triple","jt","jt*","mw","kw"))
 prob <- match.arg(prob,c("single","pair","triple"))
 alternative <- match.arg(alternative,c("smaller","greater","two.sided"))
 output <- match.arg(output,c("min","full"))

 # Adjust the group labels
 g <- relabelGroups(g)

 # Get flags, if the input X will be analyzed as vector or matrix and calculate certain constants
 dimX <- dim(X)
 XisVector <- is.null(dimX)
 if(is.null(goi)) goi <- g
 goi <- unique(goi)

 # Initialze the class variables for the output
 DNAME <- paste(deparse(substitute(X)),", grouping vector:",deparse(substitute(group)))
 TEST<- switch(test,"uit"="Union Intersection Test",
                    "triple"="Triple Based Test",
		    "jt"="Jonckheere-Terpstra",
		    "jt*"="Jonckheere-Terpstra *",
		    "mw"="Mann-Whitney Test",
		    "kw"="Kruskal-Wallis Test")
 TYPE <- switch(type,"permutation"="Permutation Test",
                     "asymptotic"="Asymptotic Test",
		     "base"="Included in base/other package Test")
 ALTERNATIVE=""
 STATISTIC=""
 PVAL=""

 PARAMETERS <- list(DNAME,TEST,TYPE,ALTERNATIVE,STATISTIC,PVAL,dimX,XisVector)

 if(test=="uit"){

    res <- uit.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output)

 } else if(test=="triple"){

    res <- triple.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output)

} else if(test=="mw"){

    if(prob=="pair")
    {
      res <- mw.gmw(X=X,g=g,goi=goi,type=type,nper=nper,alternative=alternative,mc=mc,PARAMETERS,output=output,order=order)
    } else if(prob=="single"){
      res <- mw.single.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output)
    }

} else if(test=="jt"){

    res <- jt.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output)

} else if(test=="jt*"){

    res <- jt.star.gmw(X,g,goi,type,nper,alternative,mc,PARAMETERS,output)

} else if(test=="kw"){

    res <- kw.gmw(X,g,cluster,goi,type,nper,mc,PARAMETERS,output)

 } else { 

    stop("There is not such test, sorry!")
    res <- NULL

 }

 return(res)

}