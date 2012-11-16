
# This function calculates the test statistic for the triple type of our tests
# The outcommented lines are useful for the case, that we allow a general alternative
# formulated with a vector a. But here it only consumes calculation time...
getP.Rsub <- function(x,y,z){
    ST <- createST(x,y,z)
    EQ <- createEQ(x,y,z)
    #altType <- which(alt!=0)
    #obsCounts <- numeric(length(alt))    
    # Consider the case for alternative 1<2<3   
    #if (is.element(1,altType)){
      obsCounts <- sum(getSubMatrix(ST,1,2)%*%getSubMatrix(ST,2,3)) + 0.5*sum(getSubMatrix(ST,1,2)%*%getSubMatrix(EQ,2,3)) + 0.5*sum(getSubMatrix(EQ,1,2)%*%getSubMatrix(ST,2,3)) + (1/6)*sum(getSubMatrix(EQ,1,2)%*%getSubMatrix(EQ,2,3))
    #}
    # Consider the case for alternative 1<3<2
#     if (is.element(2,altType)){
# 	obsCounts[2] <- alt[2]*sum(getSubMatrix(ST,1,3)%*%getSubMatrix(ST,3,2))
#     }
#     # Consider the case for alternative 2<1<3
#     if (is.element(3,altType)){
# 	obsCounts[3] <- alt[3]*sum(getSubMatrix(ST,2,1)%*%getSubMatrix(ST,1,3))
#     }
#     # Consider the case for alternative 2<3<1
#     if (is.element(4,altType)){
# 	obsCounts[4] <- alt[4]*sum(getSubMatrix(ST,2,3)%*%getSubMatrix(ST,3,1))
#     }
#     # Consider the case for alternative 3<1<2
#     if (is.element(5,altType)){
# 	obsCounts[5] <- alt[5]*sum(getSubMatrix(ST,3,1)%*%getSubMatrix(ST,1,2))
#     }
#     # Consider the case for alternative 3<2<1
#     if (is.element(6,altType)){
# 	obsCounts[6] <- alt[6]*sum(getSubMatrix(ST,3,2)%*%getSubMatrix(ST,2,1))
#     }
    # The obsVounts vector has now the length of possible alternative (in the T=3 case its 6) and at each
    # entry it counted the value of the triple sum over the corresponding indicator function

    # This sums up over all chosen alternatives
    #obsCountsTotal <- sum(obsCounts[altType])
    #return(obsCountsTotal/(length(x)*length(y)*length(z)))
    
    return(obsCounts/(length(x)*length(y)*length(z)))
}

# # this function calculates the permutated values
perm.triple <- function(x,y,z,nper,type){
  Nx <- length(x)
  Ny <- length(y)
  Nz <- length(z)

  result <- c()
  if(type=="submat")
  {
    for (i in 1:nper)
    {
      permvalues <- sample(c(x,y,z))
      result[i] <- getP.Rsub(permvalues[1:Nx],permvalues[(Nx+1):(Nx+Ny)],permvalues[(Nx+Ny+1):(Nx+Ny+Nz)])
    }
  } else if (type=="cnaive"){
    for (i in 1:nper)
    {
      permvalues <- sample(c(x,y,z))
      result[i] <- getP.Cnaive(permvalues[1:Nx],permvalues[(Nx+1):(Nx+Ny)],permvalues[(Nx+Ny+1):(Nx+Ny+Nz)])
    }
  } else {
    stop("We can't do that!!!")
  }
  result
}