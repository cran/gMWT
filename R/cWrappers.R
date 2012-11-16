# Calculate the probabilistic indices for pairs and triples, using the C code
getP.Cnaive <- function(x, y , z = NULL){
  if(is.null(z))
  {
     result <- .C("getPR", as.double(x), as.double(y), as.integer(length(x)), as.integer(length(y)),result=as.double(1))$result
  } else {
     result <- .C("getPTripR", as.double(x), as.double(y), as.double(z), as.integer(length(x)), as.integer(length(y)), as.integer(length(z)),result=as.double(1))$result
  }

  result
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

permObs2.C <- function(x, y){
 nx <- length(x)
 ny <- length(y)
 res <- .C("permObs2R", as.double(x), as.double(y), as.integer(nx), as.integer(ny),result=double(nx+ny))$result
 xOut <- res[1:nx]
 yOut <- res[(nx+1):(nx+ny)]
 return(list(x=xOut,y=yOut))
}


uit.C <- function(x, y, z){
   .C("uitR", as.double(x), as.double(y), as.double(z), as.integer(length(x)), as.integer(length(y)), as.integer(length(z)),result=double(5))$result[1]
}

getRho.C <- function(x, y, z){
   .C("getRhoR", as.integer(length(x)), as.integer(length(y)), as.integer(length(z)),result=double(6))$result[1]
}
