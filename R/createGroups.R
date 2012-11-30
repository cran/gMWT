# Version: 30-11-2012, Daniel Fischer


# this function creates the new vector of group labels
createGroups <- function(g,desOrder){
  reG <- g
  curClass <- 1
  for(i in 1:length(desOrder))
  {
    tempPos <- g==desOrder[i]
    reG[tempPos] <- curClass
    curClass <- curClass + 1
  }
  reG
} # End of function: createGroups