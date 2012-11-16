\name{getSigTests}
\alias{getSigTests}

\title{
 Extract Significant Test Results for a Given Alpha or for an Optimal Alpha.
}

\description{
 This function finds the p-values less than or equal to \code{alpha}. Alpha is as given or an optimal alpha based on the ratio
 between observed and expected number of rejected tests. 

}
\usage{
  getSigTests(X,alpha=NULL,crit="distance")
}

\arguments{
  \item{X}{Matrix or vector of p-values.}
  \item{alpha}{Significance level.}
  \item{crit}{Criterium for optimal alpha.}
}

\details{
 The argument \code{X} is a matrix or vector of p-values. If \code{X} is a matrix (several tests), 
 then the rows correspond to different
 tests and the columns to different variables. If no \code{alpha} is given (the default)  an optimal alpha 
 is used. The option \code{crit} specifies the optimality criterion, \code{crit="ratio"} uses
 the ratio between observed and expected number of rejected test. The optimal alpha is the significance level for
 which this ratio is the largest one. The option \code{crit="distance"} uses the distance insted of the ratio
 and the optimal alpha is the one for which the distance between these two value is largest. 

 See also the function \code{\link[gMWT]{rejectionPlot}} for more details. 

 The optimal alpha is determined for each test separately. 
}
\value{
  A list object of class 're' with the values:
  \item{repMat}{The vector of p-values less than or equal to alpha. 
    In case \code{X} is a matrix, then \code{repMat} is a list with \code{length(repMat)=nrow(X)}.}
  \item{obs}{Observed number of rejections.}
  \item{exp}{Expected number of rejections.}
  \item{alpha}{Used \code{alpha}.}
  \item{opt.alpha}{Optimal alpha (with largest ratio of observed to expected).}
  \item{X}{The input matrix or vector \code{X}.} 
}

\seealso{
\code{\link{rejectionPlot}}
}

\author{
  Daniel Fischer
}

\examples{
  X <- matrix(c(rnorm(500,2,1),rnorm(600,2,1),rnorm(400,2.2,1)),byrow=TRUE, ncol=10)
  colnames(X) <- letters[1:10]
  g <- c(rep(1,50),rep(2,60),rep(3,40))

  test <- gmw(X,g,test="kw",type="base")

  cs1 <- getSigTests(test)
  
  cs1
}

\keyword{methods}
 
