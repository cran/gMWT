\name{rejectionPlot}
\alias{rejectionPlot}

\title{
 Creates a Rejection Plot. 
}

\description{
  This function plots the observed and expected numbers of rejected tests and their ratio
  as a function of alpha.
}
\usage{
  rejectionPlot(X,lCol="red",xlim=NULL,crit="distance")
}

\arguments{
  \item{X}{Matrix or vector of p-values, as obtained from \code{gmw} function.}
  \item{lCol}{Vector of of line colors}
  \item{xlim}{An interval for alpha values (for the plotting). Default is [0,1].}
  \item{crit}{Shall the distance or the ratio be plotted, see details.}
}

\details{
 The matrix \code{X} contains p-values as given from the \code{gmw} function. The rows correspond to different tests,
 and the columns to different variables. 

 The function plots, for each test separately, the observed and expected numbers of rejected tests and their ratio (\code{crit="ratio"})
 or distance (\code{crit="distance"})  as a function of alpha.

}

\references{
  Fischer, D., Oja, H., Sen, P.K., Schleutker, J., Wahlfors, T. (2012): Generalized Mann-Whitney Type Tests for Microarray Experiments, submitted article.\cr
 
  Fischer, D., Oja, H. (2012): A permutation type test for calculating generalized Mann-Whitney tests, manuscript.
}

\author{
  Daniel Fischer
}

\seealso{
\code{\link{gmw}}
}

\examples{
X <- matrix(c(rnorm(500,2,1),rnorm(600,2,1),rnorm(400,2.5,1)),byrow=TRUE, ncol=10)
colnames(X) <- letters[1:10]
g <- c(rep(1,50),rep(2,60),rep(3,40))


test1 <- gmw(X,g,test="mw",type="external",alternative="greater",goi=c(2,3))
test2 <- gmw(X,g,test="kw",type="external")
test3 <- gmw(X,g,test="jt",type="external",alternative="greater")

tests <- rbind(test1,test2,test3)

rejectionPlot(tests,lCol=c("green","red","blue"))
}
\keyword{methods}
 
