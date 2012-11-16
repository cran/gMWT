\name{gmw}
\alias{gmw}

\title{
  Generalized Mann-Whitney Type Tests
}

\description{
  The function calculates p-values for different tests as presented in the paper "Generalized Mann-Whitney Type Tests for Microarray Experiments".
}
\usage{
  gmw(X,g,goi=NULL,test="mw",type="permutation",prob="pair",nper=2000,
  alternative="greater",mc=1,output="min", cluster=NULL, order = TRUE)
}

\arguments{
  \item{X}{Data matrix, each column corresponds to a variable, each row to an individual. Can also be a vector (one variable).}
  \item{g}{Vector of length \code{nrow(X)} (respective \code{length(X)}), assigning treatment groups (numbers)
   to observations, see details.}
  \item{goi}{Vector with elements of \code{g}, defining for which treatment groups the test should be performed.}
  \item{test}{Specifies the test statistic.}
  \item{type}{Permutation test (\code{"permutation"}) or normal approximation (\code{"approximation"}) for the
   calculation of the p-values. Tests implemented
		in R-base \code{"base"} are also accessible, see details.}
  \item{prob}{This option is only for the Mann-Whitney test, see details.}
  \item{nper}{If type is \code{"permutation"} this option specifies how many permutations are used to calculate the p-value.}
  \item{alternative}{Specifies the alternative, the options are \code{"smaller"},\code{"greater"} and \code{"two.sided"}, see details.}
  \item{mc}{Multiple Cores, determines how many tests will be performed parallel (only available under Linux), see details.}
  \item{output}{Determines the level of the details in the output.}
  \item{cluster}{A vector of same length as \code{g}, giving possible cluster information for performing the permutation test.}
  \item{order}{Boolean, shall all orders be calculated or only increasing orders?}
}

\details{
  The object \code{X} is the data vector (one variable) or the data matrix.
  Each row refers to  an observation, and  each column to a variable. The tests are performed separately for all variables.

  The vector \code{g} gives the group number. The directional tests are based on this numbering of the group.


  The \code{goi} option defines, which treatment groups are used in the test constructions.
  If no groups are specified (default), all groups are used.

  The \code{test} option specifies the test statistic. Possible options are \code{'uit'}(union intersection test), \code{'triple'} (test based on triple indicator functions),
  \code{'jt'} for Jonckheere-Terpstra test, \code{'jt*'} for a modified Jonckheere-Terpstra test, 
  \code{'mw'} for the Mann-Whitney / Wilcoxon test and 
  \code{'kw'} for the Kruskal-Wallis test. See also reference [1] for further details.

  The option \code{type} is used to decide how the p-values are computed. If \code{test='mw'}, 
  \code{test='kw'} or \code{test='jt'} then possible options are 
  \code{type='base'} and \code{type='permutation'}. If  \code{test='uit'}, \code{test='triple'} or \code{test='jt*'}.
  one can choose options \code{type='permutation'} and \code{type='asymptotic'}.

  The \code{prob} option is only for the Mann-Whitney test. For the option \code{"single"}, the tests are to
  compare a single group versus all the other groups. The option \code{"pair"} makes all pair-wise comparisons
  between the groups. 

  The option \code{alternative} is used to specify whether one-sided or two-sided alternatives are used. 
  If the test is based on the PIs, the option  \code{"greater"}  for example means that, according to the alternative,
  the groups with larger group numbers tend to have larger observations as well. The function \code{createGroups}
  may be used to renumber the groups, if needed. 
 
  The \code{mc} option is only valid if \code{X} is a matrix and the used OS is Linux, because the parallelisation is based on the package parallel, and that again is
  based on the concept of forking, which is currently only supported under Linux. 

  
  The option \code{output} can be used to control how detailed the output is. The default \code{"min"}  reports just 
  the matrix of p-values in a matrix (columns=variables, rows=tests). If \code{output="full"},
  a list will be returned with items containing  full test objects of class \code{htest}.

  The option \code{cluster} is an additional object for the Kruskal-Wallis permutation test. For cluster-dependent
  observation, only the permutations within clusters are acceptable for the p-value calculation.
}
\value{
  A matrix or vector of p-values of the underlying hypothesis test(s). In case of \code{output="full"} we give a list, and each list
  item contains the \code{htest} object for the column-wise performed test.
}

\references{
  [1] Fischer, D., Oja, H., Sen, P.K., Schleutker, J., Wahlfors, T. (2012): Generalized Mann-Whitney Type Tests for Microarray Experiments, submitted article.\cr

  [2] Fischer, D., Oja, H. (2012): A permutation type test for calculating generalized Mann-Whitney tests, manuscript.
}

\author{
  Daniel Fischer
}

\examples{
  X <- c(sample(15))
  X <- c(X,101,102,103)
  g <- c(1,1,1,2,2,2,2,3,3,3,4,4,4,4,4,5,5,5)
  cluster=c(rep(c(1,2),9))

  gmw(X,g,test="kw",type="base")
  gmw(X,g,test="kw",type="permutation")
  gmw(X,g,test="kw",type="permutation",cluster=cluster)

  gmw(X,g,test="jt",type="permutation")

}
\keyword{
  methods
}
 
