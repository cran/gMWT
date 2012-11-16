\name{gMWT-package}
\alias{gMWT-package}
\alias{gMWT}
\docType{package}

\title{
  Implementations of Generalized Mann-Whitney Type Tests
}

\description{
The package provides nonparametric tools for the comparison of several groups/treatments when 
the number of variables is large. The tools are the following.

(i) Computation of the so called probabilistic index (PI) estimates for the group comparisons.
The probabilistic indices here are (a) the probability \eqn{P_t} that a random observation from group \eqn{t}
is smaller that a random observation from other groups, (b) the probability \eqn{P_{tt'}} that 
a random observation from group \eqn{t} is smaller than a random observation from group \eqn{t'},
and (c) the probability \eqn{P_{tt't''}} that random observations from groups \eqn{t,t',t''} 
appear in this same order. The tools are also given to produce the plots of PIs for set of variables.

(ii) Computation of the p-values of some classical  and some new nonparametric tests for the comparison of
several groups/treatments.  The tests are based on the use of the  probabilistic 
indices (PI) \eqn{P_t}, \eqn{P_{tt'}} and \eqn{P_{tt't''}}. Classical Wilcoxon test, Mann-Whitney test and 
Jonkheere-Terpstra tests are included.

(iii) Tools for the simultaneous testing problem. As the package is meant for the analysis of gene expression data, 
for example, the plots of the numbers of observed and expected rejected null hypotheses and their ratio as
a function of alpha are provided. 
}



\details{
  \tabular{ll}{
    Package: \tab GMW\cr
    Type: \tab Package\cr
    Version: \tab 0.1\cr
    Date: \tab 2012-11-16\cr
    License: \tab GPL\cr
    LazyLoad: \tab yes\cr
  }
}

\author{
  Daniel Fischer, Hannu Oja

  Maintainer: Daniel Fischer <daniel.fischer@uta.fi>
}

\references{
  Fischer, D., Oja, H., Sen, P.K., Schleutker, J., Wahlfors, T. (2012): Generalized Mann-Whitney Type Tests for Microarray Experiments, submitted article.

  Fischer, D., Oja, H. (2012): A permutation type test for calculating generalized Mann-Whitney tests, manuscript.
}

\keyword{package}