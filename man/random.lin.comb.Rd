\name{random.lin.comb}
\alias{random.lin.comb}
\alias{random.lin.comb.default}
\alias{random.lin.comb.lm}
\title{Compute a random linear combination of the columns of a matrix or data frame }
\description{
Computes \eqn{Xa} where \eqn{X} is an \eqn{n \times p}{n by p} matrix, and \eqn{a} is
a random vector of length \eqn{p}.
}
\usage{
random.lin.comb(X, seed = NULL)

\method{random.lin.comb}{lm} (X, seed=NULL)
}
\arguments{
  \item{X}{ An \eqn{n \times p}{n by p} matrix or data frame.  For the lm
method, X is a linear regression model, and \code{random.lin.comb} is 
applied to \code{X$model[,-1]} }
  \item{seed}{ If not NULL, use this to initialize the random number generator}
}
\details{
  The matrix \eqn{X} is first scaled and centered.  \eqn{a} is computed to have random
  uniform components on \eqn{(-1,1)}.
}
\value{
A vector of length \eqn{n}.
}
\author{Sanford Weisberg, \email{sandy@stat.umn.edu} }
\seealso{ \code{\link{pure.error.anova}} }
\examples{
random.lin.comb(matrix(rnorm(9),ncol=3))
}
\keyword{ manip}% at least one, from doc/KEYWORDS
