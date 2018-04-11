\name{randomLinComb}
\alias{randomLinComb}
\alias{randomLinComb.default}
\alias{randomLinComb.lm}
\title{Compute a random linear combination of the columns of a matrix or data frame }
\description{
Computes \eqn{Xa} where \eqn{X} is an \eqn{n \times p}{n by p} matrix, and \eqn{a} is
a random vector of length \eqn{p}.
}
\usage{
randomLinComb(X, seed = NULL)

\S3method{randomLinComb}{lm} (X, seed=NULL)
}
\arguments{
  \item{X}{ An \eqn{n \times p}{n by p} matrix or data frame.  For the lm
method, X is a linear regression model, and \code{randomLinComb} is 
applied to \code{model.matrix(X)} }
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
randomLinComb(matrix(rnorm(9),ncol=3))
}
\keyword{ manip}% at least one, from doc/KEYWORDS
