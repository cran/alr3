\name{boot.case}
\alias{boot.case}
\alias{boot.case.default}
\alias{boot.case.nls}
\alias{next.boot}
\alias{next.boot.default}
\alias{next.boot.nls}
\title{ Case bootstrap for regression models }
\description{
This routine does a case bootstrap resampling for regression models.  It
returns a matrix of the estimated coefficients from each of the bootstrap samples.
}
\usage{
boot.case(object, f=coef, B=999)
next.boot(object,sample)
}
\arguments{
  \item{object}{ Any regression object that works with \code{update} and with
the \code{subset} argument to \code{update}. }
  \item{f}{A function that will be applied to the updated regression object to compute
  the statistics of interest.  The default is \code{coef}, to return to regression
  coefficient estimates.}
  \item{B}{Number of bootstrap samples.}
  \item{sample}{A permutation of the integers from 1 to n that defines a
  bootstrap sample.}
}
\details{
This routine does the case-bootstrap described in the reference below.  Begin with a 
regression object.  For each of B bootstrap samples, sample the rows of the data matrix
with replacement, and recompute and save estimates.  For nls objects there may be
convergence problems in the bootstrap.  The routine will continue until convergence
is attained B times, or until there are 25 consecutive failures to converge. 
next.boot is an internal function that will update a model correctly, depending
on the class of the model object.
}
\value{
A matrix with B rows and rank(object) columns giving the bootstrap estimates.  These
can be summarized as needed using standard R/S-plus tools.
}
\references{
S. Weisberg (2005). \emph{Applied Linear Regression}, third edition.  New York:  Wiley,
Chapters 4 and 11.}
\author{ Sanford Weisberg, \email{sandy@stat.umn.edu}.  The error checking was written by
Lexin Li.}

\seealso{See Also  \code{\link{update}}}


\examples{
data(transact)
m1 <- lm(Time~ T1 + T2, data=transact)
betahat <- coef(m1)
betahat.boot <- boot.case(m1,B=99) # 99 bootstrap samples--too small to be useful
summary(betahat.boot)  # default summary
# bootstrap standard errors 
apply(betahat.boot,2,sd)
# bootstrap 95% confidence intervals
cl <- function(x) quantile(x,c(.025,.975))
apply(betahat.boot,2,cl)

}
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{internal}% __ONLY ONE__ keyword per line
