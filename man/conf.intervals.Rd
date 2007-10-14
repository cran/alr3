\name{conf.intervals}
\alias{conf.intervals}
\alias{conf.intervals.default}
\alias{conf.intervals.lm}
\title{Compute marginal confidence intervals for regression estimates}
\description{
Thus function computes confidence intervals coefficient estimates
for a regression object.
}
\usage{
conf.intervals(...)
}
\arguments{
  \item{...}{Arguments are passed to \code{confint} in the stats package}
  }


\value{
This function calls \code{confint} in the stats package to compute
confidence intervals of coefficient estimates.}

\seealso{See also \code{\link{confint}}}
\author{Sanford Weisberg, sandy@stat.umn.edu}

\keyword{ regression }% at least one, from doc/KEYWORDS
