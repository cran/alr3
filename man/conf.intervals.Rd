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
conf.intervals(object,level=.95,f=qnorm((1-level)/2))
}
\arguments{
  \item{object}{A regression object}
  \item{level}{one minus significance level, default = .95}
  \item{f}{multiplier for the intervals.  For \code{conf.intervals.default}, the
default choice of \code{f} is \code{qnorm((1-level)/2)}, using the normal distribution
to get intervals, and t-intervals for \code{conf.invtervals.lm}.
  }
}

\value{
  Returns a matrix, with one row for each estimated parameter, and three 
  columns, the first giving the coefficient estimate, then the
  lower limit of the confidence interval, 
  and then the upper limit.}
}
\author{Sanford Weisberg, sandy@stat.umn.edu}

\examples{
data(forbes)
m1 <- lm(Lpres ~ Temp, data = forbes)
conf.intervals(m1, level = .95) # .95 is the default.
}
\keyword{ regression }% at least one, from doc/KEYWORDS
