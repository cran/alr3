\name{sigma.hat}
\alias{sigma.hat}
\alias{sigma.hat.default}
\alias{sigma.hat.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Return the scale estimate for a regression model }
\description{
  This function provides a consistent method to return the estimated scale
  from a linear, generalized linear, nonlinear, or other model.
}
\usage{
sigma.hat(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A regression object of type lm, glm or nls}
}
\details{
For an lm or nls object, the returned quantity is the square root of the
estimate of \eqn{\sigma}{sigma}.  For a glm object, the returned quantity is
the square root of the estimated dispersion parameter.
}
\value{
  A nonnegative number
}
\author{Sanford Weisberg, \email{sandy@stat.umn.edu} }

\examples{
data(forbes)
m1 <- lm(Lpres ~ Temp, data=forbes)
sigma.hat(m1)
}
\keyword{ regression }% at least one, from doc/KEYWORDS
