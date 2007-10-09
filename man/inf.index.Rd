\name{inf.index}
\alias{inf.index}
\alias{inf.index.lm}
\title{Influence index plot}
\description{
Provides index plots of Cook's distances,
leverages, Studentized residuals and outlier significance levels for a 
regression object.
}
\usage{
inf.index(m,cooks.distance,rstudent,outlier.t.test,leverages,...)

\method{inf.index}{lm}(m, cooks.distance = TRUE, rstudent = TRUE, 
    outlier.t.test = TRUE, leverages = TRUE, ...) 
}
\arguments{
  \item{m}{A regression object}
  \item{cooks.distance}{If TRUE, plot Cook's distances}
  \item{rstudent}{If TRUE, plot the Studentized residuals}
  \item{outlier.t.test}{If TRUE, index plot n*2*qt(-abs(rstandard),m\$df-1), the
two-tailed significance levels for the outlier test based on the Bonferroni
inequality}
  \item{leverages}{If TRUE, plot the leverages}
  \item{\dots}{Arguments passed to \code{plot}}
}

\value{
Used for its side effect of producing a graph.
}
\references{S. Weisberg (2005), \emph{Applied
Linear Regression}, third edition, Wiley, Chapter 8}
\author{Sanford Weisberg, \email{sandy@stat.umn.edu}}

\seealso{ \code{\link{cooks.distance}} }
\examples{
data(rat)
names(rat)
m1 <- lm(y~BodyWt+LiverWt+Dose,data=rat)
inf.index(m1)
}

\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
