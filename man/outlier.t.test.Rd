\name{outlier.t.test}
\alias{outlier.t.test}
\alias{outlier.t.test.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bonferroni test for outliers in linear models }
\description{
  Compute significance levels for the mean-shift outlier model using the
  Bonferroni inequality
}
\usage{
outlier.t.test(m, order=TRUE, bound=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{ A model of type lm. }
  \item{order}{If TRUE, order the cases according to the p-value.  If
  FALSE, don't order.}
  \item{bound}{Ignore cases with p-value bigger or equal to this value. }
}
\details{
  Returns length(res)*2*(1-pt(abs(res),df)), where res = rstandard(m) is
  the vector of Studentized residuals.  These are two-sided Bonferroni
  significance levels for testing a single outlier.
}
\value{
  A data frame with columns giving the value of the studentized residual 
  and corresponding Bonferroni p-value, and one row for each case for
  which the bound is satisfied.
}
\references{ Weisberg, S. (2005). \emph{Applied Linear Regression}, third
edition, Wiley.}
\author{Sanford Weisberg, sandy@stat.umn.edu}
\examples{
data(cloud)
m1 <- lm(logb(Rain,2)~S+logb(C,2)+logb(P,2)+E+A,data=cloud,subset=-2)
outlier.t.test(m1)
}
\keyword{ regression }% at least one, from doc/KEYWORDS
