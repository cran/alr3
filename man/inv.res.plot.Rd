\name{inv.res.plot}
\alias{inv.res.plot}
\alias{inverse.response.plot}
%\alias{inverse.response.plot.formula}
\alias{inverse.response.plot.default}
\title{Inverse response plots to transform the response}
\description{
For a \code{lm} object, draws an inverse.response plot with the response \eqn{Y}{Y} on the
vertical axis and the fitted values \eqn{\hat{Y}}{Yhat}
on the horizontal axis.  Uses \code{nls} to
estimate  \eqn{\lambda}{lambda} in the function 
\eqn{\hat{Y}=b_0+b_1Y^{\lambda}}{Yhat = b0 + b1(Y)^(lambda)}.
Adds the fitted curve to the plot.
}
\usage{
inverse.response.plot(m, lambda=c(0,1),maxiter=100,xlab=NULL,...)

%inverse.response.plot.formula(m, data, subset, na.action,...)

inv.res.plot is an alias for inverse.response.plot.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{A lm regression object}
  \item{lambda}{A vector of values for lambda.  A plot will be produced 
  with curves corresponding to these lambdas and to the least squares 
  estimate of lambda}
  \item{xlab}{The horizontal axis label.  If NULL, it is constructed by the function.}
  \item{maxiter}{Passed to \code{nls}, the maximum number of iterations to
  be used to compute the optimal value of lambda}
%  \item{data}
%{See documentation for \code{\link{lm}}.}
%  \item{subset}
%{See documentation for \code{\link{lm}}.}
%  \item{na.action}
%{See documentation for \code{\link{lm}}.}
  \item{\dots}{Other arguments passed to \code{inv.tran.plot} and then to \code{plot}.}
}

\value{
As a side effect, a plot is produced with the response on the horizontal
axis and fitted values on the vertical axis.  
Several lines are added to be plot as
the ols estimates of the regression of \eqn{\hat{Y}}{Yhat} on
\eqn{Y^{\lambda}}{Y^(lambda)},
interpreting \eqn{\lambda}{lambda} = 0 to be natural logarithms.

Numeric output is a list with elements
  \item{lambda}{Estimate of transformation parameter for the response}
  \item{se }{Standard error of the estimate}
  \item{RSS}{The residual sum of squares at the minimum}
}
\seealso{\code{\link{inv.tran.plot}} and \code{\link{inv.tran.estimate}}, for which
this is just a convenient front-end, and \code{\link{nls}}.}

\references{S. Weisberg (2005), \emph{Applied
Linear Regression}, third edition, Wiley, Chapter 7}
\author{Sanford Weisberg, \code{sandy@stat.umn.edu}}

\examples{
data(highway)
highway$Sigs1 <- (round(highway$Sigs*highway$Len)+1)/highway$Len
attach(highway)
d <- data.frame(Rate=Rate,logLen=logb(Len,2),
                logADT=logb(ADT,2),logTrks=logb(Trks,2),
                Slim=Slim,Shld=Shld,logSigs1=logb(Sigs1,2))
attach(d)
m2 <- lm(Rate~logLen+logADT+logTrks+Slim+Shld+logSigs1,d)
inv.res.plot(m2,key=c(6,2))
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ regression}
