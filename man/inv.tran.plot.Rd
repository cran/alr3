\name{inv.tran.plot}
\alias{inv.tran.plot}
\alias{inv.tran.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Choose a predictor transformation visually or numerically }
\description{
\code{inv.tran.plot}
draws a two-dimensional scatterplot of \eqn{Y}{Y} versus 
\eqn{X}{X}, along with the OLS
fit from the regression of \eqn{Y}{Y} on 
\eqn{(X^{\lambda}-1)/\lambda}{(Y^(lam)-1)/lam}.  \code{inv.tran.estimate}
find the nonlinear least squares estimate of \eqn{\lambda}{lambda} and its
standard error.
}
\usage{
inv.tran.plot(x,y,lambda=c(-1,0,1),lty=1:(1+length(lambda)),
        col=rainbow(length(lambda)+1),xlab=deparse(substitute(x)),
        ylab=deparse(substitute(y)),family="box.cox",optimal=TRUE,
        key="topleft",...)

inv.tran.estimate(x,y,family="box.cox",...)
}

\arguments{
  \item{x}{The predictor variable }
  \item{y}{The response variable }
  \item{lambda}{The powers used in the plot.  The optimal power than minimizes
  the residual sum of squares is always added unless optimal is FALSE. }
  \item{family}{The transformation family to use, \code{"box.cox"},
  \code{"yeo.johnson"}, or a user-defined family.}
  \item{optimal}{Include the optimal value of lambda?}
  \item{lty}{line types corresponding to the powers}
  \item{col}{color corresponding to the powers}
  \item{key}{The default is \code{"topleft"}, in which case a legend is added to
the top left corner of the plot; other choices include \code{"bottomright"}.
If \code{key} is a vector of two coordinates, the legend is drawn at the coordinates
specified. If \code{key} is FALSE, no key is drawn; if TRUE, you can place the
key interactively by clicking on the plot.}
  \item{xlab}{Label for the horizontal axis.}

  \item{ylab}{Label for the vertical axis.}

  \item{\dots}{ additional arguments passed to other methods.}
}
\value{
\code{inv.tran.plot}
  returns a graph and a data.frame with \eqn{\lambda}{lam} in the 
  first column, and the residual sum of squares from the regression
  for that \eqn{\lambda}{lam} in the second column.

  \code{inv.tran.estimate} returns a list with elements \code{lambda} for the
  estimate, \code{se} for its standard error, and \code{RSS}, the minimum
  value of the residual sum of squares.  
}
\seealso{\code{\link{powtran}}, \code{\link{inverse.response.plot}},
\code{\link{optimize}}, \code{\link{nls}}}

\references{
Weisberg, S. (2005).  \emph{Applied Linear Regression}, third edition.
New York:  Wiley. }
\author{Sanford Weisberg, \email{sandy@stat.umn.edu} }

\examples{
data(baeskel)
attach(baeskel)
inv.tran.plot(Sulfur,Tension,key=c(.6,450))
ans <-inv.tran.estimate(Sulfur,Tension)
# redraw the plot, including the nls estimate
inv.tran.plot(Sulfur,Tension,lambda=c(ans$lambda,-1,0,1),key=c(.6,450))
}    
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{regression}
