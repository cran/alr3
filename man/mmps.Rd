\name{mmps}
\alias{mmps}
\alias{mmp}
\alias{mmp.lm}
\alias{mmp.glm}
\alias{marginal.model.plot}
\title{Marginal model plot }
\description{
For a regression object, plots the response on the vertical axis versus
a linear combination \eqn{u} of terms in the mean function on the horizontal
axis.  Added to the plot are a \code{loess} smooth for the graph, along with
a loess smooth from the plot of the fitted values on \eqn{u}.
}
\usage{
mmps(object,exclude=NULL,layout=NULL,ask,...)

mmp.lm(object, u = predict(object), mean = TRUE, sd = FALSE, 
   label = deparse(substitute(u)), degree = 1, span = 2/3, 
   colors = c("blue", "red"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A regression object, of type either lm or glm,
   for which there is a \code{predict}
  method defined. }
  \item{exclude}{A vector of indices of variables for which a plot is not
  required.  For example, \code{c(2,4)} would exclude the second and fourth
  terms in the model, excluding the intercept.}
  \item{layout}{A reasonable layout for the plots in the window is
  determined by the program.  If you don't like the default you can set your
  own layout:  \code{c(2,3)} means two rows and three columns.}
  \item{ask}{If TRUE, ask before clearing the graph window to draw more plots.}
  \item{u}{ The quantity to be plotted on the horizontal axis.  The
  default is the predicted values \code{predict(object)}.  Can be any other
  vector of length equal to the number of observations in the object. }
  \item{mean}{ If TRUE, compare mean smooths}
  \item{sd}{ If TRUE, compare sd smooths }
  \item{label}{ label for horizontal axis }
  \item{degree}{Degree of the local polynomial, passed to \code{loess}.  The
  usual default for \code{loess} is 2, but the default here is one.}
  \item{span}{Span, the smoothing parameter for \code{loess}.  }
  \item{colors}{ colors for data and model smooth, respectively }
  \item{\dots}{ additional arguments passed to \code{loess} and \code{plot} }
}   
\details{
\code{mmp} draws one marginal model plot.  \code{mmps} draws all marginal model plots,
versus each of the terms in the model and versus fitted values.  \code{mmps}
should not be used if some of the terms in the model are factors or
interactions.
}
\value{
Used for its side effect of producing plots. 
}

\seealso{\code{\link{loess}}, \code{\link{plot}} }

\references{ S. Weisberg (2005), \emph{Applied
Linear Regression}, third edition, Wiley, Chapter 8}
\author{Sanford Weisberg, \email{sandy@stat.umn.edu}}

\examples{
data(ufcwc)
c1 <- lm(Height~Dbh,ufcwc)
mmp(c1,ufcwc$Dbh,label="Diameter, Dbh",color=c("black","black"))
mmps(c1)

} 
\keyword{hplot }% at least one, from doc/KEYWORDS
\keyword{regression}
