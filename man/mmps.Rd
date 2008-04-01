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
mmps(object,vars=~.,fitted=TRUE,layout=NULL,ask,...)

\method{mmp}{lm}(object, u , mean = TRUE, sd = FALSE, 
   xlab = deparse(substitute(u)), degree = 1, span = 2/3,key="topleft", 
   lineColors = c("blue", "red"), ...)
   
\method{mmp}{glm}(object, u, mean = TRUE, sd = FALSE, 
   xlab = deparse(substitute(u)), degree = 1, span = 2/3,key="topleft", 
   lineColors = c("blue", "red"), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{A regression object, of type either lm or glm,
   for which there is a \code{predict}
  method defined. }
\item{vars}{A one-sided formula.  A marginal model plot will be drawn for
  each variable on the right-side of this formula that is not a factor.  The
  default is \code{~.}, which specifies that all the terms in 
  \code{formula(object} will be used.  See examples below.}
\item{fitted}{If the default TRUE, then a marginal model plot in the direction
  of the fitted values or linear predictor of a generalized linear model will
  be drawn.}
\item{layout}{A reasonable layout for the plots in the window is
  determined by the program.  If you don't like the default you can set your
  own layout:  \code{c(2,3)} means two rows and three columns.}
\item{ask}{If TRUE, ask before clearing the graph window to draw more plots.}
\item{\dots}{Additional arguments passed from \code{mmps} to \code{mmp} and
  then to \code{plot}.
  Users should generally use \code{mmps}.}
\item{u}{ The quantity to be plotted on the horizontal axis.  The
  default is the predicted values \code{predict(object)}.  Can be any other
  vector of length equal to the number of observations in the object. }
\item{mean}{ If TRUE, compare mean smooths}
\item{sd}{ If TRUE, compare sd smooths }
\item{xlab}{ label for horizontal axis }
\item{degree}{Degree of the local polynomial, passed to \code{loess}.  The
  usual default for \code{loess} is 2, but the default here is 1.}
\item{span}{Span, the smoothing parameter for \code{loess}.  }
\item{key}{If not NULL, this gives the location of the legend, chosen from
  the default "topleft", or "bottomleft","bottomright" or "topright".}
\item{lineColors}{ colors for data and model smooth, respectively }
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
