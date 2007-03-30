\name{plot.bctrans}
\alias{plot.bctrans}
\alias{summary.bctrans}
\alias{coef.bctrans}
\alias{vcov.bctrans}
\title{ Plot and summary helpers for bctrans objects}
\description{
  Extract, view, and plot useful summaries of a \sQuote{bctrans} object for
  transformations to normality.
}
\usage{

plot.bctrans(x,y=x$Y, lambda=coef(x),family=x$family,plot=pairs,...)

summary.bctrans(object, digits=4, ...)

coef.bctrans(object,...)

vcov.bctrans(object,...)
}
\arguments{
\item{x}{An object created by \code{bctrans}.}
\item{y}{If not NULL, add this variable to the scatterplot matrix}
\item{object}{An object created by \code{bctrans}.}
\item{lambda}{A vector of transformation parameters.}
\item{family}{A family of transformations, probably \code{"box.cox"} or
\code{"yeo.johnson"}}
\item{plot}{Name of the plotting function to view the results, usually
\code{pairs} or similar.}
\item{digits}{Accuracy of summary output for \sQuote{bctrans} objects.}
\item{...}{Arguments passed to other methods.}
}

\value{
The \code{plot} method is used for its side effect of producing a scatterplot matrix.  
The helper \code{summary}
summarizes a bctrans object. This helper calls 
\code{\link{lrt.bctrans}}, but the latter can be
called directly.  The helpers \code{coef} and \code{vcov}
return, respectively, the estimated value of lambda, and an estimate of
its covariance matrix (based on a profile log-likelihood).
}

\references{
Weisberg, S. (2005) \emph{Applied Linear Regression}, third edition. Wiley.
}
\author{ Sanford Weisberg, \email{sandy@stat.umn.edu} }

\seealso{
\code{\link{bctrans}}, \code{\link{lrt.bctrans}},
\code{\link{pairs}},
%\code{\link{scatterplot.matrx}},
\code{\link{inv.res.plot}}
 }

\examples{
data(highway)
b <- highway[,c(8,1,2,10,5)] # select interesting columns
summary(ans <- bctrans(b,family="yeo.johnson")) # zeros ==> use yeo.johnson
b$Sigs <- (round(b$Sigs*b$Len)+1)/b$Len # redefine so no zeroes
summary(ans <- bctrans(b)) # fit with box.cox
b$Sigs <- (round(b$Sigs*b$Len)+1)/b$Len # redefine so no zeroes
summary(ans <- bctrans(b)) # fit with box.cox
lrt.bctrans(ans,lrt=list(c(0,-1,1,0)))
plot(ans,family="power") # plot, but use ordinary powers
b <- cbind(b,powtran(ans)) # add transformed variables to data frame
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
