\name{residual.plots}
\alias{residual.plots}
\alias{residual.plots.lm}
\alias{resplot}
\alias{resid.curv.test}
\alias{tukey.nonadd.test}
\title{Residual plots and curvature tests for linear model fits}
\description{
Plots the residuals versus each term in a mean function and versus
fitted values.  Also computes a curvature test for each of the plots
by adding a quadratic term and testing the quadratic to be zero.  This
is Tukey's test for nonadditivity when plotting against fitted values.
}
\usage{
### This is a generic function with only one required argument:

residual.plots (m, ...)

### When the first argument is a linear model (of class lm), the form of the
### function is

residual.plots.lm(m,tukey=TRUE,exclude=NULL,plot=TRUE,
     layout=NULL,ask,...)

### The following are three related functions:    

resplot(m,varname="tukey",type="pearson",
                    plot=TRUE,add.quadratic=TRUE,...)
                    
resid.curv.test(m,varname)

tukey.nonadd.test(m)
}
\arguments{
  \item{m}{ \code{lm} regression object }
  \item{tukey}{If TRUE, draw plot of residuals versus fitted values and compute
  Tukey's test of non-additivity.}
  \item{exclude}{A list of predictor numbers for which plots will NOT be drawn.  If
  NULL, all plots are drawn.}
  \item{layout}{ If set to a value like \code{c(1,1)} or \code{c(4,3)}, the layout
  of the graph will have this many rows and columns.  If not set, the program will
  select an appropriate layout.  If the number of graphs exceed nine, you must 
  select the layout yourself, or you will get a maximum of nine per page.}
  \item{ask}{If TRUE, ask the user before drawing the next plot; FALSE if don't ask.}
  \item{\dots}{\code{residual.plots} passes these arguments to \code{resplot}.  
  \code{resplot}
  passes them to \code{plot}. }
  \item{varname}{Quoted variable name for the horizontal axis, 
  \code{"tukey"} by
  default for Tukey's test and the plot versus fitted values.}
  \item{type}{Type of residuals to be used.  Pearson residuals are 
  appropriate for \code{lm} objects since there are equivalent to ordinary residuals
  with ols and correctly weighted residuals with wls.}
  \item{add.quadratic}{if TRUE, fits the quadratic regression of the
  vertical axis on the  horizontal axis.}
  \item{plot}{If TRUE, draw the plot(s).}
  }
\details{
\code{resplot} draws ONE residual plot and returns the curvature test and its
significance level.  \code{residual.plots} draws all residuals plots, versus
each term in the model and versus fitted values,  and returns all the
curvature tests.  \code{residual.plots} should not be used with factors.  
\code{resid.curv.test} computes the curvature test only.
}
\value{
Returns a data.frame with one row for each plot drawn, one column for 
the curvature test statistic, and a second column for the corresponding
p-value.  This function is used primarily for its side effect of drawing
residual plots.
}
\references{S. Weisberg (2005), \emph{Applied
Linear Regression}, third edition, Wiley, Chapter 8}
\author{Sanford Weisberg, \email{sandy@stat.umn.edu}}

\seealso{See Also \code{\link{lm}}}


\examples{
data(highway)
highway$Sigs <- (round(highway$Sigs*highway$Len)+1)/highway$Len
attach(highway)
d <- data.frame(Rate=Rate,logLen=logb(Len,2),
                logADT=logb(ADT,2),logTrks=logb(Trks,2),
                Slim=Slim,Shld=Shld,logSigs1=logb(Sigs,2))
attach(d)
m2 <- lm(Rate~logLen+logADT+logTrks+Slim+Shld+logSigs1,d)
residual.plots(m2)
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
