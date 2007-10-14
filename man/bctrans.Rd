\name{bctrans}
\alias{bctrans}
\alias{bctrans1}
\alias{lrt.bctrans}
\alias{print.bctrans}
\title{ Univariate and multivariate transformations to normality }
\description{
Estimates multivariate power transformations to multinormality 
by a maximum likelihood-like method. The univariate case is obtained when only 
one variable is specified. 
}
\usage{
bctrans(formula, data = NULL, subset, na.action = na.omit, ...)

## Use when you have a matrix or data.frame:
bctrans1(X, Y = NULL, start = NULL, family = "box.cox", call=NULL,...) 

## S3 methods for class 'bctrans'
lrt.bctrans(object, lrt=NULL, ones=TRUE, zeroes=TRUE)
}
\arguments{
\item{formula}{A formula, giving the variables to be transformed. 
The formula can be \emph{one-sided},
of the form \code{~ X1+X2+X3}, or \emph{two-sided}, of
the form \code{Y~X1+X2+X3}.  In the latter case, the response is not 
used in the
transformation, but it will be used in the \code{plot} method.  If you
have previously computed a linear model fit, say \code{m1}, then you
can use \code{formula(m1)} to provide this argument.}
\item{data}{a data.frame (or list) from which the variables in the formula
          should be taken.}
\item{subset}{an optional vector specifying a subset of observations to be used.}
\item{na.action}{If set to na.omit, the default, missing values are permitted.  If
set to na.fail, missing values are  not permitted.}
\item{...}{In bctrans, these are additional arguments passed to \code{bctrans1},
and described below. In bctrans1, these are additional arguments passed
to the function maximizer \code{optim}.} 
\item{X}{A vector, matrix, or data frame whose
columns are to be transformed.}
\item{Y}{If present this vector will be part of the object created, and will
be used in drawing plots.  It is not used for finding transformations.}
\item{start}{Starting values for the power transformation parameters; 
if NULL (the default), univariate transformations will be computed and 
used as the start values. }
\item{family}{The family of transformations.  The most common is \code{"box.cox"}
  for the Box-Cox transformation.  The \code{"yeo.johnson"} transformations are used
  if some elements of X are negative or zero.  The family \code{"power"} is used only
  in the \code{plot.bctrans} and \code{add.trans} functions to give basic power
  transformations, which can't be normalized to have Jacobian one.  The argument
  lambda indexes the family of transformations \eqn{\psi(X,\lambda)}{psi(X,lambda)}}
\item{object}{In \code{lrt.bctrans}, the name of a \sQuote{bctrans} object.}
\item{lrt}{In the \code{lrt.bctrans} command, a list of vectors each of length equal
  to the number of columns in X.  A Likelihood ratio test that the transformation
  parameters equal each of these vectors will be performed.}
\item{ones}{In \code{lrt.bctrans}, if TRUE test all the transformation parameters equal to
  one against a general alternative.}
\item{zeroes}{In \code{lrt.bctrans}, if TRUE test all the transformation parameters equal to
  zero against a general alternative.}
\item{call}{Not to be set by the user.}
}
\details{
Given a matrix \eqn{X}{X} with columns 
\eqn{X_1\ldots,X_p}{X1, ..., Xp}, this routine selects transformation 
parameters \eqn{\lambda_1,\ldots,\lambda_p}{lambda1,...,lambdap} 
from a one-parameter family of transformations such 
that the transformed variables 
\eqn{\psi(X_1,\lambda_1),\ldots,\psi(X_p,\lambda_p)}{psi(X1,lambda1),\dots,psi(Xp,lambdap)}
are as close to multivariate normal as possible.

The function uses the family of transformation you specify. 
If you use the family \code{box.cox} to select a transformation,
it is usual to use standard power transformations in further calculations.
}
\value{
  \code{bctrans} returns an object of class \code{bctrans},
  which may be printed or summarized.  It is a list with components \code{X} containing
  the input data, \code{family} the family used, \code{start} the starting values, and, 
  \code{optim}, the
  results from a call to \code{\link{optim}}, the function optimizer used in the 
  routine.
}
\references{
Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations. \emph{Journal
of the Royal Statisistical Society, Series B}. 26 211-46.

Cook, R. D. and Weisberg, S. (1999).  \emph{Applied Regression Including Computing
and Graphics}.  Wiley.

Fox, J. (2002).  \emph{R and S-Plus Companion to Applied Regression}.  Sage. 

Velilla, S. (1993).  A note on the multivariate Box-Cox transformation to 
normality.  \emph{Statistics and Probability Letters}, 17, 259-263.

Weisberg, S. (2005) \emph{Applied Linear Regression}, third edition. Wiley.

Yeo, I. and Johnson, R. (2000).  A new family of
power transformations to improve normality or symmetry.  
\emph{Biometrika}, 87, 954-959.

}
\author{ A substantial part of this code is borrowed from the function
\code{box.cox.powers} in the \code{car} package, written by John Fox, and documented in
Fox (2002).  It is based on a similar function in Arc; see Cook and Weisberg (1999).
It was modified by Sanford Weisberg, \email{sandy@stat.umn.edu} and renamed 
\code{bctrans}. }

\seealso{
\code{\link{powtran}},
\code{\link{optim}},
\code{\link{pairs}},
%\code{\link{scatterplot.matrx}},
\code{\link{inv.res.plot}},
\code{\link{plot.bctrans}}
 }
\examples{
data(highway)
b <- highway[,c(8,1,2,10,5)] # select interesting columns
summary(ans <- bctrans1(b,family="yeo.johnson")) # zeros ==> use yeo.johnson
# or, compute using a formula and get the same answer.
summary(ans2 <- bctrans(~Len+ADT+Trks+Shld+Sigs,data=highway,family="yeo.johnson"))
# or, first fit an lm, and extract the formula
m1 <- lm(Rate~Len+ADT+Trks+Shld+Sigs,data=highway)
summary(ans3 <- bctrans(formula(m1),data=highway,family="yeo.johnson"))
# work with the response
b$Sigs <- (round(b$Sigs*b$Len)+1)/b$Len # redefine so no zeroes
summary(ans <- bctrans1(b)) # fit with box.cox
lrt.bctrans(ans,lrt=list(c(0,0,-1,1,0)))
plot(ans,family="power") # plot, but use ordinary powers
b <- cbind(b,powtran(ans)) # add transformed variables to data frame
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
