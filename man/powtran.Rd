\name{powtran}
\alias{powtran}
\alias{powtran.default}
\alias{powtran.matrix}
\alias{powtran.data.frame}
\alias{powtran.bctrans}

\title{ Power transformations }
\description{
  \code{powtran} computes members of families of transformations indexed by one
parameter, the Box-Cox power family, or the Yeo and Johnson (2000) family, or the
basic power family, interpreting zero power as logarithmic. 
The family can be modified to have Jacobian one, or not, except for the basic
power family.
}
\usage{
powtran(U, lambda,family,modified)

powtran.default(U, lambda,family="box.cox",modified=TRUE)

powtran.bctrans(U, lambda=coef(U), family=U$family,modified=FALSE)

powtran.matrix(U, lambda, family="power",modified=FALSE)

powtran.data.frame(U, lambda, family="power",modified=FALSE)
}

\arguments{
  \item{U}{U can be a vector, a matrix, a data.frame, or an object of type
  bctrans.  Missing values are permitted.  U must 
  be strictly positive unless family=\dQuote{yeo.johnson}.}
  \item{lambda}{ The transformation parameter, typically a value between -2 and 2.
  For U of type matrix, data.frame, or bctrans, lambda must be a vector with the
  same number of elements as U has columns.  lambda has no default, except for
  objects for type bctrans where the default is the value computed by bctrans.}
  \item{family}{ The transformation family to be used. See the details below.}
  \item{modified}{TRUE to divide by the Jacobian, as need to transform the response, 
  FALSE for no division.  The default is TRUE except for objects of type bctrans,
  where the default is FALSE.}
}

\details{
  The Box-Cox family of \emph{scaled power transformations}, \code{family="box.cox"}, 
  equals \eqn{(U^{\lambda}-1)/\lambda}{(U^(lambda)-1)/lambda} 
  for \eqn{\lambda}{lambda} \eqn{\neq}{not equal to} zero, and 
  \eqn{\log(U)}{log(U)} if \eqn{\lambda =0}{lambda = 0}.
  
  If \code{family="yeo.johnson"} then the Yeo-Johnson transformations are used.  
  This is is Box-Cox transformation of \eqn{U+1} for nonnegative values, 
  and of \eqn{|U|+1} with parameter \eqn{2-\lambda}{2-lambda} for \eqn{U} negative.
  
  The basic power transformation is simply \eqn{U^{\lambda}}{U^lambda} 
  if lambda not zero, and \eqn{\log(U)}{log(U)} otherwise.
  
  If \code{modified=TRUE}, then the scaled transformations are divided by the 
  Jacobian, which is a function of the geometric mean.
  
}
\value{
  For U a vector, matrix or data.frame, returns an object of the same type as U of
  transformed values.  For class bctrans, returns an object of the same class and
  shape as U\$X.
}
\references{Weisberg, S. (2005).  \emph{Applied Linear Regression}, third edition.
New York:  Wiley, Chapter 7.

Yeo, In-Kwon and Johnson, Richard (2000).  A new family of power
transformations to improve normality or symmetry.  \emph{Biometrika}, 87,
954-959.
}
\author{ Sanford Weisberg, sandy@stat.umn.edu }

\seealso{\code{\link{inv.tran.plot}}, \code{\link{bctrans}},
\code{link{inv.res.plot}},\code{\link{tran.family}}. }
\examples{
data(ufcwc)
attach(ufcwc)
powtran(Height,lambda=0)
inv.tran.plot(Dbh,Height, lam = c(-1, 0, 1), family="box.cox")
summary(ans <- bctrans1(cbind(Dbh,Height))) # bctrans1 is used for a matrix input
plot(ans)
}
\keyword{ regression}% at least one, from doc/KEYWORDS
