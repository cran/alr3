\name{tran.family.box.cox}
\alias{tran.family.box.cox}
\alias{tran.family}
\alias{tran.family.default}
\alias{tran.family.boxcox}
\alias{tran.family.yeo.johnson}
\alias{tran.family.YJ}
\alias{tran.family.power}
\alias{tran.family.basic}

\title{ Compute transformation families  }
\description{
  These are \emph{internal} functions that compute families of transformations.  The user
  accesses these functions using \code{\link{powtran}}.
}
\usage{
tran.family.box.cox(U, lambda, modified = TRUE)
tran.family.yeo.johnson(U,lambda,modified=TRUE) 
tran.family.power(U,lambda,modified=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{U}{ The variable to be transformed. Missing values are permitted.  U must 
  be strictly positive unless \code{family="yeo.johnson"}. }
  \item{lambda}{ The transformation parameter, typically a value between \eqn{-2} and 2. }
  \item{modified}{ TRUE to divide by the Jacobian, as need to transform the response, 
  FALSE for no division.}
}
\details{
  Uses either the Box-Cox power family, the Yeo-Johnson generalization of this family
  for non-strictly-positive variables, or basic power transformations. 
  
  You can add additional families by writing a function \code{tran.family.myfamily}
  with the same arguments as \code{tran.family.box.cox}, and then set
  \code{family="myfamily"} in the call to \code{powtran}.
  
    \code{tran.family.YJ} is an alias for \code{tran.family.yeo.johnson}.  
  \code{tran.family.boxcox} is an alias for \code{tran.family.box.cox}.
  \code{tran.family.basic} is an alias for \code{tran.family.power}.
}
\value{
Returns a vector of the same length as U with the transformed values.
}
\references{Weisberg, S. (2005).  \emph{Applied Linear Regression}, third edition.
New York:  Wiley, Chapter 7.

Yeo, In-Kwon and Johnson, Richard (2000).  A new family of power
transformations to improve normality or symmetry.  \emph{Biometrika}, 87,
954-959.}
\author{ Sanford Weisberg, sandy@stat.umn.edu }

\seealso{\code{\link{powtran}}. }
\examples{
data(ufcwc)
attach(ufcwc)
powtran(Height,lambda=0)
}
\keyword{ internal }% at least one, from doc/KEYWORDS
