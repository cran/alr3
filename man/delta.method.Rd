\name{delta.method}
\alias{delta.method}
\alias{delta.method.default}
\alias{delta.method.nls}
\alias{compute.delta.method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standard error of a nonlinear function of regression 
coefficients}
\description{
These functions use the delta method to get a first-order approximate 
standard error for a nonlinear function of a vector of random variables
with known or estimated covariance matrix.
}
\usage{
delta.method(object, g, parameter.prefix="b",print=TRUE)
compute.delta.method(Var, g, values, para.name, print=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A regression object of type lm, aov, glm or nls.  More
  generally, object can be any object for which \code{vcov(object)} returns
  a covariance matrix and \code{coef(object)} returns a vector of parameter
  estimates. }
  \item{g}{A quoted string that is the nonlinear function to be
  evaluated; see the details below. }
  \item{values}{The estimated values of the parameters.}
  \item{para.name}{Names of the parameters, as used in the expression
  \code{g}.}
  \item{parameter.prefix}{If \code{para.name} is not set, then parameter names will
  be generated to be the \code{parameter.prefix} followed by a number.  See details
  below.}
  
  \item{Var}{Estimated or population covariance matrix of the estimates}
  \item{print}{If TRUE, print estimates, else return a list}
}

\details{The delta method computes the estimate and large-sample standard
error of a nonlinear function of parameter estimates.  Suppose the names of
the parameters are \dQuote{b0}, \dQuote{b1},...,\dQuote{bp}, as would be the
case for an object of type \code{lm} or \code{glm} with parameter.prefix
left at its default value of \dQuote{b}.

If the value of \code{g} were
\code{"(b2+b6*350)/(4*b4)"}, then the function delta.method would extract
the coefficient estimates from \code{object}, and the covariance matrix from
\code{object}, and  pass them to \code{compute.delta.method}.

For compute.delta.method, the vector values of length \eqn{p}{p} has estimated
covariance matrix \eqn{V}{V}, a \eqn{p\times p}{p by} p matrix.  \code{g} 
is an expression of the form
\code{"(b2+b6*350)/(4*b4)"} that gives the function of values to be evaluated.
The routine evaluates \code{g} at the parameter values, and then computes \eqn{\dot{g}}{g'},
the derivative 
of \eqn{g}{g} with respect to the parameter names, also evaluated at the values,
and then computes the standard error, 
\deqn{[(\dot{g})'(\mbox{Var}) \dot{g}]^{1/2}}{(t(g') (Var) g')^(1/2)}

\code{delta.method} calls \code{compute.delta.method}, which can be used directly by
providing the nonlinear function, its values, and its estimated covariance matrix.
}

\value{
  \item{Estimate }{estimate of g}
  \item{se }{Standard error of g}
  \item{func}{g}
}
\references{S. Weisberg (2005), \emph{Applied
Linear Regression}, third edition, Wiley, Section 6.1.2}
\author{Sanford Weisberg, \email{sandy@stat.umn.edu}}

\examples{
# cakes is a data frame with response Y, predictors X1 X2 
data(cakes) 
m1 <- lm(Y~ X2 + I(X2^2), data = cakes) # quadratic polynomial
delta.method(m1, "-b1/(2*b2)") # X2 that maximizes the quadratic
# second order polynomial in two predictors:
m2 <- lm(Y ~ X1 + X2 + I(X1^2) + I(X2^2) + X1:X2, data=cakes)
# Find X1 to maximize Y when X2=350:
delta.method(m2,"(b1+b5*350)/(-2*b3)")
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ regression }
