\name{pureErrorAnova}
\alias{pureErrorAnova}
\alias{pureErrorAnova.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Pure Error analysis of variance }
\description{
For a linear model object, finds the sum of squares for lack of fit 
and the sum of squares for pure error.  These are added to the standard anova table
to give a test for lack of fit.  If there is no pure error, then the regular anova
table is returned.
}
\usage{
### This is a generic function.  

pureErrorAnova(mod)

\method{pureErrorAnova}{lm}(mod)

### Methods for other than models for type lm have not been defined.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{an object of type \code{lm}}
}
\details{
  For regression models with one predictor, say \code{y ~ x}, this method fits
  \code{y ~ x + factor(x)} and prints the anova table.  With more than one predictor, it
  computes a random linear combination \eqn{L} of the terms in the mean function
  and then gives the anova table for \code{update(mod, ~.+factor(L))}.
}

\value{
  Returns an analsis of variance table.
}
\references{ Weisberg, S. (2005).  \emph{Applied Linear Regression}, third edition,
New York:  Wiley, Chapter 5.}
\author{ Sanford Weisberg, sandy@stat.umn.edu }

\seealso{  \code{\link{lm}} }
\examples{
x <- c(1,1,1,2,3,3,4,4,4,4)
y <- c(2.55,2.75,2.57,2.40,4.19,4.70,3.81,4.87,2.93,4.52)
m1 <- lm(y~x)
anova(m1)  # ignore pure error
pureErrorAnova(m1)  # include pure error

head(forbes)
m2 <- lm(Lpres~Temp, data=forbes)
pureErrorAnova(m2)  # function does nothing because there is no pure error

}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
