\name{transact}
\alias{transact}
\docType{data}
\title{Transaction data}
\description{
Data on transaction times in branch officies of a large Australian bank.
}
\format{This data frame contains the following columns:
 \describe{
    \item{T1}{
 number of type 1 transactions
    }
    \item{T2}{
 number of type 2 transactions
    }
    \item{Time}{
 total transaction time, minutes
    }
 }
}

\source{
Cunningham, R. and Heathcote, C.
(1989), Estimating a non-Gaussian regression model with multicollinearity.
Australian Journal of Statistics, 31,12-17.}
\references{Weisberg, S. (2005). \emph{Applied Linear Regression}, 3rd edition. New York: Wiley, Section 4.6.1.}
\examples{
head(transact)
}
\keyword{datasets}
