\name{lathe1}
\alias{lathe1}
\docType{data}
\title{Lathe data}
\description{
These data are the results of an experiment to study the performance of
cutting-tool material in cutting steel on a lathe.  The two factors are
revolution speed and feed rate.  The response is tool life in minutes.
}
\format{This data frame contains the following columns:
 \describe{
    \item{Feed}{
Coded feed rate, coded as (actual feed rate -13)/6.  Feed is in thousandths
of an inch per revolution.
    }
    \item{Speed}{
Coded speed, coded as (actual speed -900)/300.  Speed is in feet per minute.
    }
    \item{Life}{
Life of tool until failure, minutes 
    }
 }
}

\source{
M. R. Delozier
}
\references{Weisberg, S. (2005). \emph{Applied Linear Regression}, 3rd edition. New York: Wiley, Problem 6.3.}
\examples{
head(lathe1)
}
\keyword{datasets}
