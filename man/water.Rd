\name{water}
\alias{water}
\non_function{}
\docType{data}
\title{California water}
\description{
Can Southern California's water supply in future years be
predicted from past data? One factor affecting water availability
is stream runoff.  If runoff could be predicted, engineers, planners
and policy makers could do their jobs more efficiently.
Multiple linear regression models have been used in this regard.
This dataset contains 43 years worth of precipitation measurements
taken at six sites in the Owens Valley ( labeled APMAM, APSAB,
APSLAKE, OPBPC, OPRC, and OPSLAKE), and stream runoff volume
at a site near Bishop, California.
}
\format{This data frame contains the following columns:
 \describe{
    \item{Year}{
 collection year
    }
    \item{APMAM}{
 Snowfall measurement site
    }
    \item{APSAB}{
 Snowfall measurement site
    }
    \item{APSLAKE}{
 Snowfall measurement site
    }
    \item{OPBPC}{
 Snowfall measurement site
    }
    \item{OPRC}{
 Snowfall measurement site
    }
    \item{OPSLAKE}{
 Snowfall measurement site
    }
    \item{BSAAM}{
 Stream runoff near Bishop, CA
    }
 }
}

\source{
Source:  http://www.stat.ucla.edu.}
\references{Weisberg, S. (2005). \emph{Applied Linear Regression}, 3rd edition. New York: Wiley, Problem 1.5.}
\examples{
data(dwp)
}
\keyword{datasets}
