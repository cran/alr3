\name{allshoots}
\alias{longshoots}
\alias{shortshoots}
\alias{allshoots}
\docType{data}
\title{Apple shoots data}
\description{
  Bland's Apple Shoot data.   allshoots includes all the data,
shortshoots just the short shoot data, and longshoots includes long shoots
only.
}
\format{This data frame contains the following columns:
 \describe{
    \item{Day}{
 days from dormancy
    }
    \item{n}{
 number of shoots sampled
    }
    \item{ybar}{
 average number of stem units
    }
    \item{SD}{
 within-day standard deviation
    }
    \item{Type}{
 1 if long shoots, 0 if shortshoots.
    }
 }
}

\source{Bland, J. (1978). A comparisonof certain aspects of ontogeny in the long
and short shoots of McIntosh apple during one annual growth cycle.  Unpublished 
Ph. D. dissertation, University of Minnesota, St. Paul, Minnesota.}
\references{Weisberg, S. (2005). \emph{Applied Linear Regression}, 3rd edition.
New York: Wiley, Section 5.3.}

\examples{
head(longshoots)
}
\keyword{datasets}
