\name{salary}
\alias{salary}
\docType{data}
\title{Salary data}
\description{
Salary of faculty in a small Midwestern college in the early 1980s.
}
\format{This data frame contains the following columns:
 \describe{
    \item{Degree}{
 Degree, 1 if PhD, 0 if Masters
    }
    \item{Rank}{
 Rank (1 = Asst Prof, 2 Assoc Prof, 3 Prof)
    }
    \item{Sex}{
 1 if female, 0 if male
    }
    \item{Year}{
 Years in current rank
    }
    \item{YSdeg}{
 Years since highest degree earned
    }
    \item{Salary}{
 dollars per year
    }
 }
}

\source{Sanford Weisberg}
\references{Weisberg, S. (2005). \emph{Applied Linear Regression}, 3rd edition. New York: Wiley, Problem 6.10.}
\examples{
head(salary)
}
\keyword{datasets}
