\name{fuel2001}
\alias{fuel2001}
\non_function{}
\docType{data}
\title{Fuel consumption}
\description{
Data on motor fuel consumption and related variables, for the year 2001.
The unit is a state in the
United States or the District of Columbia.  Data are for 2001, unless noted.
}
\format{This data frame contains the following columns:
 \describe{
    \item{Drivers}{
 Number of Licensed drivers in the state
    }
    \item{FuelC}{
 Gasoline sold for road use (1000s of gal.)
    }
    \item{Income}{
 Per capita personal income (year 2000)
    }
    \item{Miles}{
 Miles of Federal-aid highway miles in the state
    }
    \item{MPC}{
 Estimated miles driven per capita
    }
    \item{Pop}{
 Population age 16 and over
    }
    \item{Tax}{
 Gasoline state tax rate, cents per gallon
    }
    \item{State}{
 State name
    }
 }
}
\source{
http://www.fhwa.dot.gov/ohim/hs01/index.htm}

\references{
Weisberg, S. (2005). \emph{Applied Linear Regression}, third edition.  New
York:  Wiley.}

\examples{
data(fuel2001)
fuel2001$Dlic <- 1000*fuel2001$Drivers/fuel2001$Pop
fuel2001$Fuel <- 1000*fuel2001$FuelC/fuel2001$Pop
fuel2001$Income <- fuel2001$Income/1000
fuel2001$logMiles <- log(fuel2001$Miles,2)
pairs(Fuel~Tax+Dlic+Income+logMiles,data=fuel2001)
}
\keyword{datasets}
