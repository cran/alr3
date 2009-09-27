\name{alrweb}
\alias{alrweb}
\alias{alrprimer}
\alias{alrerrata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Access to the Applied Linear Regression website }
\description{
  This function will access the website for Applied Linear Regression} 

\usage{
alrweb(what = c("webpage", "errata", "primer"))

alrprimer()

alrerrata()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{what}{ A character string indicating what page to open. The default 
          "webpage" will open the main webpage, "errata" displays the
          Errata sheet for the book, and "primer" fetches and displays
          the primer for R.
 }
}

\value{
  Either a webpage or a pdf document is displayed.  alrprimer() and alrerrata()
call alrweb() with the appropriate arguments.
}

\author{ Sanford Weisberg, based on the function UsingR in the UsingR
package by John Verzani }


\examples{

\dontrun{alrweb()}
}
\keyword{ interface }% at least one, from doc/KEYWORDS
