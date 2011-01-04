\name{alrWeb}
\alias{alrWeb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Access to the Applied Linear Regression website }
\description{
  This function will access the website for Applied Linear Regression} 

\usage{
alrWeb(page = c("webpage", "errata", "primer"), script)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{page}{ A character string indicating what page to open. The default 
          "webpage" will open the main webpage, "errata" displays the
          Errata sheet for the book, and "primer" fetches and displays
          the primer for R.}
  \item{script}{If set, this will open an R script for one of the chapters
          in the book in a browser window, or a script from one of the
          chapters in the R Primer.  Possible values are "chapter1",
          "chapter2", ... "chapter12", "appendix" for chapters and
          appendix from the book, or
          "primer0", "primer1", ... "primer12", "primerappendix"
          for material corresponding to the R primer.} 
}

\value{
  Either a webpage or a pdf document is displayed.  This function gives quick
access to the website for the book and in particular to the R primer.
}

\author{ Sanford Weisberg, based on the function UsingR in the UsingR
package by John Verzani }


\examples{

\dontrun{alrWeb()}
}
\keyword{ interface }% at least one, from doc/KEYWORDS
