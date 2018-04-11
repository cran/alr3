\name{alr3-deprecated}
\alias{alr3-deprecated}
%\alias{sigma.hat}
\alias{conf.intervals}
\alias{boot.case}
\alias{pure.error.anova}
\alias{delta.method}
\alias{powtran}
\alias{inv.tran.plot}
\alias{inv.tran.estimate}
\alias{inverse.response.plot}
\alias{inv.res.plot}
\alias{bctrans}
\alias{bctrans1}
\alias{lrt.bctrans}
%\alias{resid.curv.test}
%\alias{tukey.nonadd.test}
\alias{resplot}
\alias{residual.plots}
\alias{marginal.model.plot}
\alias{marginal.model.plots}
\alias{inf.index}
\alias{outlier.t.test}


\title{Deprecated Functions in the alr3 Package}

\description{
  These functions are provided for compatibility with older versions of alr3 only, and
may be removed eventually.
}

\usage{
%sigma.hat(...)

conf.intervals(...)

boot.case(...)

pure.error.anova(...)

delta.method(...)

powtran(...)

inv.tran.plot(...)

inv.tran.estimate(...)

inverse.response.plot(...)

inv.res.plot(...)

bctrans(...)

bctrans1(...)

lrt.bctrans(...)

%resid.curv.test(...)

%tukey.nonadd.test(...)

resplot(...)

residual.plots(...)

marginal.model.plot(...)

marginal.model.plots(...)

inf.index(...)

outlier.t.test(...)
}

\arguments{
  \item{\dots}{pass arguments down.}
}

\details{
\code{random.lin.comb} \code{pure.error.anova} are now synonyms for
\code{sigmaHat}, \code{randomLinComb} and \code{pureErrorAnova}, all of
which are still part of the alr3 package.

\code{conf.intervals} is a synonym for \code{confint} in the stats package.

All the remaining functions now shadow newer functions in
the car package.  In several cases, the arguments to the function
have changed as well, and so there is no guarantee that any
of these will work as expected.  Please substitute the corresponding
function from car.
}
