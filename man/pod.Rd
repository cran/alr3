\name{pod}
\alias{pod}
\alias{pod.lm}
\alias{pod.formula}
\alias{print.pod}
\alias{summary.pod}
\alias{coef.pod}
\alias{anova.pod} 
%\alias{podanova}
\alias{deviance.pod} 
\alias{resid.pod}
\alias{formula.pod}
\alias{fitted.pod}
\alias{podresponse}
\alias{df.residual.pod}
\alias{predict.pod}
\alias{plot.pod}
\alias{plot.pod.lm}
\alias{podnls.fit}
\alias{residuals.pod}
\alias{vcov.pod}

\title{ Fit partial one-dimensional, or POD models, based on a linear model }
\description{
A partial one-dimensional model, or a POD model, provides a concise description of
a regression model with many predictors and one grouping variable.  It requires a 
nonlinear regression fit.
}
\usage{
## This is a generic function with different arguments depending on the class of the
## first argument.  The generic form is

pod (x, ...)

## If the first argument to pod is a formula suitable for defining a linear model,

\method{pod}{formula}(formula, data = sys.parent(), group, subset, weights, 
    na.action, mean.function = c("pod", "common", "parallel", 
        "general"), singular.ok = FALSE, contrasts = NULL, offset, 
    control = nls.control(), ...) 


## If the first argument to pod is the result of a linear model fit, the following
## function call is used.  All the arguments for pod.formula can also be passed to
## pod.lm.

\method{pod}{lm} (x, group, mean.function, control, ...)   
    
## The following related function require a pod fit as the first argument:
    
\method{anova}{pod}(object,scale=0,test="F",...)

\method{plot}{pod}(x, colors=1:nlevels(x$group),
  pch=1:nlevels(x$group),key="topleft",identify=FALSE,
  xlab="Linear Predictor", ylab=as.character(c(formula(x)[[2]])),
  ...)
 
\method{plot}{pod.lm}(x, colors=1:nlevels(x$group),
      pch=1:nlevels(x$group), key="topleft", identify=FALSE,
      xlab="Linear Predictor", ylab=as.character(c(formula(x)[[2]])),
      ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A linear regression formula, as used in lm, or the results of a
call to pod or a lm fit.  See details below.}
  \item{x}{The result of a lm fit from which the formula and the
arguments data, subset, weights, na.action, singular.ok, contrasts and offset will
be taken.}
  \item{data}{ An optional data frame for the data to be used}
  \item{group}{ The name of a grouping variable (not quoted) that defines the
groups; see details below.}
  \item{mean.function}{Which mean function should be fit?  
The default is \dQuote{pod}, that
fits the partial one-dimensional mean function.  The other options are 
\dQuote{common}, which
fits a linear model with no grouping effects; \dQuote{parallel} fits a parallel
within-group regression, and \dQuote{general}, available in pod but not plot.pod,
fits a separate coefficeint for each
term in the model for each level of the grouping variable.}
  \item{subset}{ As in lm, the subset of cases to be used in the analysis }
  \item{weights}{Weights will be used in fitting non-pod models.  Since this
  argument is not supported for nls models, weights are ignored for fitting
  pod models.  If nls is ever updated, then pod models will correctly use
  the weights.}
  \item{na.action}{ At present, only na.omit is likely to work. }
  \item{singular.ok}{ The default equal to \code{FALSE} is recommended.}
  \item{contrasts}{ Same as in lm }
  \item{offset}{ Same as in lm }
  \item{control}{A pod model is fit using the nonlinear least squares routine nls.
This routine is very sensitive to starting values and other parameters set in the
algorithm.  This routine selects starting values for you that will often work very
well, but in some problems the user may want to change the defaults to nls program 
using the \code{\link{nls.control}} function.}
  \item{\dots}{In pod, other arguments passed to nls, such as control 
  parameters. In pod.anova, there are two additional arguments }
  \item{object}{The result of a call to pod.}
  \item{scale}{Used for test in anova}
  \item{test}{Default is to compute F tests.}
  \item{colors}{Colors for groups in the pod plot}
  \item{pch}{Ploting symbol for the groups in the pod plot}
  \item{identify}{If TRUE, clicking the mouse on a graph will print the case name
  of the nearest point.  This continues until turned off (by pushing the escape
  key, among other ways of doing this).}
  \item{key}{The default is \code{"topleft"}, in which case a legend is added to
the top left corner of the plot; other choices include \code{"bottomright"}.
If \code{key} is a vector of two coordinates, the legend is drawn at the coordinates
specified. If \code{key} is FALSE, no key is drawn; if TRUE, you can place the
key interactively by clicking on the plot.}
  \item{xlab}{Horizontal label, optional}
  \item{ylab}{Vertical label, optional}
}
\details{
Suppose we start with a linear mean function specified by \code{y ~ x1 + x2 + x3}, 
where the right-side variables can be any valid R variables, such as 
transformations, but NOT
factors or interactions (if you want to include these, you need to create the
dummy variables yourself).  The right-hand side variables must also be linearly
independent.  We also specify a grouping variable z with, say, g levels.  Let 
\eqn{G_i}{Gi} be a dummy variable with values equal to 1 when \eqn{z=i}, and 
zero otherwise.  The pod mean function is then has the nonlinear mean function
\deqn{E(y|x,z) = \beta_0 + \beta'x + \sum_{j=2}^g G_j(\theta_{0j} + \theta_{ij}\beta'x)}{
b0+ b'x + sum(j=2,g) Gj(th0j + th1j(b'x))}
This is a nonlinear mean function that specifies that the response depends on the predictors
only through one linear combination, that the dependence is linear, but the slope
and intercept may be different for each group.

The pod mean function is easily fit using \code{nls}.  For example, if z has two levels,
a nonlinear mean function that would work is
\code{y ~ b0 + b1*x1 + b2*x2 + b3*x3 + G2*(th02 + th12*(b1*x1 + b2*x2 + b3*x3))}.
Starting values can be determined by fitting the parallel mean function, for this is a
linear mean function and easily fit.

The function \code{pod} automates this process; you need only input the linear
part of the mean function, and \code{pod} does the rest.

The function \code{anova.pod} is a little different from the generic anova 
function.  It  
fits and compares all four mean functions that can be fit with the
pod function:  (1) no group effect; (2) parallel
within group method; (3) pod mean function; (4) mean function with all main effects and 
two-factor
interactions with the grouping variable.

The function \code{plot} will plot the response versus the fitted linear 
combination for either the common, parallel or pod models.  There is no 2D plot for
the general model.

Objects of class pod have methods defined for the generic functions
\code{print},
\code{summary},
\code{coef},
\code{anova},
\code{deviance},
\code{vcov},
\code{resid},
\code{formula},
\code{fitted},
\code{df.residual} and
\code{predict}.  \code{podresponse} returns the values of the response.
}
\value{
If \code{mean.function="pod"}, a list of class pod, with the following components
  \item{nls.fit}{The usual output from a \code{nls} fit.}
  \item{linear.part }{The estimated linear predictor \eqn{\beta'x}{bhat'x}.}
  \item{call}{The original call to pod}
  \item{group}{The grouping variable}
  
If \code{mean.function} has any other value, an object of class c("pod.lm", "lm")
is returned.  Since the only method for pod.lm objects is a plotting method, 
except for plotting these objects behave like lm objects.  
}
\references{Cook, R. D. and Weisberg, S. (2004).  Partial-one dimensinal models.
\emph{American Statistician}, 58, 110-116.

Weisberg, S. (2005) \emph{Applied Linear Regression}, third edition.  New York: 
Wiley.}
\author{Sanford Weisberg, sandy@stat.umn.edu }

\seealso{See Also \code{\link{nls}}, \code{\link{lm}}, \code{\link{nls.control}} }
\examples{
head(ais)
m1 <- pod(LBM ~ Ht + Wt + RCC, data= ais, group= Sex)
anova(m1) # compare four models
plot(m1) # draw the plot
m2 <- update(m1, mean.function="parallel")
plot(m2)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ regression}
\keyword{hplot}
