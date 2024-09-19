\name{numDeltaMethod}
\Rdversion{1.1}
\alias{numDeltaMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate numerical delta method for non-linear predictions.
}
\description{
  Given a regression object and an independent prediction function (as a
  function of the coefficients), calculate the point estimate and
  standard errors 
}
\usage{
numDeltaMethod(object, fun, gd=NULL, conf.int=FALSE, level=0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A regression object with methods \code{coef} and \code{vcov}.
}
  \item{fun}{
    An independent prediction function with signature
    \code{function(coef, ...)}.
}
  \item{gd}{
    Specified gradients
}
  \item{conf.int}{
    Logical for whether to also calculate the confidence interval
}
  \item{level}{
    Numeric for the level of the confidence interval
}
  \item{\dots}{
    Other arguments passed to \code{fun}.
}
}
\details{
  A more user-friendly interface is provided by \code{predictnl}.
}
\value{
  \item{fit}{Point estimates}
  \item{se.fit}{Standard errors}
  \item{Estimate}{Point estimates}
  \item{SE}{Standard errors}
  \item{conf.low}{Lower confidence interval (if conf.int=TRUE)}
  \item{conf.high}{Upper confidence interval (if conf.int=TRUE)}
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  See Also \code{\link{predictnl}}.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, fun, ...) 
{
    coef <- coef(object)
    est <- fun(coef, ...)
    Sigma <- vcov(object)
    gd <- grad(fun, coef, ...)
    se.est <- as.vector(sqrt(colSums(gd * (Sigma \%*\% gd))))
    data.frame(Estimate = est, SE = se.est)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
