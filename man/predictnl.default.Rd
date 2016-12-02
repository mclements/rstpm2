\name{predictnl.default}
\Rdversion{1.1}
\alias{predictnl.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Default implementation of the \code{predictnl} generic function. 
}
\description{
  Given a regression \code{object} and a function \code{fun} that
  predicts values from the object, return the predicted values and the
  standard errors using the numeric delta method.
}
\usage{
predictnl.default(object, fun, newdata = NULL, gd=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    Regression object, that includes methods for \code{coef} and
    \code{vcov} and either (i) a \code{coefficients} component, (ii) a
    \code{coef} component or (iii) a \code{coef<-} method (checked in that order). 
}
  \item{fun}{
    Function that has a signature \code{function(object, ...)}
}
  \item{newdata}{
    A list or data-frame that is passed to \code{fun}. If
    \code{newdata} is \code{NULL} and \code{object} has a non-NULL
    \code{data} component, then \code{newdata <- object$data}. 
}
  \item{gd}{
    specified gradients (default=NULL)
}
  \item{\dots}{
    Other arguments to \code{fun}.
}
}
\details{
  The work is done by \code{rstpm:::numDeltaMethod}.

  One potential issue for some regression objects is that predictions on
  the fitted data may use values from the regression object, so that the
  calculated standard errors are zero. The default work-around
  provided here is define \code{newdata} from \code{object$data}; other
  work-arounds include (i) always passing the original data to
  \code{newdata} and (ii) define a prediction function \code{fun} that
  always uses the original data.
}
\value{
  Returns a data-frame with components \code{Estimate} for the point
  estimate and \code{SE} for the standard errors.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, fun, newdata = NULL, ...) 
{
    if (is.null(newdata) && !is.null(object$data)) 
        newdata <- object$data
    localf <- function(coef, ...) {
        if ("coefficients" \%in\% names(object)) {
            object$coefficients <- coef
        } else if ("coef" \%in\% names(object)) {
            object$coef <- coef
        } else coef(object) <- coef
        fun(object, ...)
    }
    numDeltaMethod(object, localf, newdata = newdata, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
