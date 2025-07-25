\name{predictnl}
\alias{predictnl}
\alias{predictnl.default}
\alias{predictnl.lm}
\alias{predict.formula}
\alias{confint.predictnl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of standard errors using the numerical delta method.
}
\description{
A simple, yet exceedingly useful, approach to estimate the variance of a
function using the numerical delta method. A number of packages provide
functions that analytically calculate the gradients; we use numerical
derivatives, which generalises to models that do not offer analytical
derivatives (e.g. ordinary differential equations, integration), or to
examples that are tedious or error-prone to calculate (e.g. sums of
predictions from GLMs).
}
\usage{
\method{predictnl}{default}(object, fun, newdata=NULL, gd=NULL, ...)
\method{predictnl}{lm}(object, fun, newdata=NULL, ...)
\method{predict}{formula}(object,data,newdata,na.action,type="model.matrix",...)
\method{confint}{predictnl}(object, parm, level=0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object with \code{coef}, \code{vcov} and \code{`coef<-`}
methods (required). 
}
  \item{fun}{
A function that takes \code{object} as the first argument, possibly with
\code{newdata} and other arguments (required). See notes for why it is
often useful to include \code{newdata} as an argument to the function.
}
  \item{newdata}{
An optional argument that defines newdata to be passed to \code{fun}.
}
  \item{gd}{
An optional matrix of gradients. If this is not specified, then the
gradients are calculated using finite differences.
}
  \item{parm}{
    currently ignored
}
  \item{level}{
    significance level for 2-sided confidence intervals
}
  \item{data}{
    object used to define the model frame
}
  \item{na.action}{
    passed to \code{model.frame}
}
  \item{type}{
    currently restricted to \code{"model.matrix"}
}
  \item{\dots}{
    Other arguments that are passed to \code{fun}.
}
}
\details{
The signature for \code{fun}
is either \code{fun(object, ...)} or \code{fun(object, newdata=NULL,
  ...)}.

The different \code{predictnl} methods call the utility function
\code{numDeltaMethod}, which in turn calls the \code{grad} function for
numerical differentiation. The \code{numDeltaMethod} function calls the
standard \code{coef} and \code{vcov} methods, and the non-standard
\code{`coef<-`} method for changing the coefficients in a regression
object. This non-standard method has been provided for several
regression objects and essentially mirrors the \code{coef} method.

One potential issue is that some \code{predict} methods do not
re-calculate their predictions for the fitted dataset (i.e. when
\code{newdata=NULL}). As the \code{predictnl} function changes the
fitted coefficients, it is required that the predictions are
re-calculated. One solution is to pass \code{newdata} as an argument to
both \code{predictnl} and \code{fun}; alternatively, \code{newdata} can
be specified in \code{fun}. These approaches are described in the examples
below. The \code{numDeltaMethod} method called by \code{predictnl}
provides a warning when the variance estimates are zero, which may be
due to this cause.

For completeness, it is worth discussing why the example
\code{predictnl(fit,predict)} does not work for when \code{fit} is a
\code{glm} object. First, \code{predict.glm} does not update the
predictions for the fitted data. Second, the default \code{predict}
method has a signature \code{predict(object, ...)}, which does not
include a \code{newdata} argument. We could then either (i) require that
a \code{newdata} argument be passed to the \code{fun} function for all
examples, which would make this corner case work, or (ii) only pass the
\code{newdata} argument if it is non-null or in the formals for the
\code{fun} function, which would fail for this corner case. The current
API defaults to the latter case (ii).  To support this approach, the
\code{predictnl.lm} method replaces a null \code{newdata} with
\code{object$data}. We also provide a revised
\code{numdelta:::predict.lm} method that performs the same operation,
although its use is not encouraged due to its clumsiness.
}
\value{ Returns an object of class
an object with class \code{c("predictnl","data.frame")} elements
\code{c("fit","se.fit","Estimate","SE")} and with methods \code{print}
and \code{confint}. Note that the Estimate and SE fields are deprecated
and their use is discouraged, as we would like to remove them from future releases.
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{
Mark Clements
}
%% \note{
%% %%  ~~further notes~~
%% }

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{

df <- data.frame(x=0:1, y=c(10, 20))
fit <- glm(y ~ x, df, family=poisson)

predictnl(fit,
          function(obj,newdata)
          diff(predict(obj,newdata,type="response")))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
