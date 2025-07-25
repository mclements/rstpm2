\name{predictnl-methods}
\docType{methods}
\alias{predictnl-methods}
\alias{predictnl,mle2-method}
\alias{predictnl,stpm2-method}
\alias{predictnl,pstpm2-method}
\alias{predictnl,aft-method}
\title{Estimation of standard errors using the numerical delta method.}
\description{
A simple, yet exceedingly useful, approach to estimate the variance of a
function using the numerical delta method. A number of packages provide
functions that analytically calculate the gradients; we use numerical
derivatives, which generalises to models that do not offer analytical
derivatives (e.g. ordinary differential equations, integration), or to
examples that are tedious or error-prone to calculate (e.g. sums of
predictions from GLMs).
}
\section{Methods}{
  \describe{
    \item{object= "stpm2"}{an \code{stpm2} fit}
}}

\usage{

\S4method{predictnl}{mle2}(object, fun, newdata=NULL, gd=NULL, ...)
\S4method{predictnl}{stpm2}(object, fun, newdata=NULL,
    link = c("I", "log", "cloglog", "logit"), gd=NULL, ...)
\S4method{predictnl}{pstpm2}(object, fun, newdata=NULL,
    link = c("I", "log", "cloglog", "logit"), gd=NULL, ...)
\S4method{predictnl}{aft}(object, fun, newdata=NULL,
    link = c("I", "log", "cloglog", "logit"), gd=NULL, ...)

}
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
  \item{link}{
    A character string to represent the link on which to calculate the
    variance and confidence intervals.
}
  \item{gd}{
An optional matrix of gradients. If this is not specified, then the
gradients are calculated using finite differences.
}
  \item{\dots}{
    Other arguments that are passed to \code{fun}.
}

}

\value{
  A data-frame with components \code{Estimate}, \code{lower} and
  \code{upper}, with an attribute "newdata" for the \code{newdata}
  data-frame.  
}

\details{
  The confidence interval estimation is based on the delta method using
  numerical differentiation.
}

\seealso{\code{\link{stpm2}}}
\keyword{methods}
%%\keyword{ ~~ other possible keyword(s)}
