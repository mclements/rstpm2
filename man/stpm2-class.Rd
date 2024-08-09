\name{stpm2-class}
\Rdversion{1.1}
\docType{class}
\alias{stpm2-class}
\alias{plot,stpm2,missing-method}
\alias{lines,stpm2-method}
%\alias{predict,stpm2-method}
\alias{predictnl,stpm2-method}
\alias{summary,stpm2-method}
\alias{eform,stpm2-method}

\title{Class "stpm2" ~~~ }
\description{Regression object for \code{stpm2}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stpm2", ...)} and
\code{stpm2( ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{xlevels}:}{Object of class \code{"list"} ~~ }
    \item{\code{contrasts}:}{Object of class \code{"listOrNULL"} ~~ }
    \item{\code{terms}:}{Object of class \code{"terms"} ~~ }
    \item{\code{logli}:}{Object of class \code{"function"} ~~ }
    \item{\code{lm}:}{Object of class \code{"lm"} ~~ }
    \item{\code{timeVar}:}{Object of class \code{"character"} ~~ }
    \item{\code{time0Var}:}{Object of class \code{"character"} ~~ }
    \item{\code{timeExpr}:}{Object of class \code{"nameOrcall"} ~~ }
    \item{\code{time0Expr}:}{Object of class \code{"nameOrcall"} ~~ }
    \item{\code{delayed}:}{Object of class \code{"logical"} ~~ }
    \item{\code{frailty}:}{Object of class \code{"logical"} ~~ }
    \item{\code{interval}:}{Object of class \code{"logical"} ~~ }
    \item{\code{model.frame}:}{Object of class \code{"list"} ~~ }
    \item{\code{call.formula}:}{Object of class \code{"formula"} ~~ }
    \item{\code{x}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{xd}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{termsd}:}{Object of class \code{"terms"} ~~ }
    \item{\code{Call}:}{Object of class \code{"character"} ~~ }
    \item{\code{y}:}{Object of class \code{"Surv"} ~~ }
    \item{\code{link}:}{Object of class \code{"list"} ~~ }
    \item{\code{call}:}{Object of class \code{"language"} ~~ }
    \item{\code{call.orig}:}{Object of class \code{"language"} ~~ }
    \item{\code{coef}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{fullcoef}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{min}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{details}:}{Object of class \code{"list"} ~~ }
    \item{\code{minuslogl}:}{Object of class \code{"function"} ~~ }
    \item{\code{method}:}{Object of class \code{"character"} ~~ }
    \item{\code{data}:}{Object of class \code{"list"} ~~ }
    \item{\code{formula}:}{Object of class \code{"character"} ~~ }
    \item{\code{optimizer}:}{Object of class \code{"character"} ~~ }
    \item{\code{args}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{\link[bbmle]{mle2}}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "stpm2", y = "missing")}: ... }
    \item{lines}{\code{signature(x = "stpm2", ...)}: ... }
    %\item{predict}{\code{signature(object = "stpm2")}: ... }
    \item{predictnl}{\code{signature(object = "stpm2", ...)}: ... }
    \item{summary}{\code{signature(object = "stpm2", ...)}: ... }
    \item{eform}{\code{signature(object = "stpm2", ...)}: ... }
	 }
}
%% \references{ ~put references to the literature/web site here ~ }
%% \author{ ~~who you are~~ }
%% \note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% 	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%% 	or \code{\linkS4class{CLASSNAME}} for links to other classes
%% }
\examples{
showClass("stpm2")
}
\keyword{classes}
