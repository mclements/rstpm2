\name{pstpm2-class}
\Rdversion{1.1}
\docType{class}
\alias{pstpm2-class}
\alias{plot,pstpm2,missing-method}
\alias{lines,pstpm2-method}
%\alias{predict,pstpm2-method}
\alias{predictnl,pstpm2-method}
\alias{AIC,pstpm2-method}
\alias{AICc,pstpm2-method}
\alias{BIC,pstpm2-method}
\alias{anova,pstpm2-method}
\alias{qAICc,pstpm2-method}
\alias{summary,pstpm2-method}
\alias{eform,pstpm2-method}

\title{Class "pstpm2" }
\description{Regression object for \code{pstpm2}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pstpm2", ...)} and
\code{pstpm2( ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{xlevels}:}{Object of class \code{"list"} ~~ }
    \item{\code{contrasts}:}{Object of class \code{"listOrNULL"} ~~ }
    \item{\code{terms}:}{Object of class \code{"terms"} ~~ }
    \item{\code{gam}:}{Object of class \code{"gam"} ~~ }
    \item{\code{logli}:}{Object of class \code{"function"} ~~ }
    \item{\code{timeVar}:}{Object of class \code{"character"} ~~ }
    \item{\code{time0Var}:}{Object of class \code{"character"} ~~ }
    \item{\code{time0Expr}:}{Object of class \code{"nameOrcall"} ~~ }
    \item{\code{timeExpr}:}{Object of class \code{"nameOrcall"} ~~ }
    \item{\code{like}:}{Object of class \code{"function"} ~~ }
    \item{\code{model.frame}:}{Object of class \code{"list"} ~~ }
    \item{\code{fullformula}:}{Object of class \code{"formula"} ~~ }
    \item{\code{delayed}:}{Object of class \code{"logical"} ~~ }
    \item{\code{frailty}:}{Object of class \code{"logical"} ~~ }
    \item{\code{x}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{xd}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{termsd}:}{Object of class \code{"terms"} ~~ }
    \item{\code{Call}:}{Object of class \code{"character"} ~~ }
    \item{\code{y}:}{Object of class \code{"Surv"} ~~ }
    \item{\code{sp}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nevent}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{link}:}{Object of class \code{"list"} ~~ }
    \item{\code{edf}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{edf_var}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{df}:}{Object of class \code{"numeric"} ~~ }
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
Class \code{"\linkS4class{mle2}"}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "pstpm2", y = "missing")}: ... }
    \item{lines}{\code{signature(x = "pstpm2", ...)}: ... }
    %\item{predict}{\code{signature(object = "pstpm2")}: ... }
    \item{anova}{\code{signature(object = "pstpm2",...)}: ... }
    \item{AIC}{\code{signature(object = "pstpm2",...,k=2)}: ... }
    \item{AICc}{\code{signature(object = "pstpm2",...,nobs=NULL, k=2)}: ... }
    \item{BIC}{\code{signature(object = "pstpm2",..., nobs = NULL)}: ... }
    \item{qAICc}{\code{signature(object = "pstpm2",..., nobs = NULL, dispersion = 1, k = 2)}: ... }
    \item{qAIC}{\code{signature(object = "pstpm2",..., dispersion = 1, k = 2)}: ... }
    \item{summary}{\code{signature(object = "pstpm2",...)}: ... }
    \item{eform}{\code{signature(object = "pstpm2",...)}: ... }
    \item{predictnl}{\code{signature(object = "pstpm2",...)}: ... }
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
showClass("pstpm2")
}
\keyword{classes}
