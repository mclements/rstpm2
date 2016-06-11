\name{tvcCoxph-class}
\Rdversion{1.1}
\docType{class}
\alias{tvcCoxph-class}
\alias{plot,tvcCoxph,missing-method}

\title{Class \code{"tvcCoxph"}}
\description{
Experimental approach to modelling time-dependent effects in Cox regression.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("tvcCoxph", ...)}
  or \code{cox.tvc(...)}. See the \code{"\linkS4class{mle2}"}
  documentation.
}
\section{Slots}{
  \describe{
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
  }
}
\section{Extends}{
Class \code{"\linkS4class{mle2}"}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "tvcCoxph", y = "missing")}: ... }
	 }
}
%% \references{
%% %%  ~~put references to the literature/web site here~~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%% %%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%% }
\examples{
showClass("tvcCoxph")
}
\keyword{classes}
