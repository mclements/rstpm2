\name{aft-class}
\Rdversion{1.1}
\docType{class}
\alias{aft-class}
\alias{plot,aft,missing-method}
\alias{lines,aft-method}
\alias{predict,aft-method}
%\alias{predictnl,aft-method}
%\alias{summary,aft-method}

\title{Class "stpm2" ~~~ }
\description{Regression object for \code{aft}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("aft", ...)} and
\code{aft( ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{args}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class for \code{\link[bbmle]{mle2}}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "aft", y = "missing")}: ... }
    \item{lines}{\code{signature(x = "aft")}: ... }
    \item{predict}{\code{signature(object = "aft")}: ... }
    %\item{predictnl}{\code{signature(object = "aft")}: ... }
    %\item{summary}{\code{signature(object = "aft", ...)}: ... }
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
showClass("aft")
}
\keyword{classes}
