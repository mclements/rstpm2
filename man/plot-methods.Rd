\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,stpm2-method}
\alias{plot,pstpm2-method}
\title{ plots for an stpm2 fit}
\description{
Given an \code{stpm2} fit, return a plot
}
\section{Methods}{
\describe{
\item{x = "stpm2", y = "missing"}{an \code{stpm2} fit}
}}

\usage{

    \S4method{plot}{stpm2}(x,y,newdata,type="surv",
                      xlab="Time",line.col=1,ci.col="grey",
                      add=FALSE,ci=TRUE,rug=TRUE,
                      var=NULL,exposed=incrVar(var),times=NULL,...)
    \S4method{plot}{pstpm2}(x,y,newdata,type="surv",
                      xlab="Time",line.col=1,ci.col="grey",
                      add=FALSE,ci=TRUE,rug=TRUE,
                      var=NULL,exposed=incrVar(var),times=NULL,...)
}
\arguments{
  \item{x}{an \code{stpm2} object}
  \item{y}{not used (for generic compatibility)}
  \item{newdata}{required list of new data. This defines the unexposed
               newdata (\emph{excluding} the event times).}
  \item{type}{specify the type of prediction
  }
  \item{xlab}{x-axis label}
  \item{line.col}{line colour}
  \item{ci.col}{confidence interval colour}
  \item{ci}{whether to plot the confidence interval band (default=TRUE)}
  \item{add}{whether to add to the current plot (\code{add=TRUE}) or make a new
                      plot (\code{add=FALSE}) (default=FALSE)}
  \item{rug}{whether to add a rug plot of the event times to the current
                      plot (default=TRUE)}
  \item{var}{specify the variable name or names for the
    exposed/unexposed (names are given as characters)}
  \item{exposed}{function that takes newdata and returns the exposed
    dataset. By default, this increments \code{var}}
  \item{times}{specifies the times. By default, this uses a span of the
                      observed times.}
  \item{\dots}{additional arguments (add to the \code{plot} command)}
}

%% \examples{


%%   }
\seealso{\code{\link{stpm2}}}
\keyword{methods}
%%\keyword{ ~~ other possible keyword(s)}


