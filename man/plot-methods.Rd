\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,stpm2-method}
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
                      var=NULL,...)
}
\arguments{
  \item{x}{an \code{stpm2} object}
  \item{y}{not used (for generic compatibility)}
  \item{newdata}{required list of new data. This defines the unexposed
               newdata (\emph{excluding} the event times).}
  \item{type}{specify the type of prediction:
    \itemize{
      \item{"surv"}{survival probabilities}
      \item{"cumhaz"}{cumulative hazard}
      \item{"hazard"}{hazard}
      \item{"hr"}{hazard ratio}
      \item{"sdiff"}{survival difference}
      \item{"hdiff"}{hazard difference}
    }
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
  \item{\dots}{additional arguments (add to the \code{plot} command)}
}

\examples{


## The function is currently defined as
setMethod("plot", signature(x="stpm2", y="missing"),
          function(x,y,newdata,type="surv",
                      xlab="Time",line.col=1,ci.col="grey",
                      add=FALSE,ci=TRUE,rug=TRUE,
                      var=NULL,...) {
  y <- predict(x,newdata,type=type,var=var,grid=TRUE,se.fit=TRUE)
  ylab <- switch(type,hr="Hazard ratio",hazard="Hazard",surv="Survival",
                 sdiff="Survival difference",hdiff="Hazard difference")
  xx <- attr(y,"newdata")
  xx <- xx[,ncol(xx)]
  if (!add) matplot(xx, y, type="n", xlab=xlab, ylab=ylab, ...)
  if (ci) polygon(c(xx,rev(xx)), c(y[,2],rev(y[,3])), col=ci.col, border=NA)
  lines(xx,y[,1],col=line.col)
  if (rug) {
      Y <- x@y
      eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
      rug(eventTimes)
    }
})

  }
\seealso{\code{\link{stpm2}}}
\keyword{methods}
%%\keyword{ ~~ other possible keyword(s)}


