\name{residuals-methods}
\docType{methods}
\alias{residuals-methods}
\alias{residuals,stpm2-method}
\alias{residuals,pstpm2-method}
\title{Residual values for an stpm2 or pstpm2 fit}
\description{
  Given an \code{stpm2} or \code{pstpm2} fit, return residuals
}
\section{Methods}{
  \describe{
    \item{object= "stpm2"}{an \code{stpm2} fit}
}}

\usage{

    \S4method{residuals}{stpm2}(object, type=c("li","gradli"))
    \S4method{residuals}{pstpm2}(object, type=c("li","gradli"))
}
\arguments{
  \item{object}{an \code{stpm2} or \code{pstpm2} object}
  \item{type}{specify the type of residuals:
    \itemize{
      \item{"li"}{log-likelihood components (not strictly residuals)}
      \item{"gradli"}{gradient of the log-likelihood components (not
    strictly residuals)}
    }
  }
}

\value{
  A vector or matrix.
}

\details{
  The gradients are analytical.
}

\seealso{\code{\link{stpm2}}}
\keyword{methods}
%%\keyword{ ~~ other possible keyword(s)}
