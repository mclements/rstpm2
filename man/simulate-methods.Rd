\name{simulate-methods}
\docType{methods}
\alias{simulate-methods}
\alias{simulate,stpm2-method}
\alias{simulate,pstpm2-method}
\title{Simulate values from an stpm2 or pstpm2 fit}
\description{
  Given an \code{stpm2} fit and a data-frame
  of new data, return simulated values
}
\section{Methods}{
\describe{
  \item{object = "stpm2"}{an \code{stpm2} fit}
}}
\usage{
    \S4method{simulate}{stpm2}(object, nsim=1,
        seed=NULL, newdata=NULL,
        lower=1e-06, upper=1e+05, start=NULL, ...)
    \S4method{simulate}{pstpm2}(object, nsim=1,
        seed=NULL, newdata=NULL,
        lower=1e-06, upper=1e+05, start=NULL, ...)
}
\arguments{
  \item{object}{an stpm2 or pstpm2 object}
  \item{nsim}{number of simulations per row in newdata}
  \item{seed}{optional random number seed}
  \item{newdata}{list of new data. If not specified, then defaults to object@data}
  \item{lower}{smallest possible time}
  \item{upper}{largest possible time}
  \item{start}{left truncated entry time (assumed to be zero if NULL)}
  \item{\dots}{additional arguments (for generic compatibility)}
}
\examples{
set.seed(1002)
fit1 <- gsm(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3)
simulate(fit1, nsim=10, newdata=data.frame(hormon=1))
simulate(fit1, newdata=data.frame(hormon=0:1))

}
\keyword{methods}
