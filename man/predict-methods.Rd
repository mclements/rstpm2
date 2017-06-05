\name{predict-methods}
\docType{methods}
\alias{predict-methods}
\alias{predict,stpm2-method}
\alias{predict,pstpm2-method}
\title{ Predicted values for an stpm2 or pstpm2 fit}
\description{
  Given an \code{stpm2} fit and an optional list of new data, return predictions
}
\section{Methods}{
  \describe{
    \item{object= "stpm2"}{an \code{stpm2} fit}
}}

\usage{

    \S4method{predict}{stpm2}(object, newdata=NULL,
               type=c("surv","cumhaz","hazard","density","hr","sdiff",
               "hdiff","loghazard","link","meansurv","meansurvdiff",
               "odds","or","margsurv","marghaz","marghr","meanhaz","af",
               "fail","margfail","meanmargsurv","uncured"),
               grid=FALSE,seqLength=300,
               se.fit=FALSE,link=NULL,exposed=incrVar(var),var,
               keep.attributes=TRUE, use.gr=TRUE,...)
    \S4method{predict}{pstpm2}(object, newdata=NULL,
               type=c("surv","cumhaz","hazard","density","hr","sdiff",
               "hdiff","loghazard","link","meansurv","meansurvdiff",
               "odds","or","margsurv","marghaz","marghr","meanhaz","af",
               "fail","margfail","meanmargsurv"),
               grid=FALSE,seqLength=300,
               se.fit=FALSE,link=NULL,exposed=incrVar(var),var,
               keep.attributes=TRUE, use.gr=TRUE,...)
}
\arguments{
  \item{object}{an \code{stpm2} or \code{pstpm2} object}
  \item{newdata}{optional list of new data (required if type in
               ("hr","sdiff","hdiff","meansurvdiff","or","uncured")). For type in
               ("hr","sdiff","hdiff","meansurvdiff","or","af","uncured"), this defines the unexposed
               newdata. This can be combined with \code{grid} to get a
               regular set of event times (i.e. newdata would \emph{not}
               include the event times). }
  \item{type}{specify the type of prediction:
    \itemize{
      \item{"surv"}{survival probabilities}
      \item{"cumhaz"}{cumulative hazard}
      \item{"hazard"}{hazard}
      \item{"density"}{density}
      \item{"hr"}{hazard ratio}
      \item{"sdiff"}{survival difference}
      \item{"hdiff"}{hazard difference}
      \item{"loghazard"}{log hazards}
      \item{"meansurv"}{mean survival}
      \item{"meansurvdiff"}{mean survival difference}
      \item{"odds"}{odds}
      \item{"or"}{odds ratio}
      \item{"margsurv"}{marginal (population) survival}
      \item{"marghaz"}{marginal (population) hazard}
      \item{"marghr"}{marginal (population) hazard ratio}
      \item{"meanhaz"}{mean hazard}
      \item{"af"}{attributable fraction}
      \item{"fail"}{failure (=1-survival)}
      \item{"margfail"}{marginal failure (=1-marginal survival)}
      \item{"meanmargsurv"}{mean marginal survival, averaged over the
               frailty distribution}
      \item{"uncured"}{distribution for the uncured}
    }
  }
  \item{grid}{whether to merge newdata with a regular sequence of event
    times (default=FALSE)}
  \item{seqLength}{length of the sequence used when \code{grid=TRUE}}
  \item{se.fit}{whether to calculate confidence intervals (default=FALSE)}
  \item{link}{allows a different link for the confidence interval
    calculation (default=NULL, such that
    switch(type,surv="cloglog",cumhaz="log",hazard="log",hr="log",sdiff="I",
    hdiff="I",loghazard="I",link="I",odds="log",or="log",margsurv="cloglog",
    marghaz="log",marghr="log"))}
  \item{exposed}{a function that takes newdata and returns a transformed
               data-frame for those exposed or the counterfactual (defaults to incrementing ``var'')}
  \item{var}{specify the variable name or names for the
               exposed/unexposed (names are given as characters)}
  \item{keep.attributes}{Boolean to determine whether the output should
               include the newdata as an attribute (default=TRUE)}
  \item{use.gr}{Boolean to determine whether to use gradients in the
               variance calculations when they are available (default=TRUE)}
  \item{\dots}{additional arguments (for generic compatibility)}
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
