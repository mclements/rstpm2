\name{Rstpm2-package}
\Rdversion{1.1}
\alias{Rstpm2-package}
\alias{Rstpm2}
\docType{package}
\title{
Flexible parametric survival models.
}
\description{
  The package implements the stpm2 models from Stata. Such models use a
  flexible parametric formulation for survival models, using natural
  splines to model the log-cumulative hazard. Model predictions are
  rich, allowing for direct estimation of the hazard, survival, hazard
  ratios, hazard differences and survival differences. The models allow
  for time-varying effects, left truncation and relative survival.
  
  The R implementation departs from the Stata implementation, using the
  ns() function, which is based on a projection of B-splines, rather
  than using truncated power splines as per Stata.
}

\details{
\tabular{ll}{
Package: \tab Rstpm2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-07-06\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
Depends: \tab methods, bbmle\cr
Imports: \tab splines, survival, stats, graphics\cr
}
The package exports the \code{\link{stpm2}} object, which inherits from
the \code{\link{mle2}} object from the \code{bbmle} package. Methods are
specified for the \code{\link{stpm2}} object, including \code{predict}
and \code{plot} methods. 
}
\author{
Mark Clements and Paul Lambert.

Maintainer: <mark.clements@ki.se>
}
% \references{
% ~~ Literature or other references for background information ~~
% }
% ~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
% ~~ the R documentation directory ~~
\keyword{ package }
\seealso{
\code{\link{stpm2}}
}
\examples{
data(brcancer)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3))
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                     tvc=list(hormon=3)))
anova(fit,fit.tvc)
plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon")
}
