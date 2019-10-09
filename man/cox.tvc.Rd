\name{cox.tvc}
\alias{cox.tvc}
\title{Test for a time-varying effect in the \code{coxph} model}
\description{
  Test for a time-varying effect in the \code{coxph} model by re-fitting
  the partial likelihood including a time-varying effect, plot 
  the effect size, and return the re-fitted model. The
  main advantage of this function over the \code{tt()} special is
  that it scales well for moderate sized datasets
  (cf. \code{tt} which expands the dataset and scales very poorly).
}
\usage{
cox.tvc(obj, var=NULL, method="logt")
}
\arguments{
  \item{obj}{A \code{coxph} object. Currently restricted to right
    censoring with Breslow ties and without stratification, etc.}
  \item{var}{String for the effect name. Currently assumes simple
    continuous effects.}
  \item{method}{A string representing the possible time
    transformations. Currently only "logt".}
}
\value{
  Returns a \code{tvcCoxph} object (which inherits from the \code{mle2} class) of the re-fitted model.
}
\seealso{
  \code{\link{coxph}}, \code{\link{cox.zph}}
}
%% \references{
%% }
\examples{
## As per the example for cox.zph:
fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps,  
             data=ovarian) 
temp <- rstpm2:::cox.tvc(fit, "age") 
print(temp)                  # display the results 
plot(temp)                   # plot curves 
}
\keyword{time-varying}
\keyword{Cox}
