\name{stpm2}
\alias{stpm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Flexible parametric survival model.
}
\description{
This implements the Royston-Parmar model.
}
\usage{
stpm2(formula, data, df = 3, cure = FALSE, logH.args = NULL,
      logH.formula = NULL, tvc = NULL, tvc.formula =
      NULL, control = list(parscale = 0.1, maxit = 300),
      init = NULL, coxph.strata = NULL, weights = NULL,
      robust = FALSE, baseoff = FALSE, bhazard = NULL,
      timeVar = "", time0Var = "", use.gr = TRUE,
      use.rcpp= TRUE, reltol=1.0e-8, trace = 0,
      type=c("PH","PO","probit","AH"), frailty = FALSE,
      cluster = NULL, logtheta=0, contrasts = NULL,
      subset = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object, with the response on the left of a \code{~} operator, and
the regression terms (excluding time) on the right.  The response must be a survival object as
returned by the \code{\link{Surv}} function. The terms should include linear
terms for any time-varying coefficients. [required]
}
  \item{data}{
a data.frame in which to interpret the variables named in
the \code{formula} argument. [at present: required]
}
  \item{df}{
an integer that describes the degrees of freedom for the \code{ns}
function for modelling the baseline log-cumulative hazard (default=3). 
}
  \item{logH.args}{
a list describing the arguments for the \code{nsx} function for modelling
the baseline log-cumulative hazard (default=NULL). Use this or
logH.formula for
changing the knot placement and specifying cure models.
}
  \item{logH.formula}{
a formula for describing the baseline log-cumulative hazard function
(default=NULL). Only one of \code{df}, \code{logH.args} or
\code{logH.formula} is required. The default model is equal to \code{nsx(log(time),df=3)}.
}
  \item{tvc}{
a list with the names of the time-varying coefficients and the degrees
of freedom (e.g. \code{tvc=list(x=3)} specifies \code{x} as a time-varying
coefficient with 3 degrees of freedom).
}
  \item{tvc.formula}{
a formula for describing the time-varying coefficients. If a
time-varying coefficient is being model, then only one of
\code{tvc} and \code{tvc.formula} is required.
}
  \item{bhazard}{
a vector for the background hazard for relative survival estimation. At
present, this does not use \code{data} and it is required for
all individuals - although it is only used at the event times.
}
\item{control}{\code{control} argument passed to \code{optim}.}
\item{init}{\code{init} should either be \code{FALSE}, such that initial
                  values will be determined using Cox regression, or a numeric
                  vector of initial values.}
\item{coxph.strata}{variable in the \code{data} argument for
  stratification of the \code{coxph} model fit for estimating initial values.}
\item{weights}{an optional vector of 'prior weights' to be used in the
  fitting process. Should be \code{NULL} or a numeric vector.}
\item{robust}{Boolean used to determine whether to use a robust variance
  estimator.}
\item{baseoff}{Boolean used to determine whether fully define the model
using \code{tvc.formula} rather than combining \code{logH.formula} and \code{tvc.formula}}
\item{timeVar}{variable defining the time variable. By default, this is
determined from the survival object, however this may be ambiguous if
two variables define the time}
\item{use.gr}{Boolean to determine whether to use the gradient in the optimisation}
\item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}.
}
\item{subset}{an optional vector specifying a subset of observations to
  be used in the fitting process.}
\item{cure}{logical for whether to estimate a cure model.}
\item{time0Var}{string variable to determine the entry variable; useful
  for when more than one data variable is used in the entry time.}
\item{type}{type of link function.}
\item{use.rcpp}{logical for whether to use C++/Rcpp for estimation}
\item{reltol}{relative tolerance for the model convergence}
\item{trace}{logical for whether to provide trace information}
\item{frailty}{logical for whether to fit a shared frailty model (experimental)}
\item{cluster}{string for the data variable that determines the cluster
  for the frailty (experimental)}
\item{logtheta}{initial value for log-theta used in the gamma shared frailty
  model (experimental)}
\item{\dots}{
  additional arguments to be passed to the \code{\link{mle2}} .
}
}
\details{
  
  
  The implementation extends the \code{mle2} object from the
  \code{bbmle} package.

  Note that the default baseline log-cumulative hazard function uses
  natural splines for log of time with three degrees of freedom. The
  design matrix is calculated for the event times rather than all
  observations.
  
}
\value{
An \code{stpm2-class} object that inherits from \code{mle2-class}. 
  %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
data(brcancer)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3))

## some predictions
head(predict(fit,se.fit=TRUE,type="surv"))
head(predict(fit,se.fit=TRUE,type="hazard"))

## some plots
plot(fit,newdata=data.frame(hormon=0),type="hazard")
plot(fit,newdata=data.frame(hormon=0),type="surv")

## the same model using logH.formula
summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,logH.formula=~ns(log(rectime),df=3)))


## time-varying coefficient
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                     tvc=list(hormon=3)))
anova(fit,fit.tvc) # compare with and without tvc

## some more plots
plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon", ylim=c(0,2))
                                        # no lines method: use add=TRUE
plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon",
     add=TRUE,ci=FALSE,line.col=2)

plot(fit.tvc,newdata=data.frame(hormon=0),type="sdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hazard")
plot(fit.tvc,newdata=data.frame(hormon=1),type="hazard",line.col=2,ci=FALSE,add=TRUE)

## compare number of knots
hormon0 <- data.frame(hormon=0)
plot(fit,type="hazard",newdata=hormon0)
AIC(fit)
for (df in 4:6) {
    fit.new <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=df)
    plot(fit.new,type="hazard",newdata=hormon0,add=TRUE,ci=FALSE,line.col=df)
    print(AIC(fit.new))
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
