\name{stpm2}
\alias{stpm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fully parametric generalised survival model
}
\description{
  This implements the generalised survival model g(S(t|x)) = eta, where
g is a link function, S is survival, t is time, x are covariates and eta
is a linear predictor. The main model assumption is that the time
effects in the linear predictor are smooth. This extends the class of
flexible parametric survival models developed by Royston and
colleagues. The model has been extended to include relative survival,
Gamma frailties and normal random effects.
}
\usage{
stpm2(formula, data, 
      smooth.formula = NULL, smooth.args = NULL,
      df = 3, cure = FALSE, logH.args = NULL,
      logH.formula = NULL, tvc = NULL, tvc.formula =
      NULL, control = list(parscale = 1, maxit = 300),
      init = NULL, coxph.strata = NULL, weights = NULL,
      robust = FALSE, baseoff = FALSE, bhazard = NULL,
      timeVar = "", time0Var = "", use.gr = TRUE,
      optimiser=c("BFGS","NelderMead"),
      reltol=1.0e-8, trace = 0,
      link.type=c("PH","PO","probit","AH","AO"), theta.AO=0, 
      frailty = !is.null(cluster) & !robust, cluster = NULL, logtheta=-6, nodes=9,
      RandDist=c("Gamma","LogN"), recurrent = FALSE, adaptive=TRUE,
      maxkappa=1000, Z=~1,
      contrasts = NULL,
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
  \item{smooth.formula}{
a formula for describing the time effects for the linear predictor,
including the baseline and the time-dependent effects
(default=NULL). Only one of \code{df}, \code{smooth.formula},
\code{smooth.args}, \code{logH.args} or
\code{logH.formula} is required. The default model is equal to \code{nsx(log(time),df=3)}.
}
  \item{smooth.args}{
a list describing the arguments for the \code{nsx} function for modelling
the baseline time effect on the linear predictor scale (default=NULL). Use this or
smooth.formula for
changing the knot placement and specifying cure models.
}
  \item{df}{
an integer that describes the degrees of freedom for the \code{ns}
function for modelling the baseline log-cumulative hazard (default=3). 
}
  \item{logH.args}{
as per \code{smooth.args}. Deprecated.
}
\item{logH.formula}{
  as per \code{smooth.formula}. Deprecated.
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
\item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}.
}
\item{subset}{an optional vector specifying a subset of observations to
  be used in the fitting process.}
\item{cure}{logical for whether to estimate a cure model.}
\item{time0Var}{string variable to determine the entry variable; useful
  for when more than one data variable is used in the entry time.}
\item{use.gr}{logical indicating whether to use gradients in the calculation}
\item{optimiser}{select which optimiser is used}
\item{link.type}{type of link function. For "PH" (generalised
  proportional hazards), g(S)=log(-log(S)); for
"PO" (generalised proportional odds), g(S)=-logit(S); for "probit"
(generalised probit), g(S)=-probit(S); for "AH" (generalised additive
hazards), g(S)=-log(S); for "AO" (generalised Aranda-Ordaz), g(S)=log((S^(-theta.AO)-1)/theta.AO).}
\item{theta.AO}{theta parameter for the Aranda-Ordaz link type.}
\item{reltol}{relative tolerance for the model convergence}
\item{trace}{logical for whether to provide trace information}
\item{frailty}{logical for whether to fit a shared frailty model}
\item{cluster}{string for the data variable that determines the cluster
  for the frailty}
\item{nodes}{number of integration points for Gaussian quadrature}
\item{RandDist}{type of distribution for the random effect or frailty}
\item{recurrent}{logical for whether clustered, left truncated data are
  recurrent or for first event (where the latter requires an adjustment
  for the frailties or random effects)}
\item{logtheta}{initial value for log-theta used in the gamma shared frailty
  model}
\item{adaptive}{logical for whether to use adaptive or non-adaptive
  quadrature}
\item{maxkappa}{double float value for the maximum value of the weight
  used in the constraint}
\item{Z}{formula for the design matrix for the random effects}
\item{\dots}{
  additional arguments to be passed to the \code{\link{mle2}} .
}
}
\details{
  
  The implementation extends the \code{mle2} object from the
  \code{bbmle} package. The model inherits all of the methods from the
  \code{mle2} class.

  The default linear predictor includes a time effect modelled using
  natural splines for log(time) with three degrees of freedom.
  
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
\author{
 Mark Clements, Xing-Rong Liu.
}
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
