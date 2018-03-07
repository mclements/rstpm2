\name{pstpm2}
\Rdversion{1.1}
\alias{pstpm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Penalised generalised survival model
}
\description{
  This implements the generalised survival model g(S(t|x)) = eta, where
g is a link function, S is survival, t is time, x are covariates and eta
is a linear predictor. The linear predictor can include penalised
smoothers for the time effects, for time:covariate interactions and for
covariate effects using the mgcv smoothers. The main model assumption is that the time
effects in the linear predictor are smooth. This extends the class of
flexible parametric survival models developed by Royston and
colleagues. The model has been extended to include relative survival,
Gamma frailties and normal random effects.
}
\usage{
pstpm2(formula, data, smooth.formula = NULL, smooth.args = NULL,
       logH.args = NULL, 
       tvc = NULL, 
       control = list(parscale = 1, maxit = 300), init = NULL,
       coxph.strata = NULL, coxph.formula = NULL,
       weights = NULL, robust = FALSE, 
       bhazard = NULL, bhazinit = 0.1, timeVar = "", time0Var = "",
       sp=NULL, use.gr = TRUE, 
       criterion=c("GCV","BIC"), penalty = c("logH","h"),
       smoother.parameters = NULL,
       alpha=if (is.null(sp)) switch(criterion,GCV=1,BIC=1) else 1,
       sp.init=1, trace = 0,
       link.type=c("PH","PO","probit","AH","AO"), theta.AO=0,
       optimiser = c("BFGS", "NelderMead", "Nlm"), log.time.transform=TRUE,
       recurrent = FALSE, frailty=!is.null(cluster) & !robust,cluster = NULL,
       logtheta=-6, nodes=9,
       RandDist=c("Gamma","LogN"), adaptive = TRUE, maxkappa=1000, Z = ~1, 
       reltol = list(search = 1.0e-10, final = 1.0e-10, outer=1.0e-5),outer_optim=1,
       contrasts = NULL, subset = NULL, robust_initial = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object, with the response on the left of a \code{~} operator, and
the parametric terms on the right.  The response must be a survival object as
returned by the \code{\link{Surv}} function. [required]
}
  \item{data}{
a data.frame in which to interpret the variables named in
the \code{formula} argument.
}
  \item{smooth.formula}{
a \code{mgcv::gam} formula for describing the time effects and
time-dependent effects and smoothed covariate effects on the linear
predictor scale
(default=NULL). The default model is equal to \code{~s(log(time),k=-1)}
where \code{time} is the time variable.
}
  \item{smooth.args}{
a list describing the arguments for the \code{s} function for modelling
the baseline time effect on the linear predictor scale (default=NULL). 
}
  \item{logH.args}{
as per \code{smooth.args}. Deprecated.
}
  \item{tvc}{
a list with the names of the time-varying coefficients
(e.g. tvc=list(hormon), which is equivalent to smooth.formula=~...+s(log(time),by=hormon)).
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
\item{bhazard}{variable for the baseline hazard for relative survival}
  \item{bhazinit}{
scalar used to adjust the background cumulative hazards for calculating
initial values. Default=0.1.
}
\item{timeVar}{variable defining the time variable. By default, this is
determined from the survival object, however this may be ambiguous if
two variables define the time}
% \item{time0Var}{variable defining the entry time variable. By default, this is
% determined from the survival object, however this may be ambiguous if
% two variables define the entry time}
\item{sp}{fix the value of the smoothing parameters.}
\item{use.gr}{in R, a Boolean to determine whether to use the gradient in the optimisation}
\item{criterion}{in Rcpp, determine whether to use "GCV" or "BIC" for for the smoothing parameter selection.}
\item{penalty}{use either the "logH" penalty, which is the default penalty from mgcv, or the "h" hazard penalty.}
\item{smoother.parameters}{for the hazard penalty, a list with components which are lists with components var, transform and inverse.}
\item{alpha}{an ad hoc tuning parameter for the smoothing parameter.}
\item{sp.init}{initial values for the smoothing parameters.}
\item{trace}{integer for trace reporting; 0 represents no additional reporting.}
\item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}.
}
\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
\item{coxph.formula}{additional formula used to improve the fitting of
  initial values [optional and rarely used].}
\item{time0Var}{string variable to determine the entry variable; useful
  for when more than one data variable is used in the entry time.}
\item{link.type}{type of link function. For "PH" (generalised
  proportional hazards), g(S)=log(-log(S)); for
"PO" (generalised proportional odds), g(S)=-logit(S); for "probit"
(generalised probit), g(S)=-probit(S); for "AH" (generalised additive
hazards), g(S)=-log(S); for "AO" (generalised Aranda-Ordaz), g(S)=log((S^(-theta.AO)-1)/theta.AO).}
\item{theta.AO}{theta parameter for the Aranda-Ordaz link type.}
\item{optimiser}{select which optimiser is used}
\item{log.time.transform}{should a log-transformation be used for
  calculating the derivative of the design matrix with respect to time? (default=TRUE)}
\item{recurrent}{logical for whether clustered, left truncated data are
  recurrent or for first event (where the latter requires an adjustment
  for the frailties or random effects)}
\item{frailty}{logical for whether to fit a shared frailty model}
\item{cluster}{string for the data variable that determines the cluster
  for the frailty}
\item{logtheta}{initial value for log-theta used in the gamma shared frailty
  model}
\item{nodes}{number of integration points for Gaussian quadrature}
\item{RandDist}{type of distribution for the random effect or frailty}
\item{adaptive}{logical for whether to use adaptive or non-adaptive
  quadrature}
\item{maxkappa}{double float value for the maximum value of the weight
  used in the constraint}
\item{Z}{formula for the design matrix for the random effects}
\item{reltol}{list with components for search and final relative tolerances.}
\item{outer_optim}{Integer to indicate the algorithm for outer optimisation. If outer_optim=1, then use Neldear-Mead, otherwise use Nlm.}
\item{robust_initial}{logical for whether to use Nelder-Mead
  to find initial values (max 50 iterations). This is useful for
  ill-posed initial values.}  
\item{\dots}{
  additional arguments to be passed to the \code{\link{mle2}} .
}
}
\details{
  
  The implementation extends the \code{mle2} object from the
  \code{bbmle} package.

  The default smoother for time on the linear predictor scale is s(log(time)).
}
\value{
A \code{pstpm2-class} object. 
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
\dontrun{
data(brcancer)
## standard Kaplan-Meier curves by hormon
plot(survfit(Surv(rectime/365,censrec==1)~1,data=brcancer,subset=hormon==1),
  xlab="Recurrence free survival time (years)",
  ylab="Survival")
lines(survfit(Surv(rectime/365,censrec==1)~1,data=brcancer,subset=hormon==0),col=2,
  conf.int=TRUE)
legend("topright", legend=c("Hormonal therapy","No hormonal therapy"),lty=1,col=1:2,bty="n")

## now fit a penalised stpm2 model
fit <- pstpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer)
## no S4 generic lines() method: instead, use plot(..., add=TRUE)
plot(fit,newdata=data.frame(hormon=1),type="surv",add=TRUE,ci=FALSE,line.col="blue",lwd=2,
  rug=FALSE)
plot(fit,newdata=data.frame(hormon=0),type="surv",add=TRUE,ci=FALSE,line.col="green",lwd=2,
  rug=FALSE)

## plot showing proportional hazards
plot(fit,newdata=data.frame(hormon=1),type="hazard",line.col="blue",lwd=2,
  rug=FALSE,ylim=c(0,1e-3))
plot(fit,newdata=data.frame(hormon=0),type="hazard",add=TRUE,ci=FALSE,line.col="green",lwd=2,
  rug=FALSE)

## time-varying hazard ratios
fit.tvc <- pstpm2(Surv(rectime,censrec==1)~1,
  data=brcancer,
  smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon))
plot(fit.tvc,newdata=data.frame(hormon=1),type="hazard",line.col="blue",lwd=2,
  rug=FALSE)
plot(fit.tvc,newdata=data.frame(hormon=0),type="hazard",line.col="red",lwd=2,
  add=TRUE)

## Smooth covariate effects
fit.smoothx <- pstpm2(Surv(rectime,censrec==1)~1,
  data=brcancer,
  smooth.formula=~s(log(rectime))+s(x1))
ages <- seq(21,80,length=301)
haz <- predict(fit.smoothx,newdata=data.frame(hormon=1,rectime=365,x1=ages),
               type="hazard",se.fit=TRUE)
matplot(ages,haz/haz[150,1],type="l",log="y",ylab="Hazard ratio")

## compare with df=5 from stpm2
fit.stpm2 <- stpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer,df=7)
plot(fit,newdata=data.frame(hormon=1),type="hazard",line.col="blue",lwd=2,
  rug=FALSE,ylim=c(0,1e-3))
plot(fit.stpm2,newdata=data.frame(hormon=1),type="hazard",line.col="orange",lwd=2,
  rug=FALSE,add=TRUE,ci=FALSE)

## time-varying coefficient
##summary(fit.tvc <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
##                     tvc=list(hormon=3)))
##anova(fit,fit.tvc) # compare with and without tvc (unclear whether this is valid)

## some more plots
## plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon")
                                        # no lines method: use add=TRUE
## plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon",
##     add=TRUE,ci=FALSE,line.col=2)

## plot(fit.tvc,newdata=data.frame(hormon=0),type="sdiff",var="hormon")

## plot(fit.tvc,newdata=data.frame(hormon=0),type="hdiff",var="hormon")

## plot(fit.tvc,newdata=data.frame(hormon=0),type="hazard")
## plot(fit.tvc,newdata=data.frame(hormon=1),type="hazard",line.col=2,ci=FALSE,add=TRUE)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
