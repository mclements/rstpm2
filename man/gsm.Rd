\name{gsm}
\alias{gsm}
\alias{pstpm2}
\alias{stpm2}
\title{
  Parametric and penalised generalised survival models
}
\description{
  This implements the generalised survival model g(S(t|x)) = eta, where
  g is a link function, S is survival, t is time, x are covariates and eta
  is a linear predictor. The linear predictor can include either
  parametric or penalised
  smoothers for the time effects, for time:covariate interactions and for
  covariate effects. The main model assumption is that the time
  effects in the linear predictor are smooth. This extends the class of
  flexible parametric survival models developed by Royston and
  colleagues. The model has been extended to include relative survival
  (excess hazards),
  Gamma frailties and normal random effects.
}
\usage{
gsm(formula, data, smooth.formula = NULL, smooth.args = NULL,
                df = 3, cure = FALSE,
                tvc = NULL, tvc.formula = NULL,
                control = list(), init = NULL,
                weights = NULL, robust = FALSE, baseoff = FALSE,
                timeVar = "", time0Var = "", use.gr = NULL,
                optimiser=NULL, log.time.transform=TRUE,
                reltol=NULL, trace = NULL,
                link.type=c("PH","PO","probit","AH","AO"), theta.AO=0,
                contrasts = NULL, subset = NULL,
                robust_initial=NULL,
                coxph.strata = NULL, coxph.formula = NULL,
                logH.formula = NULL, logH.args = NULL,
                bhazard = NULL, bhazinit=NULL, copula=FALSE,
                frailty = !is.null(cluster) & !robust & !copula,
                cluster = NULL, logtheta=NULL,
                nodes=NULL, RandDist=c("Gamma","LogN"), recurrent = FALSE,
                adaptive = NULL, maxkappa = NULL,
                sp=NULL, criterion=NULL, penalty=NULL,
                smoother.parameters=NULL, Z=~1, outer_optim=NULL,
                alpha=1, sp.init=1,
                penalised=FALSE,
                \dots)
stpm2(formula, data, weights=NULL, subset=NULL, coxph.strata=NULL, \dots)
pstpm2(formula, data, weights=NULL, subset=NULL, coxph.strata=NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object, with the response on the left of a \code{~} operator, and
the parametric terms on the right.  The response must be a survival object as
returned by the \code{\link{Surv}} function. Specials include
\code{cluster} and \code{bhazard}. [required]
}
  \item{data}{
a data.frame in which to interpret the variables named in
the \code{formula} argument.
}
  \item{smooth.formula}{
    either a parametric formula or a penalised \code{mgcv::gam} formula for
    describing the time effects and
    time-dependent effects and smoothed covariate effects on the linear
    predictor scale
    (default=NULL). The default model is equal to \code{~s(log(time),k=-1)}
    where \code{time} is the time variable.
}
  \item{df}{
an integer that describes the degrees of freedom for the \code{ns}
function for modelling the baseline log-cumulative hazard
(default=3). Parametric model only.
}
  \item{smooth.args}{
a list describing the arguments for the \code{s} function for modelling
the baseline time effect on the linear predictor scale (default=NULL). 
}
  \item{tvc}{
a list with the names of the time-varying coefficients. For a parametric
model, this uses natural splines (e.g. \code{tvc=list(hormon=3)} is equivalent
to
\code{smooth.formula=~...+hormon:nsx(log(time),df=3)}), which by default
does \emph{not} include an intercept (or main effect) term. For a penalised
model, this uses cubic splines
(e.g. \code{tvc=list(hormon=-1)} is equivalent to
\code{smooth.formula=~...+s(log(time),by=hormon,k=-1))}, which by
default \emph{does} include an intercept (or main effect) term (and this
code will remove any main effect from \code{formula}). 
}
  \item{tvc.formula}{
separate formula for the time-varying effects. This is combined with
\code{smooth.formula} or the default \code{smooth.formula}.
}
\item{baseoff}{Boolean used to determine whether fully define the model
using \code{tvc.formula} rather than combining \code{logH.formula} and \code{tvc.formula}}
  \item{logH.args}{
as per \code{smooth.args}. Deprecated.
}
  \item{logH.formula}{
    as per \code{smooth.formula}. Deprecated.
  }
\item{cure}{logical for whether to estimate a cure model (parametric
  model only).}
\item{control}{list of arguments passed to \code{\link{gsm.control}}.}
\item{init}{\code{init} should either be \code{NULL}, such that initial
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
initial values. Default=0.1. Deprecated argument: use of the
\code{control} argument is preferred.
}
\item{copula}{logical to indicate whether to use a copula model (experimental)}
\item{timeVar}{variable defining the time variable. By default, this is
determined from the survival object, however this may be ambiguous if
two variables define the time}
% \item{time0Var}{variable defining the entry time variable. By default, this is
% determined from the survival object, however this may be ambiguous if
% two variables define the entry time}
\item{sp}{fix the value of the smoothing parameters.}
\item{use.gr}{in R, a Boolean to determine whether to use the gradient
  in the optimisation. Default=TRUE, Deprecated argument: use of the
\code{control} argument is preferred.}
\item{criterion}{in Rcpp, determine whether to use "GCV" or "BIC" for for the smoothing parameter selection.}
\item{penalty}{use either the "logH" penalty, which is the default
  penalty from mgcv, or the "h" hazard penalty. Default="logH". Deprecated argument: use of the
\code{control} argument is preferred.}
\item{smoother.parameters}{for the hazard penalty, a list with components which are lists with components var, transform and inverse.}
\item{alpha}{an ad hoc tuning parameter for the smoothing parameter.}
\item{sp.init}{initial values for the smoothing parameters.}
\item{trace}{integer for trace reporting; 0 represents no additional
  reporting. Default=0. Deprecated argument: use of the
\code{control} argument is preferred.}
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
\item{optimiser}{select which optimiser is used. Default="BFGS". Deprecated argument: use of the
\code{control} argument is preferred.}
\item{log.time.transform}{should a log-transformation be used for
  calculating the derivative of the design matrix with respect to time? (default=TRUE)}
\item{recurrent}{logical for whether clustered, left truncated data are
  recurrent or for first event (where the latter requires an adjustment
  for the frailties or random effects)}
\item{frailty}{logical for whether to fit a shared frailty model}
\item{cluster}{variable that determines the cluster
  for the frailty. This can be a vector, a string for the column, or a
  name. This can also be specified using a special.}
\item{logtheta}{initial value for log-theta used in the gamma shared frailty
  model (defaults to value from a \code{coxph} model fit)}
\item{nodes}{number of integration points for Gaussian
  quadrature. Default=9. Deprecated argument: use of the
\code{control} argument is preferred.}
\item{RandDist}{type of distribution for the random effect or frailty}
\item{adaptive}{logical for whether to use adaptive or non-adaptive
  quadrature, Default=TRUE. Deprecated argument: use of the
\code{control} argument is preferred.}
\item{maxkappa}{double float value for the maximum value of the weight
  used in the constraint. Default=1000. Deprecated argument: use of the
\code{control} argument is preferred.}
\item{Z}{formula for the design matrix for the random effects}
\item{reltol}{list with components for search and final relative
  tolerances. Default=list(search=1e-10, final=1e-10, outer=1e-5). Deprecated argument: use of the
\code{control} argument with arguments \code{reltol.search},
\code{reltol.final} and \code{reltol.outer} is preferred.}
\item{outer_optim}{Integer to indicate the algorithm for outer
  optimisation. If outer_optim=1 (default), then use Neldear-Mead,
  otherwise use Nlm.} 
\item{robust_initial}{logical for whether to use Nelder-Mead
  to find initial values (max 50 iterations). This is useful for
  ill-posed initial values. Default= FALSE. Deprecated argument: use of the
  \code{control} argument is preferred.}
\item{penalised}{logical to show whether to use penalised models with
  \code{pstpm} (\code{penalised=TRUE}) or parametrics models with
  \code{stpm2} (\code{penalised=FALSE}).}
\item{\dots}{ additional arguments to be passed to the \code{\link{mle2}}.}
}
\details{
  
  The implementation extends the \code{mle2} object from the
  \code{bbmle} package.

  The default smoothers for time on the linear predictor scale are
  \code{nsxs(log(time),df=3)} for the parametric model and
  \code{s(log(time))} for the penalised model.

  A frequently asked question is: why does rstpm2 give different spline
  estimates to flexsurv and Stata's stpm2? The short answer is that
  rstpm2 uses a different natural spline basis compared with flexsurv
  and Stata's stpm2 and slightly different
  knot placement than Stata's stpm2. If the knot
  placement is the same, then the predictions and other coefficients are
  expected to be very similar. As a longer answer, the default smoother
  in rstpm2 is to use an extension of the \code{splines::ns} function
  (\code{rstpm2::nsx}), which uses a QR projection of B-splines for
  natural splines. In contrast, flexsurv and Stata's stpm2 use truncated
  power splines for the natural spline basis (also termed 'restricted
  cubic splines'). The B-splines are known to
  have good numerical properties, while Stata's stpm2 implementation
  defaults to using matrix orthogonalisation to account for any
  numerical instability in the truncated power basis. Furthermore,
  rstpm2 allows for any smooth parametric function to be used as a
  smoother in \code{stpm2}/\code{gsm}, which is an extension over
  flexsurv and Stata's stpm2. Finally, it may be difficult to get rstpm2 and
  Stata's stpm2 to return the same estimates: although
  \code{nsx} includes an argument \code{stata.stpm2.compatible = FALSE} (change to \code{TRUE} for
  compatibility), the design matrix for rstpm2 is based on individuals
  with events, while Stata's stpm2 determines the spline knots from the
  individuals with events and the design matrix is otherwise based on all individuals.
}
\value{
Either a \code{stpm2-class} or \code{pstpm2-class} object. 
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
  Mark Clements, Xing-Rong Liu, Benjamin Christoffersen.
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
    summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3))
    
    ## some predictions
    head(predict(fit,se.fit=TRUE,type="surv"))
    head(predict(fit,se.fit=TRUE,type="hazard"))
    
    ## some plots
    plot(fit,newdata=data.frame(hormon=0),type="hazard")
    plot(fit,newdata=data.frame(hormon=0),type="surv")

    ## time-varying coefficient
    summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                             tvc=list(hormon=3)))
    anova(fit,fit.tvc) # compare with and without tvc
    
    ## some more plots
    plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon", ylim=c(0,2))
    lines(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon",
          col=2)

    plot(fit.tvc,newdata=data.frame(hormon=0),type="sdiff",var="hormon")
    
    plot(fit.tvc,newdata=data.frame(hormon=0),type="hdiff",var="hormon")

    library(scales)
    cols <- c(alpha("red",alpha=0.2), alpha("blue",alpha=0.2))
    plot(fit.tvc,newdata=data.frame(hormon=0),type="hazard",ci.col=cols[1])
    lines(fit.tvc,newdata=data.frame(hormon=1),type="hazard",lty=2,ci.col=cols[2],
          ci=TRUE)
    legend("topright",legend=c("No hormonal treatment", "(95% CI)", "Hormonal treatment", "(95% CI)"),
           lty=c(1,1,2,1), lwd=c(1,10,1,10), col=c("black",cols[1],"black",cols[2]), bty="n")

    
    ## compare number of knots
    hormon0 <- data.frame(hormon=0)
    plot(fit,type="hazard",newdata=hormon0)
    AIC(fit)
    for (df in 4:6) {
        fit.new <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=df)
        plot(fit.new,type="hazard",newdata=hormon0,add=TRUE,ci=FALSE,line.col=df)
        print(AIC(fit.new))
    }

    ## compatibility with Stata's stpm2 using the smooth.formula argument (see Details)
    summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                  smooth.formula=~nsx(log(rectime),df=3,stata.stpm2.compatible=TRUE)))
    summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                  smooth.formula=~nsx(log(rectime),df=3,stata=TRUE)+
                  hormon:nsx(log(rectime),df=3,stata=TRUE)))

    }
}
\keyword{ survival }
