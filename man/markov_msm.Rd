\name{markov_msm}
\alias{markov_msm}
\alias{vcov.markov_msm}
\alias{as.data.frame.markov_msm}
\alias{as.data.frame.markov_msm_diff}
\alias{as.data.frame.markov_msm_ratio}
\alias{standardise}
\alias{standardise.markov_msm}
\alias{plot.markov_msm}
\alias{subset.markov_msm}
\alias{diff}
\alias{diff.markov_msm}
\alias{ratio_markov_msm}
\alias{rbind.markov_msm}
\alias{transform.markov_msm}
\alias{collapse_markov_msm}
\alias{zeroModel}
\alias{hrModel}
\alias{aftModel}
\alias{addModel}
\alias{hazFun}
\alias{splineFun}
\title{
Predictions for continuous time, nonhomogeneous Markov multi-state
models using parametric and penalised survival models.
}
\description{
A numerically efficient algorithm to calculate predictions from a
continuous time, nonhomogeneous Markov multi-state model. The main
inputs are the models for the transition intensities, the initial values,
the transition matrix and the covariate patterns. The predictions
include state occupancy probabilities (possibly with discounting and
utilities), length of stay and costs. Standard errors are calculated
using the delta method. Includes, differences, ratios and standardisation.
}
\usage{
markov_msm(x, trans, t = c(0,1), newdata = NULL, init=NULL,
              tmvar = NULL, 
              sing.inf = 1e+10, method="adams", rtol=1e-10, atol=1e-10, slow=FALSE,
              min.tm=1e-8,
              utility=function(t) rep(1, nrow(trans)),
              utility.sd=rep(0,nrow(trans)),
              use.costs=FALSE,
              transition.costs=function(t) rep(0, sum(!is.na(trans))), # per transition
              transition.costs.sd=rep(0,sum(!is.na(trans))),
              state.costs=function(t) rep(0,nrow(trans)), # per unit time
              state.costs.sd=rep(0,nrow(trans)),
              discount.rate = 0,
              block.size=500,
              spline.interpolation=FALSE,
              debug=FALSE,
              \dots)
\method{vcov}{markov_msm}(object, \dots)
\method{as.data.frame}{markov_msm}(x, row.names=NULL, optional=FALSE,
                                   ci=TRUE,
                                   P.conf.type="logit", L.conf.type="log",
				   C.conf.type="log",
                                   P.range=c(0,1), L.range=c(0,Inf),
				   C.range=c(0,Inf),
                                   state.weights=NULL, obs.weights=NULL,
                                   \dots)
\method{as.data.frame}{markov_msm_diff}(x, row.names=NULL, optional=FALSE,
                                   P.conf.type="plain", L.conf.type="plain",
				   C.conf.type="plain",
                                   P.range=c(-Inf,Inf), L.range=c(-Inf,Inf),
				   C.range=c(-Inf,Inf),
                                   \dots)
\method{as.data.frame}{markov_msm_ratio}(x, row.names=NULL, optional=FALSE, ...)
standardise(x, \dots)
\method{standardise}{markov_msm}(x,
                                 weights = rep(1,nrow(x$newdata)),
                                 normalise = TRUE, \dots)
\method{plot}{markov_msm}(x, y, stacked=TRUE, which=c('P','L'),
                          xlab="Time", ylab=NULL, col=2:6, border=col,
                          ggplot2=FALSE, lattice=FALSE, alpha=0.2,
                          strata=NULL,
                          \dots)
\method{subset}{markov_msm}(x, subset, \dots)
\method{diff}{markov_msm}(x, y, \dots)
ratio_markov_msm(x, y, \dots)
\method{rbind}{markov_msm}(\dots, deparse.level=1)
\method{transform}{markov_msm}(`_data`, \dots)
collapse_markov_msm(object, which=NULL, sep="; ")
zeroModel(object)
hrModel(object,hr=1,ci=NULL,seloghr=NULL)
aftModel(object,af=1,ci=NULL,selogaf=NULL)
addModel(\dots)
hazFun(f, tmvar="t", \dots)
splineFun(time,rate,method="natural",scale=1,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
For \code{markov_msm}:

  \item{x}{
list of functions or parametric or penalised survival models. Currently
the models include
combinations of \code{\link{stpm2}}, \code{\link{pstpm2}}, \code{\link[stats]{glm}},
\code{\link[mgcv]{gam}},
\code{\link[survPen]{survPen}} or an object of class \code{"zeroModel"} from
\code{\link{zeroModel}} based on one of the other classes. The order in
the list matches the indexing in the \code{trans} argument. The
functions can optionally use a \code{t} argument for time and/or a
\code{newdata} argument. Uncertainty in the models are incorporated into
the gradients, while uncertainty in the functions are currently not modelled.
}
  \item{trans}{Transition matrix describing the states and transitions
  in the multi-state model. If S is the number of states in the
  multi-state model, \code{trans} should be an S x S matrix,
  with (i,j)-element a positive integer if a transition from i to j
  is possible in the multi-state model, \code{NA} otherwise. In particular,
  all diagonal elements should be \code{NA}. The
  integers indicating the possible transitions in the multi-state
  model should be sequentially numbered, 1,\ldots,K, with K the number
  of transitions. See \code{\link[mstate]{msprep}}}
  \item{t}{
numerical vector for the times to evaluation the predictions. Includes
the start time
}
  \item{newdata}{
\code{\link{data.frame}} of the covariates to use in the predictions
}
  \item{init}{
vector of the initial values with the same length as the number of states. Defaults to the first state having an
initial value of 1 (i.e. \code{"[<-"(rep(0,nrow(trans)),1,1)}).
}
  \item{tmvar}{
specifies the name of the time variable. This should be set for
regression models that do not specify this (e.g. \code{\link{glm}}) or
where the time variable is ambiguous
}
  \item{sing.inf}{If there is a singularity in the observed hazard,
    for example a Weibull distribution with \code{shape < 1} has infinite
    hazard at \code{t=0}, then as a workaround, the hazard is assumed to
    be a large finite number, \code{sing.inf}, at this time.   The
    results should not be sensitive to the exact value assumed, but
    users should make sure by adjusting this parameter in these cases.
  }
  \item{method}{
For \code{markov_msm}, the method used by the ordinary differential equation solver. Defaults to
Adams method (\code{"adams"}) for non-stiff differential equations.

For \code{splineFun}, the method jused for spline interpolation; see \code{splinefun}.
}
  \item{rtol }{relative error tolerance, either a
    scalar or an array as long as the number of states. Passed to \code{\link{lsode}}
  }
  \item{atol }{absolute error tolerance, either a scalar or an array as
    long as the number of states. Passed to \code{\link{lsode}}
  }
  \item{slow}{
logical to show whether to use the slow \code{R}-only
implementation. Useful for debugging. Currently needed for costs.
}
  \item{min.tm}{
  Minimum time used for evaluations. Avoids log(0) for some models.
}
  \item{utility}{
  a function of the form \code{function(t)} that returns a utility for
each state at time \code{t} for the length of stay values
}
  \item{utility.sd}{
  a function of the form \code{function(t)} that returns the standard
  deviation for the utility for
each state at time \code{t} for the length of stay values
}
  \item{use.costs}{
logical for whether to use costs. Default: FALSE
}
  \item{transition.costs}{
a function of the form \code{function(t)} that returns the cost for each transition
}
  \item{transition.costs.sd}{
a function of the form \code{function(t)} that returns the standard deviation
for the cost for each transition
}
  \item{state.costs}{
a function of the form \code{function(t)} that returns the cost per unit
time for each state
}
  \item{state.costs.sd}{
a function of the form \code{function(t)} that returns the standard deviation
 for the cost per unit time for each state
}
  \item{discount.rate}{
  numerical value for the proportional reduction (per unit time) in the length of stay
and costs
}
  \item{block.size}{
  divide \code{newdata} into blocks. Uses less memory but is slower. Reduce this number if the function call runs out of memory. 
}
  \item{spline.interpolation}{
  logical for whether to use spline interpolation for the transition
  hazards rather than the model predictions directly (default=TRUE).
}
  \item{debug}{
  logical flag for whether to keep the full output from the ordinary differential equation in the \code{res} component (default=\code{FALSE}). 
}
  \item{\dots}{
other arguments. For \code{markov_msm}, these are passed to the \code{\link{ode}} solver from the
\code{deSolve} package. For \code{plot.markov_msm}, these arguments are passed to \code{\link{plot.default}}
}

For \code{as.data.frame.markov_msm}:

\item{row.names}{add in row names to the output data-frame}
\item{optional}{(not currently used)}
\item{ci}{logical for whether to include confidence intervals. Default:
TRUE}
\item{P.conf.type}{type of transformation for the confidence interval
calculation for the state occupancy probabilities. Default: log-log transformation. This is changed for
\code{\link{diff}} and \code{\link{ratio_markov_msm}} objects}
\item{L.conf.type}{type of transformation for the confidence interval
calculation for the length of stay calculation. Default: log transformation. This is changed for
\code{\link{diff}} and \code{\link{ratio_markov_msm}} objects}
\item{C.conf.type}{type of transformation for the confidence interval
calculation for the length of stay calculation. Default: log transformation. This is changed for
\code{\link{diff}} and \code{\link{ratio_markov_msm}} objects}
\item{P.range}{valid values for the state occupancy probabilities. Default: (0,1). This is changed for
\code{\link{diff}} and \code{\link{ratio_markov_msm}} objects}
\item{L.range}{valid values for the state occupancy probabilities. Default: (0,Inf). This is changed for
\code{\link{diff}} and \code{\link{ratio_markov_msm}} objects}
\item{C.range}{valid values for the state occupancy probabilities. Default: (0,Inf). This is changed for
\code{\link{diff}} and \code{\link{ratio_markov_msm}} objects}
\item{state.weights}{Not currently documented}
\item{obs.weights}{Not currently documented}

For \code{standardise.markov_msm}:

\item{weights}{numerical vector to use in standardising the state
occupancy probabilities, length of stay and costs. Default: 1 for each observation. }
\item{normalise}{logical for whether to normalise the weights to
1. Default: TRUE}

For \code{plot.markov_msm}:

\item{y}{(currently ignored)}
\item{stacked}{logical for whether to stack the plots. Default: TRUE}
\item{xlab}{x-axis label}
\item{ylab}{x-axis label}
\item{col}{colours (ignored if \code{ggplot2=TRUE})}
\item{border}{border colours for the \code{\link{polygon}} (ignored if \code{ggplot=TRUE})}
\item{ggplot2}{use \code{ggplot2}}
\item{alpha}{alpha value for confidence bands (ggplot)}
\item{lattice}{use \code{lattice}}
\item{strata}{formula for the stratification factors for the plot}

For \code{subset.markov_msm}:
\item{subset}{expression that is evaluated on the \code{newdata}
component of the object to filter (or restrict) for the covariates used
for predictions}

For \code{transform.markov_msm}:
\item{_data}{an object of class \code{"markov_msm"}}

For \code{rbind.markov_msm}:
\item{deparse.level}{not currently used}

For \code{collapse.states}:
\item{which}{either an index of the states to collapse or a character vector of the state names to collapse}
\item{sep}{separator to use for the collapsed state names}

For \code{zeroModel} to predict zero rates:
\item{object}{survival regression object to be wrapped}

For \code{hrModel} to predict rates times a hazard ratio:
\item{hr}{hazard ratio}
\item{seloghr}{alternative specification for the se of the log(hazard ratio); see also \code{ci} argument}

For \code{aftModel} to predict accelerated rates:
\item{af}{acceleration factor}
\item{selogaf}{alternative specification for the se of the log(acceleration factor); see also \code{ci} argument}

\code{addModel} predict rates based on adding rates from different models

\code{hazFun} provides a rate function without uncertainty:
\item{f}{rate function, possibly with \code{tmvar} and/or \code{newdata} as arguments}

\code{splineFun} predicts rates using spline interpolation:
\item{time}{exact times}
\item{rate}{rates as per \code{time}}
\item{scale}{rate multiplier (e.g. \code{scale=365.25} for converting from daily rates to yearly rates)}

}
\details{

The predictions are calculated using an ordinary differential equation
solver. The algorithm uses a single run of the solver to calculate the
state occupancy probabilities, length of stay, costs and their partial
derivatives with respect to the model parameters. The predictions can also be combined to
calculate differences, ratios and standardised.

The current implementation supports a list of models for each
transition. 

The current implementation also only allows for a vector of initial
values rather than a matrix. The predictions will need to be re-run for
different vectors of initial values.

For \code{as.data.frame.markov_msm_ratio}, the data are provided in
log form, hence the default transformations and bounds are as per
\code{as.data.frame.markov_msm_diff}, with untransformed data on the
real line.

TODO: allow for one model to predict for the different transitions.

}
\value{
  \code{markov_msm} returns an object of \code{\link{class}} \code{"markov_msm"}.

  The function \code{summary} is used to
  obtain and print a summary and analysis of variance table of the
  results.  The generic accessor functions \code{coef} and \code{vcov} extract
  various useful features of the value returned by \code{markov_msm}.

  An object of class \code{"markov_msm"} is a list containing at least the
  following components:


 \item{time}{a numeric vector with the times for the predictions}
 \item{P}{an \code{\link{array}} for the predicted state occupancy
 probabilities. The array has three dimensions: time, state, and observations.}
 \item{L}{an \code{\link{array}} for the predicted sojourn times (or
 length of stay). The array has three dimensions: time, state, and observations.}
 \item{Pu}{an \code{\link{array}} for the partial derivatives of the
 predicted state occupancy probabilities with respect to the model coefficients. The array has
 four dimensions: time, state, coefficients, and observations.}
 \item{Lu}{an \code{\link{array}} for the partial derivatives of the predicted sojourn times (or
 length of stay) with respect to the model coefficients. The array has
 four dimensions: time, state, coefficients, and observations.}
 \item{newdata}{a \code{\link{data.frame}} with the covariates used for
 the predictions}
 \item{vcov}{the variance-covariance matrix for the models of the
 transition intensities}
 \item{trans}{copy of the \code{trans} input argument}
 \item{call}{the call to the function}

 For debugging:
 \item{res}{data returned from the ordinary differential equation
 solver. This may include more information on the predictions}
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{
Mark Clements
}
%% \note{
%% }

\seealso{
\code{\link[flexsurv]{pmatrix.fs}}, \code{\link[mstate]{probtrans}}
}
\examples{
\dontrun{
library(readstata13)
library(mstate)
library(ggplot2)
library(survival)
## Two states: Initial -> Final
## Note: this shows how to use markov_msm to estimate survival and risk probabilities based on
## smooth hazard models.
two_states <- function(model, ...) {
    transmat = matrix(c(NA,1,NA,NA),2,2,byrow=TRUE)
    rownames(transmat) <- colnames(transmat) <- c("Initial","Final")
    rstpm2::markov_msm(list(model), ..., trans = transmat)
}
## Note: the first argument is the hazard model. The other arguments are arguments to the
## markov_msm function, except for the transition matrix, which is defined by the new function.
death = gsm(Surv(time,status)~factor(rx), data=survival::colon, subset=(etype==2), df=3)
cr = two_states(death, newdata=data.frame(rx="Obs"), t = seq(0,2500, length=301))
plot(cr,ggplot=TRUE)

## Competing risks
## Note: this shows how to adapt the markov_msm model for competing risks.
competing_risks <- function(listOfModels, ...) {
    nRisks = length(listOfModels)
    transmat = matrix(NA,nRisks+1,nRisks+1)
    transmat[1,1+(1:nRisks)] = 1:nRisks
    rownames(transmat) <- colnames(transmat) <- c("Initial",names(listOfModels))
    rstpm2::markov_msm(listOfModels, ..., trans = transmat)
}
## Note: The first argument for competing_risks is a list of models. Names from that list are
## used for labelling the states. The other arguments are as per the markov_msm function,
## except for the transition matrix, which is defined by the competing_risks function.
recurrence = gsm(Surv(time,status)~factor(rx), data=survival::colon, subset=(etype==1), df=3)
death = gsm(Surv(time,status)~factor(rx), data=survival::colon, subset=(etype==2), df=3)
cr = competing_risks(list(Recurrence=recurrence,Death=death),
                     newdata=data.frame(rx=levels(survival::colon$rx)),
                     t = seq(0,2500, length=301))
## Plot the probabilities for each state for three different treatment arms
plot(cr, ggplot=TRUE) + facet_grid(~ rx)
## And: differences in probabilities
cr_diff = diff(subset(cr,rx=="Lev+5FU"),subset(cr,rx=="Obs"))
plot(cr_diff, ggplot=TRUE, stacked=FALSE)

## Extended example: Crowther and Lambert (2017)
## library(rstpm2); library(readstata13); library(ggplot2)
mex.1 <- read.dta13("http://fmwww.bc.edu/repec/bocode/m/multistate_example.dta")
transmat <- rbind("Post-surgery"=c(NA,1,2), 
                  "Relapsed"=c(NA,NA,3),
                  "Died"=c(NA,NA,NA))
colnames(transmat) <- rownames(transmat)
mex.2 <- transform(mex.1,osi=(osi=="deceased")+0)
levels(mex.2$size)[2] <- ">20-50 mm" # fix typo
mex <- mstate::msprep(time=c(NA,"rf","os"),status=c(NA,"rfi","osi"),
                      data=mex.2,trans=transmat,id="pid",
                      keep=c("age","size","nodes","pr_1","hormon"))
mex <- transform(mex,
                 size2=(unclass(size)==2)+0, # avoids issues with TRUE/FALSE
                 size3=(unclass(size)==3)+0,
                 hormon=(hormon=="yes")+0,
                 Tstart=Tstart/12,
                 Tstop=Tstop/12)
##
c.ar <- stpm2(Surv(Tstart,Tstop,status) ~ age + size2 + size3 + nodes + pr_1 + hormon,
              data = mex, subset=trans==1, df=3, tvc=list(size2=1,size3=1,pr_1=1))
c.ad <- stpm2(Surv(Tstart, Tstop, status) ~ age + size + nodes + pr_1 + hormon,
              data = mex, subset=trans==2, df=1)
c.rd <- stpm2( Surv(Tstart,Tstop,status) ~ age + size + nodes + pr_1 + hormon,
              data=mex, subset=trans==3, df=3, tvc=list(pr_1=1))
##
nd <- expand.grid(nodes=seq(0,20,10), size=levels(mex$size))
nd <- transform(nd, age=54, pr_1=3, hormon=0,
                size2=(unclass(size)==2)+0,
                size3=(unclass(size)==3)+0)
## Predictions
system.time(pred1 <- rstpm2::markov_msm(list(c.ar,c.ad,c.rd), t = seq(0,15,length=301),
                                        newdata=nd, trans = transmat)) # ~2 seconds
pred1 <- transform(pred1, Nodes=paste("Nodes =",nodes), Size=paste("Size",size))
## Figure 3
plot(pred1, ggplot=TRUE) + facet_grid(Nodes ~ Size) + xlab("Years since surgery")
plot(pred1, ggplot=TRUE, flipped=TRUE) +
    facet_grid(Nodes ~ Size) + xlab("Years since surgery")
plot(pred1, strata=~nodes+size, xlab="Years since surgery", lattice=TRUE)
## Figure 4
plot(subset(pred1, nodes==0 & size=="<=20 mm"), stacked=FALSE, ggplot=TRUE) +
    facet_grid(. ~ state) +
    xlab("Years since surgery")
## Figure 5
a <- diff(subset(pred1,nodes==0 & size=="<=20 mm"),
          subset(pred1,nodes==0 & size==">20-50 mm"))
a <- transform(a, label = "Prob(Size<=20 mm)-Prob(20mm<Size<50mm)")
b <- ratio_markov_msm(subset(pred1,nodes==0 & size=="<=20 mm"),
                      subset(pred1,nodes==0 & size==">20-50 mm"))
b <- transform(b,label="Prob(Size<=20 mm)-Prob(20mm<Size<50mm)")
##
c <- diff(subset(pred1,nodes==0 & size=="<=20 mm"),
          subset(pred1,nodes==0 & size==">50 mm"))
c <- transform(c, label = "Prob(Size<=20 mm)-Prob(Size>=50mm)")
d <- ratio_markov_msm(subset(pred1,nodes==0 & size=="<=20 mm"),
                      subset(pred1,nodes==0 & size==">50 mm"))
d <- transform(d,label= "Prob(Size<=20 mm)-Prob(Size>=50mm)")
##
e <- diff(subset(pred1,nodes==0 & size==">20-50 mm"),
          subset(pred1,nodes==0 & size==">50 mm"))
e <- transform(e,label="Prob(20mm<Size<50 mm)-Prob(Size>=50mm)")
f <- ratio_markov_msm(subset(pred1,nodes==0 & size==">20-50 mm"),
                      subset(pred1,nodes==0 & size==">50 mm"))
f <- transform(f, label = "Prob(20mm<Size<50 mm)-Prob(Size>=50mm)")
## combine
diffs <- rbind(a,c,e)
ratios <- rbind(b,d,f)
## Figure 5
plot(diffs, stacked=FALSE, ggplot2=TRUE) + xlab("Years since surgery") +
    ylim(c(-0.4, 0.4)) + facet_grid(label ~ state)
##
plot(ratios, stacked=FALSE, ggplot2=TRUE) + xlab("Years since surgery") +
    ylim(c(0, 3)) + facet_grid(label ~ state)
## Figure 6
plot(subset(pred1, nodes==0 & size=="<=20 mm"), stacked=FALSE, which="L", ggplot2=TRUE) +
    facet_grid(. ~ state) + xlab("Years since surgery")
## Figure 7
plot(diffs, stacked=FALSE, which="L", ggplot2=TRUE) + xlab("Years since surgery") +
    ylim(c(-4, 4)) + facet_grid(label ~ state)
plot(ratios, stacked=FALSE, which="L", ggplot2=TRUE) + xlab("Years since surgery") +
    ylim(c(0.1, 10)) + coord_trans(y="log10") + facet_grid(label ~ state)
}
}
\keyword{ survival }

