\name{markov_msm}
\alias{markov_msm}
\alias{coef.markov_msm}
\alias{vcov.markov_msm}
\alias{as.data.frame.markov_msm}
\alias{as.data.frame.markov_msm_diff}
\alias{as.data.frame.markov_msm_ratio}
\alias{standardise.markov_msm}
\alias{plot.markov_msm}
\alias{subset.markov_msm}
\alias{diff.markov_msm}
\alias{ggplot.markov_msm}
\alias{ggplot.markov_msm_diff}
\alias{ggplot.markov_msm_ratio}
\alias{ratio_markov_msm}
\alias{cbind.markov_msm}
\alias{zeroRate}
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
utilities), length of stay and costs. 
}
\usage{
markov_msm(x, trans, t = c(0, 1), newdata = NULL, init = `[<-`(rep(0,
nrow(trans)), 1, 1), tmvar = NULL, sing.inf = 1e+10, method = "adams",
rtol = 1e-10, atol = 1e-10, slow = FALSE, min.tm = 1e-08, utility =
function(t) rep(1, nrow(trans)), use.costs = FALSE, transition.costs =
lapply(1:sum(!is.na(trans)), function(i) function(t) 1000), state.costs
= function(t) rep(0, nrow(trans)), discount.rate = 0, ...)
\method{vcov}{markov_msm}(object, ...)
\method{coef}{markov_msm}(object, ...)
\method{as.data.frame}{markov_msm}(x, row.names=NULL, optional=FALSE,
                                     ci=TRUE,
                                     P.conf.type="log-log", L.conf.type="log",
                                     P.range=c(0,1), L.range=c(0,Inf),
                                     ...)
\method{as.data.frame}{markov_msm_diff}(x, ...)
\method{as.data.frame}{markov_msm_ratio}(x, ...)
\method{standardise}{markov_msm}(x,
                                   weights = rep(1,nrow(x$newdata)),
                                   normalise = TRUE)
\method{plot}{markov_msm}(x, y, type="stacked",
                            xlab="Time", ylab="Probability", col=2:6, border=col,
                            ...)
\method{subset}{markov_msm}(x, subset, ...)
\method{diff}{markov_msm}(x, y, ...)
\method{gglot}{markov_msm}(data, mapping=NULL,
                              which=c("P","L"), 
                              stacked = TRUE, alpha=0.2,
                              ..., environment=parent.frame())
\method{gglot}{markov_msm_diff}(data, mapping=NULL, which=c("P","L"), ...)
\method{gglot}{markov_msm_ratio}(data, mapping=NULL, which=c("P","L"), ...)
ratio_markov_msm(x, y, ...)
\method{cbind}{markov_msm}(..., deparse.level=1)
\method{transform}{markov_msm}(`_data`, ...)
zeroRate(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
For \code{markov_msm}:

  \item{x}{
list of parametric or penalised survival models. Currently this includes
combinations of \code{\link{stpm2}}, \code{\link{pstpm2}}, \code{\link{glm}},
\code{\link{survPen}} or an object of class \code{"zeroRate"} from
\code{\link{zeroRate}} based on one of the other classes. The order in
the list matches the indexing in the \code{trans} argument
}
  \item{trans}{Transition matrix describing the states and transitions
  in the multi-state model. If S is the number of states in the
  multi-state model, \code{trans} should be an S x S matrix,
  with (i,j)-element a positive integer if a transition from i to j
  is possible in the multi-state model, \code{NA} otherwise. In particular,
  all diagonal elements should be \code{NA}. The
  integers indicating the possible transitions in the multi-state
  model should be sequentially numbered, 1,...,K, with K the number
  of transitions. See \code{\link{msprep}}}
  \item{t}{
numerical vector for the times to evaluation the predictions. Includes
the start time
}
  \item{newdata}{
\code{\link{data-frame}} of the covariates to use in the predictions
}
  \item{init}{
vector of the initial values with the same length as the number of states. Defaults to the first state having an
initial value of 1. 
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
method used by the ordinary differential equation solver. Defaults to
Adams method (\code{"adams"}) for non-stiff differential equations
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
  \item{use.costs}{
logical for whether to use costs. Default: FALSE
}
  \item{transition.costs}{
a function of the form \code{function(t)} that returns the cost for each transition
}
  \item{state.costs}{
a function of the form \code{function(t)} that returns the cost per unit
time for each state
}
  \item{discount.rate}{
  numerical value for the proportional reduction (per unit time) in the length of stay
and costs
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
\item{P.range}{valid values for the state occupancy probabilities. Default: (0,1). This is changed for
\code{\link{diff}} and \code{\link{ratio_markov_msm}} objects}
\item{L.range}{valid values for the state occupancy probabilities. Default: (0,Inf). This is changed for
\code{\link{diff}} and \code{\link{ratio_markov_msm}} objects}

For \code{standardise.markov_msm}:

\item{weights}{numerical vector to use in standardising the state
occupancy probabilities, length of stay and costs. Default: 1 for each observation. }
\item{normalise}{logical for whether to normalise the weights to
1. Default: TRUE}

For \code{plot.markov_msm}:

\item{x}{object of class \code{"markov_msm"}}
\item{y}{(currently ignored)}
\item{type}{type of plot. Currently only \code{"stacked"}. See
\code{ggplot.markov_msm} for a richer set of plots}
\item{xlab}{x-axis label}
\item{ylab}{x-axis label}
\item{col}{colours}
\item{border}{border colours for the \code{\link{polygon}}}

For \code{subset.markov_msm}:
\item{x}{object of class \code{"markov_msm"}}
\item{subset}{expression that is evaluated on the \code{newdata}
component of the object to filter (or restrict) for the covariates used
for predictions}

For \code{diff.markov_msm}:
\item{x}{object of class \code{"markov_msm"}}
\item{y}{object of class \code{"markov_msm"}}

For \code{zeroRate}:
\item{object}{survival regression object to be wrapped to give predict
hazards of zero}

}
\details{

The predictions are calculated using an ordinary differential equation
solver. The algorithm uses a single run of the solver to calculate the
state occupancy probabilities, length of stay, costs and their partial
derivatives with respect to the model parameters. The predictions can also be combined to
calculate differences, ratios and standardised.

The current implementation supports a list of models for each
transition. It would be useful to also allow for one model to predict
for the different transitions.

The current implementation also only allows for a vector of initial
values rather than a matrix. The predictions will need to be re-run for
different vectors of initial values.

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
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mark Clements
}
%% \note{
%% }

\seealso{
\code{\link{pmatrix.fs}}, \code{\link{probtrans}}
}
\examples{

}
\keyword{ survival }

