\name{markov_sde}
\alias{markov_sde}
\alias{standardise.markov_sde}
\alias{plot.markov_sde}
\alias{as.data.frame.markov_sde}
\title{
Predictions for continuous time, nonhomogeneous Markov multi-state
models using Aalen's additive hazards models.
}
\description{
A numerically efficient algorithm to calculate predictions from a
continuous time, nonhomogeneous Markov multi-state model. The main
inputs are are a list of Aalen's additive hazards models, the initial values,
the transition matrix and the covariate patterns. The predictions
include state occupancy probabilities and length of stay. Standard errors are calculated
using the delta method. Includes differences and standardisation.
}
\usage{
markov_sde(models, trans, newdata, init = NULL, nLebesgue = 10000 + 1, los = FALSE,
           nOut = 300, weights = 1)
\method{standardise}{markov_sde}(x, \dots)
\method{plot}{markov_sde}(x, y, stacked=TRUE, which=c("P","L"), index=NULL,
    xlab="Time", ylab=NULL, col=2:6, border=col,
    ggplot2=FALSE, lattice=FALSE, alpha=0.2,
    strata=NULL,
    \dots)
\method{as.data.frame}{markov_sde}(x, row.names=NULL, optional=NULL, ci=TRUE,
    P.conf.type="logit", L.conf.type="log",
    P.range=c(0,1), L.range=c(0,Inf),
    \dots)
}
\arguments{
  \item{models}{
list of models. Currently allows only for \code{\link[timereg]{aalen}} regression models.
}
  \item{trans}{
Transition matrix describing the states and transitions
  in the multi-state model. If S is the number of states in the
  multi-state model, \code{trans} should be an S x S matrix,
  with (i,j)-element a positive integer if a transition from i to j
  is possible in the multi-state model, \code{NA} otherwise. In particular,
  all diagonal elements should be \code{NA}. The
  integers indicating the possible transitions in the multi-state
  model should be sequentially numbered, 1,\ldots,K, with K the number
  of transitions. See \code{\link[mstate]{msprep}}
}
  \item{newdata}{
\code{\link{data.frame}} of the covariates to use in the predictions
}
  \item{init}{
vector of the initial values with the same length as the number of states. Defaults to the first state having an
initial value of 1 (i.e. \code{"[<-"(rep(0,nrow(trans)),1,1)}).
}
  \item{nLebesgue}{
Number of steps for the continuous integration
}
  \item{los}{
logical variable for whether to estimate the length of stay
}
  \item{nOut}{
number of rows to represent the continuous changes
}
  \item{weights}{
numeric vector to represent differences or standardisation
}

For \code{plot.markov_sde}:

\item{y}{(currently ignored)}
\item{stacked}{logical for whether to stack the plots. Default: TRUE}
\item{index}{indicator of which row of \code{newdata} to plot}
\item{which}{character to indicate either transition probabilities (\code{"P"}) or length of stay (\code{"L"})). Default: \code{"P"}.}
\item{xlab}{x-axis label}
\item{ylab}{x-axis label}
\item{col}{colours (ignored if \code{ggplot2=TRUE})}
\item{border}{border colours for the \code{\link{polygon}} (ignored if \code{ggplot=TRUE})}
\item{ggplot2}{use \code{ggplot2}}
\item{alpha}{alpha value for confidence bands (ggplot)}
\item{lattice}{use \code{lattice}}
\item{strata}{formula for the stratification factors for the plot}

For \code{as.data.frame.markov_sde}:

\item{row.names}{add in row names to the output data-frame}
\item{optional}{(not currently used)}
\item{ci}{logical for whether to include confidence intervals. Default:
TRUE}
\item{P.conf.type}{type of transformation for the confidence interval
calculation for the state occupancy probabilities. Default: logit transformation. This is changed to \code{"identity"} if any of the weights are negative}
\item{L.conf.type}{type of transformation for the confidence interval
calculation for the length of stay calculation. Default: log transformation. \code{"identity"} if any of the weights are negative}
\item{P.range}{valid values for the state occupancy probabilities. Default: (0,1). }
\item{L.range}{valid values for the state occupancy probabilities. Default: (0,Inf). }

For \code{standardise.markov_sde}:
\item{x}{object to extract standardised values}

  \item{\dots}{
other arguments. For \code{plot.markov_sde}, these arguments are passed to \code{\link{plot.default}}. For \code{standardise.markov_sde}, these arguments are not used, as the standardisation must be done earlier in \code{markov_sde}.
  }
}
\details{
Uses an approach developed by Ryalen and colleagues. This is a re-implementation in C++.

The current implementation only allows for a vector of initial
values rather than a matrix. The predictions will need to be re-run for
different vectors of initial values.

}
\value{
  \code{markov_sde} returns an object of \code{\link{class}} \code{"markov_sde"}.
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{
Mark Clements
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{markov_msm}}
}
%% \examples{

%% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
