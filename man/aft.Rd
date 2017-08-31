\name{aft}
\alias{aft}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric accelerated failure time model with smooth time functions
}
\description{
This implements the accelerated failure time models S_0(t exp(beta x))
and S_0(int_0^t exp(beta x(u)) du). The baseline function S_0(t*) is
modelled as exp(-exp(eta_0(log(t*)))), where eta_0(log(t*)) is a linear
predictor using natural splines.
}
\usage{
aft(formula, data, smooth.formula = NULL, df = 3,
    control = list(parscale = 1, maxit = 1000),
    init = NULL, weights = NULL, timeVar = "", time0Var = "",
    reltol = 1e-08, trace = 0, contrasts = NULL, subset = NULL,
    use.gr = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object, with the response on the left of a \code{~} operator,
and the regression terms (excluding time) on the right. The response
should be a survival object as returned by the \code{\link{Surv}}
function. The terms can include linear effects for any time-varying
coefficients. [required]
}
  \item{data}{
a data-frame in which to interpret the variables named in the
\code{formula} argument. [at present: required]
}
  \item{smooth.formula}{
a formula for describing the time effects for the linear predictor,
excluding the baseline S_0(t*), but including time-dependent
acceleration factors. The time-dependent acceleration factors can be
modelled with any smooth functions.
}
  \item{df}{
an integer that describes the degrees of freedom for the \code{ns}
function for modelling the baseline log-cumulative hazards function (default=3).
}
  \item{control}{
\code{control} argument passed to \code{optim}.
}
  \item{init}{
\code{init} should either be \code{FALSE}, such that initial
                  values will be determined using Cox regression, or a numeric
                  vector of initial values.
}
\item{weights}{
  an optional vector of 'prior weights' to be used in the
  fitting process. Should be \code{NULL} or a numeric vector.
}
  \item{timeVar}{
string variable defining the time variable. By default, this is
determined from the survival object, however this may be ambiguous if
two variables define the time.
}
\item{time0Var}{
  string variable to determine the entry variable; useful
  for when more than one data variable is used in the entry time.
}
  \item{reltol}{
relative tolerance for the model convergence
}
  \item{trace}{
integer for whether to provide trace information from the optim procedure
}
  \item{contrasts}{
an optional list. See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}.
}
  \item{subset}{
an optional vector specifying a subset of observations to
  be used in the fitting process.
}
  \item{use.gr}{
logical indicating whether to use gradients in the calculation
}
  \item{\dots}{
  additional arguments to be passed to the \code{\link{mle2}}.
}
}
\details{
  The implementation extends the \code{mle2} object from the
  \code{bbmle} package. The model inherits all of the methods from the
  \code{mle2} class.
}
\value{
An \code{stpm2-class} object that inherits from \code{mle2-class}. 
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{
Mark Clements.
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{survreg}}, \code{\link{coxph}}
}
\examples{
summary(aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ smooth }

