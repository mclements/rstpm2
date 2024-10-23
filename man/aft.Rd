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
    tvc = NULL, cure.formula = ~1, control = list(),
    init = NULL, weights = NULL, tvc.intercept = TRUE,
    tvc.integrated = FALSE,
    timeVar = "", time0Var = "",
    cure = FALSE, mixture = FALSE, contrasts = NULL, subset = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object, with the response on the left of a \code{~} operator,
and the regression terms (excluding time) on the right. The response
should be a survival object as returned by the \code{\link[survival]{Surv}}
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
  \item{tvc}{
a list with the names of the time-varying coefficients. This uses natural splines
(e.g. \code{tvc=list(hormon=3)} is equivalent to
\code{smooth.formula=~...+hormon:nsx(log(time),df=3)}), which by default
does \emph{not} include an intercept (or main effect) term. 
}
  \item{cure.formula}{
a formula for describing the cure fraction.
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
  \item{tvc.intercept}{
logical for whether to include an intercept in the time-varying
acceleration factor (defaults to TRUE)
}
  \item{tvc.integrated}{
logical for whether the time-varying
acceleration factor should be based on a integration, rather than a
cumulative effect (defaults to FALSE)
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
  \item{cure}{
logical for whether to model for cure using a non-mixture model (default=FALSE)
}
  \item{mixture}{
logical for whether to model for cure using a mixture model (default=FALSE)
}
  \item{contrasts}{
an optional list. See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}.
}
  \item{subset}{
an optional vector specifying a subset of observations to
  be used in the fitting process.
}
  \item{\dots}{
  additional arguments to be passed to the \code{\link[bbmle]{mle2}}.
}
}
\details{
  The implementation extends the \code{mle2} object from the
  \code{bbmle} package. The model inherits all of the methods from the
  \code{mle2} class.
}
\value{
An \code{aft-class} object that inherits from \code{mle2-class}. 
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
\code{\link[survival]{survreg}}, \code{\link[survival]{coxph}}
}
\examples{
summary(aft(Surv(rectime,censrec==1)~hormon,data=brcancer,df=4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ smooth }

