\name{gsm.control}
\alias{gsm.control}
\title{
Defaults for the gsm call
}
\description{
Set useful default and allow changes for the gsm call. This is meant to
make the gsm call simpler.
}
\usage{
gsm.control(parscale = 1, maxit = 300, optimiser = c("BFGS", "NelderMead"), trace = 0,
            nodes = 9, adaptive = TRUE, kappa.init = 1, maxkappa = 1000,
            suppressWarnings.coxph.frailty = TRUE, robust_initial = FALSE, bhazinit = 0.1,
            eps.init = 1e-5, use.gr = TRUE, penalty = c("logH", "h"), outer_optim = 1,
            reltol.search = 1e-10, reltol.final = 1e-10, reltol.outer = 1e-05,
            criterion = c("GCV", "BIC"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parscale}{
numeric vector or scalar for the scaling of the parameter values;
default 1
}
  \item{maxit}{
integer for the maximum number of iterations for the optimisation process
}
  \item{optimiser}{
which optimiser to use for the outer optimisation
}
  \item{trace}{
integer indicating the trace level for each optimiser
}
  \item{nodes}{
number of quadrature nodes
}
  \item{adaptive}{
logical for whether to use adaptive or non-adaptive
  quadrature, Default=TRUE. 
}
  \item{kappa.init}{
initial value for the quadratic penalty for inequality constraints
}
  \item{eps.init}{
initial value for epsilon
}
  \item{maxkappa}{
double float value for the maximum value of the weight
  used in the constraint.
}
  \item{suppressWarnings.coxph.frailty}{
logical
}
  \item{robust_initial}{
Not currently documented.
}
  \item{bhazinit}{
Not currently documented.
}
  \item{use.gr}{
Logical for whether to use gradients.
}
  \item{penalty}{
Not currently documented.
}
  \item{outer_optim}{
Not currently documented.
}
  \item{reltol.search}{
Relative tolerance. Not currently documented.
}
  \item{reltol.final}{
Relative tolerance. Not currently documented.
}
  \item{reltol.outer}{
Relative tolerance. Not currently documented.
}
  \item{criterion}{
Not currently documented.
}
}
%% \details{
%% }
%% \value{
%% %%  ~Describe the value returned
%% %%  If it is a LIST, use
%% %%  \item{comp1 }{Description of 'comp1'}
%% %%  \item{comp2 }{Description of 'comp2'}
%% %% ...
%% }
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
%% \examples{
%% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
