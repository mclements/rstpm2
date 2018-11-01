\name{vuniroot}
\title{Vectorised One Dimensional Root (Zero) Finding}
\alias{vuniroot}
\usage{
vuniroot(f, \dots,
        lower, upper,
        f.lower = f(lower, \dots), f.upper = f(upper, \dots),
        check.conv = FALSE,
        tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
}
\arguments{
  \item{f}{the function for which the root is sought.}
  \item{\dots}{additional named or unnamed arguments to be passed
    to \code{f}}
  \item{lower, upper}{the lower and upper end points of the interval to
    be searched.}
  \item{f.lower, f.upper}{the same as \code{f(upper)} and
    \code{f(lower)}, respectively.  Passing these values from the caller
    where they are often known is more economical as soon as \code{f()}
    contains non-trivial computations.}
  \item{check.conv}{logical indicating whether a convergence warning of the
    underlying \code{\link{vuniroot}} should be caught as an error and if
    non-convergence in \code{maxiter} iterations should be an error
    instead of a warning.}
  \item{tol}{the desired accuracy (convergence tolerance).}
  \item{maxiter}{the maximum number of iterations.}
  \item{trace}{integer number; if positive, tracing information is
    produced.  Higher values giving more details.}
}
\description{
  The function \code{vuniroot} searches the interval from \code{lower}
  to \code{upper} for a root (i.e., zero) of the vectorised function \code{f} with
  respect to its first argument.
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  Both \code{lower} and \code{upper} must be
  specified: the endpoint are re-ordered if necessary (cf. \code{uniroot}).
  The function values at the endpoints must be of opposite signs (or
  zero).  

  \code{vuniroot()} uses a C++ subroutine based on \file{"zeroin"} (from Netlib)
  and algorithms given in the reference below.  They assume a
  continuous function (which then is known to have at least one root in
  the interval).

  Convergence is declared either if \code{f(x) == 0} or the change in
  \code{x} for one step of the algorithm is less than \code{tol} (plus an
  allowance for representation error in \code{x}).

  If the algorithm does not converge in \code{maxiter} steps, a warning
  is printed and the current approximation is returned.

  \code{f} will be called as \code{f(\var{x}, ...)} for a numeric value
  of \var{x}.

  The argument passed to \code{f} has special semantics and used to be
  shared between calls.  The function should not copy it.
}

\value{
  A list with at least three components: \code{root} and \code{f.root}
  give the location of the root and the value of the function evaluated
  at that point. \code{iter} gives the number of
  iterations used.

}
\source{
  Based on \file{zeroin.c} in \url{http://www.netlib.org/c/brent.shar}.
}
\references{
  Brent, R. (1973)
  \emph{Algorithms for Minimization without Derivatives.}
  Englewood Cliffs, NJ: Prentice-Hall.
}
\seealso{
  \code{\link{uniroot}} for the standard single root solver
  \code{\link{polyroot}} for all complex roots of a polynomial;
  \code{\link{optimize}}, \code{\link{nlm}}.
}
\examples{\donttest{
require(utils) # for str

## some platforms hit zero exactly on the first step:
## if so the estimated precision is 2/3.
f <- function (x, a) x - a
str(xmin <- vuniroot(f, lower=c(0, 0), upper=c(1,1), tol = 0.0001, a = c(1/3,2/3)))

## handheld calculator example: fixed point of cos(.):
vuniroot(function(x) cos(x) - x, lower = -pi, upper = pi, tol = 1e-9)$root
}% donttest because printed output is so much platform dependent

}
\keyword{optimize}
