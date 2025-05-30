\name{voptimize}
\title{Vectorised One Dimensional Optimization}
\alias{voptimize}
\alias{voptimise}
\usage{
voptimize(f, interval, \dots,
          lower=pmin(interval[,1], interval[,2]),
          upper=pmax(interval[,1], interval[,2]),
          maximum = FALSE,
          tol = .Machine$double.eps^0.25)
voptimise(f, interval, \dots,
          lower=pmin(interval[,1], interval[,2]),
          upper=pmax(interval[,1], interval[,2]),
          maximum = FALSE,
          tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{f}{the function to be optimized.  The function is
    either minimized or maximized over its first argument
    depending on the value of \code{maximum}.}
  \item{interval}{a matrix with two columns containing the end-points of the interval
    to be searched for the minimum.}
  \item{\dots}{additional named or unnamed arguments to be passed
    to \code{f}}
  \item{lower, upper}{the lower and upper end points of the interval to
    be searched.}
  \item{maximum}{logical.  Should we maximize or minimize (the default)?}
  \item{tol}{the desired accuracy.}
}
\description{
  The function \code{voptimize} searches the interval from
  \code{lower} to \code{upper} for a minimum or maximum of
  the vectorised function \code{f} with respect to its first argument.

  \code{optimise} is an alias for \code{optimize}.
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  The method used is a combination of golden section search and
  successive parabolic interpolation, and was designed for use with
  continuous functions.  Convergence is never much slower
  than that for a Fibonacci search.  If \code{f} has a continuous second
  derivative which is positive at the minimum (which is not at \code{lower} or
  \code{upper}), then convergence is superlinear, and usually of the
  order of about 1.324.

  The function \code{f} is never evaluated at two points closer together
  than \eqn{\epsilon}{eps *}\eqn{ |x_0| + (tol/3)}, where
  \eqn{\epsilon}{eps} is approximately \code{sqrt(\link{.Machine}$double.eps)}
  and \eqn{x_0} is the final abscissa \code{optimize()$minimum}.\cr
  If \code{f} is a unimodal function and the computed values of \code{f}
  are always unimodal when separated by at least \eqn{\epsilon}{eps *}
  \eqn{ |x| + (tol/3)}, then \eqn{x_0} approximates the abscissa of the
  global minimum of \code{f} on the interval \code{lower,upper} with an
  error less than \eqn{\epsilon}{eps *}\eqn{ |x_0|+ tol}.\cr
  If \code{f} is not unimodal, then \code{optimize()} may approximate a
  local, but perhaps non-global, minimum to the same accuracy.

  The first evaluation of \code{f} is always at
  \eqn{x_1 = a + (1-\phi)(b-a)} where \code{(a,b) = (lower, upper)} and
  \eqn{\phi = (\sqrt 5 - 1)/2 = 0.61803..}{phi = (sqrt(5) - 1)/2 = 0.61803..}
  is the golden section ratio.
  Almost always, the second evaluation is at
  \eqn{x_2 = a + \phi(b-a)}{x_2 = a + phi(b-a)}.
  Note that a local minimum inside \eqn{[x_1,x_2]} will be found as
  solution, even when \code{f} is constant in there, see the last
  example.

  \code{f} will be called as \code{f(\var{x}, ...)} for a numeric value
  of \var{x}.

  The argument passed to \code{f} has special semantics and used to be
  shared between calls.  The function should not copy it.

  The implementation is a vectorised version of the \code{optimize} function.
}

\value{
  A list with components \code{minimum} (or \code{maximum})
  and \code{objective} which give the location of the minimum (or maximum)
  and the value of the function at that point.
}
\source{
  Based on R's C translation of Fortran code \url{https://netlib.org/fmm/fmin.f}
  (author(s) unstated)
  based on the Algol 60 procedure \code{localmin} given in the reference.
}
\references{
  Brent, R. (1973)
  \emph{Algorithms for Minimization without Derivatives.}
  Englewood Cliffs, NJ: Prentice-Hall.
}
\seealso{
  \code{\link{optimize}} for the standard single optimiser solver,
  \code{\link{nlm}}, \code{\link{uniroot}}.
}
\examples{
library(graphics)

f <- function (x, a) (x - a)^2
xmin <- voptimize(f, lower=c(0, 0), upper=c(1,1), tol = 0.0001, a = c(1/3,2/3))
xmin

## See where the function is evaluated:
voptimize(function(x) x^2*(print(x)-1), lower = c(0,0), upper = c(10,10))

## "wrong" solution with unlucky interval and piecewise constant f():
f  <- function(x) ifelse(x > -1, ifelse(x < 4, exp(-1/abs(x - 1)), 10), 10)
fp <- function(x) { print(x); f(x) }

plot(f, -2,5, ylim = 0:1, col = 2)
voptimize(fp, cbind(-4, 20))   # doesn't see the minimum
voptimize(fp, cbind(-7, 20))   # ok

}
\keyword{optimize}
