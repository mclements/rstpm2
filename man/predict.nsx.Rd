\name{predict.nsx}
\title{Evaluate a Spline Basis}
\alias{predict.nsx}
\description{
  Evaluate a predefined spline basis at given values.
}
\usage{
\method{predict}{nsx}(object, newx, \dots)
}
\arguments{
  \item{object}{the result of a call to \code{\link{nsx}} having attributes describing \code{knots},
    \code{degree}, etc.}
  \item{newx}{the \code{x} values at which evaluations are required.}
  \item{\dots}{Optional additional arguments.  At present no additional
    arguments are used.}
}
\value{
  An object just like \code{object}, except evaluated at the new values
  of \code{x}.

  These are methods for the generic function \code{\link{predict}} for
  objects inheriting from classes \code{"nsx"}.  See
  \code{\link{predict}} for the general behavior of this function.
}
\seealso{
  \code{\link{nsx}}.
}
\examples{
basis <- nsx(women$height, df = 5)
newX <- seq(58, 72, length.out = 51)
# evaluate the basis at the new data
predict(basis, newX)
}
\keyword{smooth}
