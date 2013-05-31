\name{nsx}
\alias{nsx}
\title{Generate a Basis Matrix for Natural Cubic Splines (with eXtensions)}
\description{
  Generate the B-spline basis matrix for a natural cubic spline (with eXtensions).
}
\usage{
nsx(x, df = NULL, knots = NULL, intercept = FALSE,
    Boundary.knots = range(x), derivs = if (cure) c(2, 1) else c(2, 2),
    log = FALSE, centre = FALSE, 
    cure = FALSE, stata.stpm2.compatible = FALSE)
}
\arguments{
  \item{x}{the predictor variable.  Missing values are allowed.}
  \item{df}{degrees of freedom. One can supply \code{df} rather than
    knots; \code{ns()} then chooses \code{df - 1 - intercept + 4 - sum(derivs)} knots at
    suitably chosen quantiles of \code{x} (which will ignore missing values).}
  \item{knots}{breakpoints that define the spline.  The default is no
    knots; together with the natural boundary conditions this results in
    a basis for linear regression on \code{x}.  Typical values are the
    mean or median for one knot, quantiles for more knots.  See also
    \code{Boundary.knots}.}
  \item{intercept}{if \code{TRUE}, an intercept is included in the
    basis; default is \code{FALSE}.}
  \item{Boundary.knots}{boundary points at which to impose the natural
    boundary conditions and anchor the B-spline basis (default the range
    of the data).  If both \code{knots} and \code{Boundary.knots} are
    supplied, the basis parameters do not depend on \code{x}. Data can
    extend beyond \code{Boundary.knots}}
  \item{derivs}{an integer vector of length 2 with values between 0 and
    2 giving the derivative constraint order at
    the left and right boundary knots; an order of 2 constrains the
    second derivative to zero (f''(x)=0); an order of 1 constrains the
    first and second derivatives to zero (f'(x)=f''(x)=0); an order of 1
  constrains the zero, first and second derivatives to zero (f(x)=f'(x)=f''(x)=0)}
  \item{log}{a Boolean indicating whether the underlying values have been
    log transformed; used to calculate derivatives}
  \item{centre}{if specified, then centre the splines at this value
    (i.e. f(centre)=0) (default=FALSE)}
  \item{cure}{a Boolean indicated whether to estimate cure; changes the
    default derivs argument, such that the right boundary has the first
    and second derivatives constrained to zero; defaults to FALSE}
  \item{stata.stpm2.compatible}{a Boolean to determine whether to use
    Stata stpm's default knot placement; defaults to FALSE}
}
\value{
    A matrix of dimension \code{length(x) * df} where either \code{df} was
    supplied or if \code{knots} were supplied,
    \code{df = length(knots) + 1 + intercept}.
    Attributes are returned that correspond to the arguments to \code{ns},
    and explicitly give the \code{knots}, \code{Boundary.knots} etc for
    use by \code{predict.nsx()}.

    \code{nsx()} is based on the functions \code{\link{ns}} and \code{\link{spline.des}}.  It
    generates a basis matrix for representing the family of
    piecewise-cubic splines with the specified sequence of
    interior knots, and the natural boundary conditions.  These enforce
    the constraint that the function is linear beyond the boundary knots,
    which can either be supplied, else default to the extremes of the
    data.  A primary use is in modeling formula to directly specify a
    natural spline term in a model.

    The extensions from \code{\link{ns}} are: specification of the
    derivative constraints at the boundary knots; whether to centre the
    knots; incorporation of cure using derivatives; compatible knots
    with Stata's stpm2; and an indicator for a log-transformation of
    \code{x} for calculating derivatives.
    
}
\seealso{
  \code{\link{ns}}, \code{\link{bs}}, \code{predict.nsx}, \code{\link{SafePrediction}}
}
\references{
  Hastie, T. J. (1992)
  Generalized additive models.
  Chapter 7 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\examples{
require(stats); require(graphics)
nsx(women$height, df = 5)
summary(fm1 <- lm(weight ~ ns(height, df = 5), data = women))

## example of safe prediction
plot(women, xlab = "Height (in)", ylab = "Weight (lb)")
ht <- seq(57, 73, length.out = 200)
lines(ht, predict(fm1, data.frame(height=ht)))

plotNsx <- function(obj,x,coef=c(-3,-4,2)) {
   ## numerical derivatives returning a function
   diff <- function(f,delta=1e-4,...) function(x) (f(x+delta,...)-f(x-delta,...))/(2*delta)
   par(mfrow=c(2,2))
   plot(x,predict(obj,x)%*%coef,type="l",ylab="y",main="Natural spline")
   ## first derivative
   plot(x,diff(function(x) predict(obj,x)%*%c(0.1,-1,1))(x),type="l",
      ylab="y", main="First derivative")
   ## second derivative
   plot(x,diff(diff(function(x) predict(obj,x)%*%c(0.1,-1,1)))(x),type="l",
        ylab="y", main="Second derivative")
}
x=seq(0,12,length=100)
plotNsx(nsx(1:10,df=3,derivs=c(2,2)),x)
plotNsx(nsx(1:10,df=3,derivs=c(2,1)),x)
plotNsx(nsx(1:10,df=3,derivs=c(1,1)),x)
plotNsx(nsx(1:10,df=3,derivs=c(1,0)),x) # NB: predict(obj,Boundary.knot[1])==0.0 if intercept==FALSE
rm(plotNsx,x)
}
\keyword{smooth}
