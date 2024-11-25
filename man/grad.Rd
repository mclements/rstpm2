\name{grad}
\Rdversion{1.1}
\alias{grad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gradient function (internal function)
  %%  ~~function to do ... ~~
}
\description{
  Numerical gradient for a function at a given value (internal).
}
\usage{
grad(func, x, ..., method=c("fast","richardson"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
    Function taking a vector argument x (returns a vector of length>=1)
}
  \item{x}{
    vector of arguments for where the gradient is wanted.
}
  \item{\dots}{
    other arguments to the function
}
  \item{method}{
    string argument to determine whether to use the fast two sided
    calculation or use a Richardson extrapolation.
}
}
\details{
  (func(x+delta,...)-func(x-delta,...))/(2 delta) where delta is the third root
  of the machine precision times pmax(1,abs(x)).
}
\value{
  A vector if func(x) has length 1, otherwise a matrix with rows for x
  and columns for func(x).
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
  numDelta()
}
%% \examples{
%% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
