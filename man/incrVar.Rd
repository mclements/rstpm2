\name{incrVar}
\Rdversion{1.1}
\alias{incrVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Utility that returns a function to increment a variable in a data-frame.
}
\description{
  A functional approach to defining an increment in one or more
  variables in a data-frame. Given a variable name and an increment
  value, return a function that takes any data-frame to return a
  data-frame with incremented values.  
}
\usage{
incrVar(var, increment = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var}{
    String for the name(s) of the variable(s) to be incremented
}
  \item{increment}{
    Value that the variable should be incremented.
}
}
\details{
  Useful for defining transformations for calculating rate ratios.
}
\value{
  A function with a single \code{data} argument that increments the variables
  in the \code{data} list/data-frame. 
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (var, increment = 1) 
{
    n <- length(var)
    if (n > 1 && length(increment)==1) 
        increment <- rep(increment, n)
    function(data) {
        for (i in 1:n) {
            data[[var[i]]] <- data[[var[i]]] + increment[i]
        }
        data
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
