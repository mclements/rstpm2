\name{smoothpwc}
\alias{smoothpwc}
\title{
Utility to use a smooth function in markov_msm based on piece-wise constant values
}
\description{
Utility to use a smooth function in markov_msm based on piece-wise constant values

}
\usage{
smoothpwc(midts, rates, tmvar = "t", offsetvar = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{midts}{
mid-point values for time in each segment
}
  \item{rates}{
rates at those mid-points (or for the interval)
}
  \item{tmvar}{
string for the time variable
}
  \item{offsetvar}{
string for a time offset variable
}
  \item{\dots}{
other arguments
}
}
\details{
  Uses splines to smooth the log-rates. This assumes that the rates are
  strictly greater than zero.
}
\value{
  a function that is used in markov_msm
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (midts, rates, tmvar = "t", offsetvar = "", ...) 
{
    log.smoother <- splinefunx(midts, log(rates), constant.right = TRUE)
    haz <- function(newdata) {
        t <- newdata[[tmvar]] + (if (offsetvar != "") 
            newdata[[offsetvar]]
        else 0)
        exp(log.smoother(t))
    }
    structure(list(haz = haz), class = "smoothpwc")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
