\name{coef<-}
\alias{coef<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic method to update the coef in an object.
}
\description{
Generic method to update the coef in an object.
}
\usage{
coef(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object to be updated
}
  \item{value}{
value of the coefficient to be updated.
}
}
\details{
This simple generic method is used for the numerical delta method.
}
\value{
The updated object is returned.
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, value) 
UseMethod("coef<-")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
