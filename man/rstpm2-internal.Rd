\name{rstpm2-internal}
\Rdversion{1.1}
\alias{lhs}
\alias{rhs}
\alias{lhs<-}
\alias{rhs<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Internal functions for the rstpm2 package.
}
\description{
  Various utility functions used internally to the rstpm2 package.
}
\usage{
lhs(formula)
rhs(formula)
lhs(formula) <- value
rhs(formula) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    A formula
}
  \item{value}{
    A symbolic value to replace the current value. 
}
}
%% \details{
%%   %%  ~~ If necessary, more details than the description above ~~
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
%% ##---- Should be DIRECTLY executable !! ----
%% ##-- ==>  Define data, use random,
%% ##--	or do  help(data=index)  for the standard data sets.

%% ## The function is currently defined as
%% function (formula) 
%% if (length(formula) == 3) formula[[2]] else NULL
%% }
%% % Add one or more standard keywords, see file 'KEYWORDS' in the
%% % R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
