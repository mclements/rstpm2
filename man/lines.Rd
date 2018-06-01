\name{lines.stpm2}
\alias{lines.stpm2}
\alias{lines.pstpm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S3 methods for lines
}
\description{
S3 methods for lines
}
\usage{
\method{lines}{stpm2}(x, newdata = NULL, type = "surv", col = 1, ci.col= "grey",
lty = par("lty"), ci = FALSE, rug = FALSE, var = NULL,
exposed = incrVar(var), times = NULL,
type.relsurv = c("excess", "total", "other"),
ratetable = survival::survexp.us, rmap, scale = 365.24, ...)
\method{lines}{pstpm2}(x, newdata = NULL, type = "surv", col = 1,
ci.col= "grey",
lty = par("lty"), ci = FALSE, rug = FALSE, var = NULL,
exposed = incrVar(var), times = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{stpm2} object}
  \item{newdata}{required list of new data. This defines the unexposed
    newdata (\emph{excluding} the event times).}
  \item{type}{specify the type of prediction
  }
  \item{col}{line colour}
  \item{lty}{line type}
  \item{ci.col}{confidence interval colour}
  \item{ci}{whether to plot the confidence interval band (default=TRUE)}
  \item{rug}{whether to add a rug plot of the event times to the current
    plot (default=TRUE)}
  \item{var}{specify the variable name or names for the
    exposed/unexposed (names are given as characters)}
  \item{exposed}{function that takes newdata and returns the exposed
    dataset. By default, this increments \code{var}}
  \item{times}{specifies the times. By default, this uses a span of the
                      observed times.}
  \item{type.relsurv}{type of predictions for relative survival models:
               either "excess", "total" or "other"}
  \item{scale}{scale to go from the days in the \code{ratetable} object
    to the analysis time used in the analysis}
  \item{rmap}{an optional list that maps data set names to the ratetable
          names. See \code{survexp}}
  \item{ratetable}{a table of event rates used in relative survival when
               \code{type.relsurv} is "total" or "other"}
  \item{\dots}{additional arguments (add to the \code{plot} command)}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
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

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
%% \examples{
%% ##---- Should be DIRECTLY executable !! ----
%% ##-- ==>  Define data, use random,
%% ##--	or do  help(data=index)  for the standard data sets.

%% }
%% % Add one or more standard keywords, see file 'KEYWORDS' in the
%% % R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
