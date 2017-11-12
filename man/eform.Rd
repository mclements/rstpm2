\name{eform.stpm2}
\alias{eform}
\alias{eform.stpm2}
\title{
S3 method for to provide exponentiated coefficents with confidence intervals. 
}
\description{
S3 method for to provide exponentiated coefficents with confidence intervals. 
}
\usage{
eform(object, ...)
\method{eform}{stpm2}(object, parm, level = 0.95, method = c("Profile"), name = "exp(beta)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
regression object
}
  \item{parm}{
not currently used
}
  \item{level}{
significance level for the confidence interval
}
  \item{method}{
Currently only the profile method is available.
}
  \item{name}{
name for the fitted value
}
  \item{...}{
    other arguments
}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
%% \examples{
%% ##---- Should be DIRECTLY executable !! ----
%% ##-- ==>  Define data, use random,
%% ##--	or do  help(data=index)  for the standard data sets.

%% ## The function is currently defined as
%% function (object, parm, level = 0.95, method = c("Profile"), 
%%     name = "exp(beta)") 
%% {
%%     method <- match.arg(method)
%%     if (missing(parm)) 
%%         parm <- TRUE
%%     estfun <- switch(method, Profile = confint)
%%     val <- exp(cbind(coef = coef(object), estfun(object, level = level)))
%%     colnames(val) <- c(name, colnames(val)[-1])
%%     val[parm, ]
%%   }
%% }
%% % Add one or more standard keywords, see file 'KEYWORDS' in the
%% % R documentation directory.
%% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
