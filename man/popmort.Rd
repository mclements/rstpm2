\name{popmort}
\Rdversion{1.1}
\alias{popmort}
\docType{data}
\title{
Background mortality rates for the colon dataset.
}
\description{
  Background mortality rates for the colon dataset.
}
\usage{data(popmort)}
\format{
  A data frame with 10600 observations on the following 5 variables.
  \describe{
    \item{\code{sex}}{Sex (1=male, 2=female)}
    \item{\code{prob}}{One year probability of survival}
    \item{\code{rate}}{All cause mortality rate}
    \item{\code{age}}{Age by single year of age through to age 105 years}
    \item{\code{year}}{Calendar period}
  }
}
%% \details{
%% %%  ~~ If necessary, more details than the __description__ above ~~
%% }
%% \source{
%% %%  ~~ reference to a publication or URL from which the data were obtained ~~
%% }
%% \references{
%% %%  ~~ possibly secondary sources and usages ~~
%% }
\examples{
data(popmort)
## maybe str(popmort) ; ...
}
\keyword{datasets}
