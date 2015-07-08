\name{brcancer}
\Rdversion{1.1}
\alias{brcancer}
\docType{data}
\title{
German breast cancer data from Stata.
}
\description{
See \url{http://www.stata-press.com/data/r11/brcancer.dta}.
}
\usage{data(brcancer)}
\format{
  A data frame with 686 observations on the following 15 variables.
  \describe{
    \item{\code{id}}{a numeric vector}
    \item{\code{hormon}}{hormonal therapy}
    \item{\code{x1}}{age, years}
    \item{\code{x2}}{menopausal status}
    \item{\code{x3}}{tumour size, mm}
    \item{\code{x4}}{tumour grade}
    \item{\code{x5}}{number of positive nodes}
    \item{\code{x6}}{progesterone receptor, fmol}
    \item{\code{x7}}{estrogen receptor, fmol}
    \item{\code{rectime}}{recurrence free survival time, days}
    \item{\code{censrec}}{censoring indicator}
    \item{\code{x4a}}{tumour grade>=2}
    \item{\code{x4b}}{tumour grade==3}
    \item{\code{x5e}}{exp(-0.12*x5)}
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
data(brcancer)
## maybe str(brcancer) ; plot(brcancer) ...
}
\keyword{datasets}
