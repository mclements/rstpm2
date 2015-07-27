\name{colon}
\Rdversion{1.1}
\alias{colon}
\docType{data}
\title{
Colon cancer.
}
\description{
  Diagnoses of colon cancer.
}
% cat(mapply(function(name,label) sprintf("\\item{\\code{%s}}{%s}",name,label), names(colon), attr(colon, "var.labels")),sep="\n")
\usage{data(colon)}
\format{
  A data frame with 15564 observations on the following 13 variables.
  \describe{
\item{\code{sex}}{Sex (1=male, 2=female))}
\item{\code{age}}{Age at diagnosis}
\item{\code{stage}}{Clinical stage at diagnosis  (1=Unknown,
      2=Localised, 3=Regional, 4=Distant)}
\item{\code{mmdx}}{Month of diagnosis}
\item{\code{yydx}}{Year of diagnosis}
\item{\code{surv_mm}}{Survival time in months}
\item{\code{surv_yy}}{Survival time in years}
\item{\code{status}}{Vital status at last contact (1=Alive,
  2=Dead: cancer, 3=Dead; other, 4=Lost to follow-up)}
\item{\code{subsite}}{Anatomical subsite of tumour (1=Coecum and
  ascending, 2=Transverse, 3=Descending and sigmoid, 4=Other and NOS)}
\item{\code{year8594}}{Year of diagnosis (1=Diagnosed 75-84, 2=Diagnosed
  85-94)}
\item{\code{agegrp}}{Age in 4 categories (1=0-44, 2=45-59, 3=60-74, 4=75+)}
\item{\code{dx}}{Date of diagnosis}
\item{\code{exit}}{Date of exit}

  }
}
\details{
  Caution: there is a \code{colon} dataset in the \code{survival}
  package. We recommend using \code{data(colon,package="rstpm2")} to
  ensure the correct dataset is used.
}
%% \source{
%% %%  ~~ reference to a publication or URL from which the data were obtained ~~
%% }
%% \references{
%% %%  ~~ possibly secondary sources and usages ~~
%% }
\examples{
data(colon,package="rstpm2") # avoids name conflict with survival::colon
## maybe str(colon) ; ...
}
\keyword{datasets}
