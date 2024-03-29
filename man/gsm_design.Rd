% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_gsm.R
\name{gsm_design}
\alias{gsm_design}
\title{Extract design information from an stpm2/gsm object and newdata
for use in C++}
\usage{
gsm_design(object, newdata, inflate = 100)
}
\arguments{
\item{object}{stpm2/gsm object}

\item{newdata}{list or data-frame used for evaluation}

\item{inflate}{double value to inflate minimum and maximum times for root finding}
}
\value{
list that can be read by `gsm ssim::read_gsm(SEX args)` in C++
}
\description{
Extract design information from an stpm2/gsm object and newdata
for use in C++
}
