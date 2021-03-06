% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_models.R
\name{format_data}
\alias{format_data}
\title{Format data for analysis}
\usage{
format_data(model, years = c(2013:2016), threshold = 0.2,
  subset = "presence", lkj_prior = 25, ...)
}
\arguments{
\item{years}{restrict analysis by years.}

\item{threshold}{restrict analysis to some proportion of total abundance.}

\item{subset}{subset species by presence or abundance.}

\item{lkj_prior}{shape value for covariance matrices in joint models.}
}
\description{
Takes the preloaded Pinnacle vegetation dataset and prepares it for
analysis in Stan. The Stan models are written such that the are able
to take largely the same data structure, with the exception of M0.
}
