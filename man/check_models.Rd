% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_models.R
\name{check_models}
\alias{check_models}
\title{Check models}
\usage{
check_models(model = "m1")
}
\arguments{
\item{path}{Directory to save model outputs, defaults to /models.}

\item{model}{m0-m3, defaults to all.}
}
\description{
Having used \code{run_models()}, we can now evaluate the convergence of
the model and a variety of posterior predictive checks. These results
are output to the console.
}
