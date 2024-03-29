% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exponentialThreshold.R
\name{exponentialThreshold}
\alias{constraintDerivET}
\alias{constraintET}
\alias{exponentialThreshold}
\title{Exponential Threshold}
\usage{
constraintET(x, k = 10)

constraintDerivET(x, k = 10)
}
\arguments{
\item{x}{A numeric vector of the log-likelihoods of the observations.}

\item{k}{The value of the deviance at which point deviance values become
adjusted.  For this particular function, (adjusted) deviance values will
never exceed 2*k}
}
\description{
These functions define an exponential threshold to be applied to the 
log-likelihood values.
}

