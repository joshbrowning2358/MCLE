% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate_poisson.R
\name{univariatePoisson}
\alias{devPsn}
\alias{gradDevPsn}
\alias{paramList2VecPsn}
\alias{paramVec2ListPsn}
\alias{univariatePoisson}
\title{Univariate Poisson}
\usage{
devPsn(x, params)

gradDevPsn(x, params)

paramVec2ListPsn(paramVec)

paramList2VecPsn(paramList)
}
\arguments{
\item{x}{A numeric vector of observations.  As the Poisson distribution is
discrete, the x values should all be integers.}

\item{params}{A numeric value, providing the value of lambda.}

\item{paramVec}{A numeric vector (of length one) containing the value of lambda.}

\item{paramList}{A list with one element: lambda (a numeric value providing 
the mean).}
}
\description{
Functions defining the univariate Poisson distribution.
}
\details{
The deviance function returns the deviance of all the observations with the 
given parameters.  The gradDev function computes the gradient of this 
deviance.  The other two functions are useful for converting parameter 
vectors to lists and vice versa, and thus they return the list or vector.
}

