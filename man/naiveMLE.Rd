% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/naiveMLE.R
\name{naiveMLE}
\alias{naiveMLE}
\title{Naive MLE}
\usage{
naiveMLE(data, distribution, initial)
}
\arguments{
\item{data}{A matrix of data.  If univariate, data should have 1 column.  If
multivariate, data should have one column for each dimension.}

\item{distribution}{A list containing objects dev and grad.  Both are 
functions which take the data and parameters and calculate the deviance 
and gradient.  See univariate_normal.R for examples of how to define these
functions.}

\item{initial}{Initial estimate for the parameters.  This should be a list 
of parameters such that it could be passed to the functions defined by
distribution.}
}
\value{
A numeric vector giving the maximum likelihood estimates.
}
\description{
This function computes a Maximum Likelihood Estimate by using optim() to 
optimize the joint log-likelihood.
}
\examples{
\dontrun{
data = rnorm(100)
uniNorm = list(dev = devUN, grad = gradDevUN,
               paramList2Vec = paramList2VecUN,
               paramVec2List = paramVec2ListUN)
initial = list(mu = 0, sigma = 1)
naiveMLE(data, distribution = uniNorm, initial)
mean(data)
sd(data)

#data = matrix(rnorm(200), nrow = 100)
#distribution = normalMult
#initial = list(mu = c(0, 0), sigma = diag(2))
#naiveMLE(data, distribution, initial)
#apply(data, mean)
#cov(data)
}

}
