##' Naive MLE
##' 
##' This function computes a Maximum Likelihood Estimate by using optim() to 
##' optimize the joint log-likelihood.
##' 
##' @param data A matrix of data.  If univariate, data should have 1 column.  If
##'   multivariate, data should have one column for each dimension.
##' @param distribution A list containing objects dev and grad.  Both are 
##'   functions which take the data and parameters and calculate the deviance 
##'   and gradient.  See univariate_normal.R for examples of how to define these
##'   functions.
##' @param initial Initial estimate for the parameters.  This should be a list 
##'   of parameters such that it could be passed to the functions defined by
##'   distribution.
##' 
##' @examples
##' \dontrun{
##' data = rnorm(100)
##' uniNorm = distribution(devUN, gradDevUN, paramList2VecUN, paramVec2ListUN)
##' initial = list(mu = 0, sigma = 1)
##' naiveMLE(data, distribution, initial)
##' mean(data)
##' sd(data)
##' 
##' data = matrix(rnorm(200), nrow = 100)
##' distribution = normalMult
##' initial = list(mu = c(0, 0), sigma = diag(2))
##' naiveMLE(data, distribution, initial)
##' apply(data, mean)
##' cov(data)
##' }
##' 
##' @return A numeric vector giving the maximum likelihood estimates.
##'   

naiveMLE = function(data, distribution, initial){
    par = sapply(initial, as.numeric)
    if(is(par, "list"))
        par = do.call("c", par)
    fn = function(par){
        distribution@dev(x = data, params = par)
    }
    gr = function(par){
        distribution$grad(x = data, params = par)
    }
    optim(par, fn = fn, gr = gr)
}