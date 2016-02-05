##' Exponential Threshold
##' 
##' These functions define an exponential threshold to be applied to the
##' log-likelihood values.
##' 
##' @param x A numeric vector of the log-likelihoods of the observations.
##' 
##' @name exponentialThreshold
NULL

##' @rdname exponentialThreshold
constraintET = function(x, k = 10){
    ifelse(x <= k, x, 2*k - k*exp(-x/k + 1))
}

##' @rdname exponentialThreshold
constraintDerivET = function(x, k = 10){
    ifelse(x <= k, 1, exp(-x/k + 1))
}