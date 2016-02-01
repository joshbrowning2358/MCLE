##' Multivariate Skew-t
##' 
##' Functions defining the multivariate skew t-distribution.  These should be 
##' used to construct a distribution class as defined in this package.
##' 
##' The deviance function returns the deviance of all the observations with the 
##' given parameters.  The gradDev function computes the gradient of this 
##' deviance.  The other two functions are useful for converting parameter 
##' vectors to lists and vice versa, and thus they return the list or vector.
##' 
##' @param x A numeric vector of observations.
##' @param params A vector, typically as created by paramList2Vec called on a 
##'   list object with four elements: beta, Omega, alpha, nu.
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with four elements: beta (a vector prodiving the
##'   center), Omega (a matrix describing the dispersion), alpha (a vector
##'   describing the skewness), and nu (a numeric providing the heaviness of the
##'   tails).
##'   
##' @name univariateNormal
NULL

##' @rdname univariateNormal
devMST = function(x, params){
    sn:::mst.pdev(x = matrix(1, nrow(x)), y = x, param = params,
                   w = rep(1, nrow(x)))
}

gradDevUN = function(x, params){
    sn:::mst.pdev.grad(x = matrix(1, nrow(x)), y = x, param = params,
                       w = rep(1, nrow(x)))
}

paramVec2ListUN = function(paramVec){
    # d + d*(d+1)/2 + d + 1 = length(paramVec)
    # d^2/2 + 5d/2 + 1 - length(paramVec) = 0
    # 
    d = -5/2 + sqrt((5/2)^2 - 4 * 1/2 * (1-length(paramVec))) / 
        (2 * 1/2)
    sn:::optpar2dplist(paramVec, p = 1, d = d)$dp
}

paramList2VecUN = function(paramList){
    sn:::dplist2optpar(paramList)
}