##' Multivariate Skew-t
##' 
##' Functions defining the multivariate skew t-distribution.
##' 
##' The deviance function returns the deviance of all the observations with the 
##' given parameters.  The gradDev function computes the gradient of this 
##' deviance.  The other two functions are useful for converting parameter 
##' vectors to lists and vice versa, and thus they return the list or vector.
##' 
##' @param x A numeric matrix of observations (one row per observation).
##' @param params A vector, typically as created by paramList2Vec called on a 
##'   list object with four elements: beta, Omega, alpha, nu.
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with four elements: beta (a vector prodiving the 
##'   center), Omega (a matrix describing the dispersion), alpha (a vector 
##'   describing the skewness), and nu (a numeric providing the heaviness of the
##'   tails).
##'   
##' @name multivariateSkewT
NULL

##' @rdname multivariateSkewT
devMT = function(x, params){
    l = length(params)
    out = apply(x, 1, function(data){
        sn:::mst.pdev(y = matrix(data, nrow = 1), x = matrix(1), param = params)
    })
    return(out)
}

gradDevMT = function(x, params){
    l = length(params)
    out = apply(x, 1, function(data){
        sn:::mst.pdev.grad(y = matrix(data, nrow = 1), x = matrix(1),
                           param = params)
    })
    return(out)
}

paramVec2ListMST = function(paramVec){
    # d + d*(d+1)/2 + d + 1 = length(paramVec)
    # d^2/2 + 5d/2 + 1 - length(paramVec) = 0
    # Quadratic formula to get:
    d = -5/2 + sqrt((5/2)^2 - 4 * 1/2 * (1-length(paramVec))) / 
        (2 * 1/2)
    sn:::optpar2dplist(paramVec, p = 1, d = d)$dp
}

paramList2VecMST = function(paramList){
    sn:::dplist2optpar(paramList)
}