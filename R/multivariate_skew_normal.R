##' Multivariate Skew Normal
##' 
##' Functions defining the multivariate skew normal distribution.  These should
##' be used to construct a distribution class as defined in this package.
##' 
##' The deviance function returns the deviance of all the observations with the 
##' given parameters.  The gradDev function computes the gradient of this 
##' deviance.  The other two functions are useful for converting parameter 
##' vectors to lists and vice versa, and thus they return the list or vector.
##' 
##' Note: the multivariate skew normal distribution is implemented as a 
##' multivariate skew t-distribution with nu parameter equal to 1e8.  This leads
##' to minimal differences between the two distributions.
##' 
##' @param x A matrix of observations with one row per observation.
##' @param params A vector, typically as created by paramList2Vec called on a 
##'   list object with three elements: xi (a numeric vector) Omega (a matrix),
##'   and alpha (a numeric vector).
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with three elements: xi (a numeric vector), Omega
##'   (a numeric matrix) and alpha (a numeric vector).
##'   
##' @name multivariateSkewNormal
NULL

##' @rdname multivariateSkewNormal
devMSN = function(x, params){
    devMST(x = x, params = params, fixed.nu = 1e8)
}

##' @rdname multivariateSkewNormal
gradDevMSN = function(x, params){
    gradDevMST(x = x, params = params, fixed.nu = 1e8)
}

##' @rdname multivariateSkewNormal
paramVec2ListMSN = function(paramVec){
    ## d + d(d+1)/2 + d = len
    ## d^2 + 5d - 2*len = 0
    d = (-5 + sqrt(25 + 8 * length(paramVec))) / 2
    out = sn:::optpar2dplist(paramVec, p = 1, d = d)$dp
    names(out) = c("xi", "Omega", "alpha")
    return(out)
}

##' @rdname multivariateSkewNormal
paramList2VecMSN = function(paramList){
    stopifnot(names(paramList) %in% c("xi", "Omega", "alpha"))
    paramList$nu = 1e8
    paramVec = sn:::dplist2optpar(paramList)
    # Remove nu from vec:
    paramVec = paramVec[-length(paramVec)]
    return(paramVec)
}
