##' Multivariate Normal
##' 
##' Functions defining the multivariate normal distribution.  These should be 
##' used to construct a distribution class as defined in this package.
##' 
##' The deviance function returns the deviance of all the observations with the 
##' given parameters.  The gradDev function computes the gradient of this 
##' deviance.  The other two functions are useful for converting parameter
##' vectors to lists and vice versa, and thus they return the list or vector.
##' 
##' @param x A matrix of observations with one row per observation.
##' @param params A vector, typically as created by paramList2Vec called on a 
##'   list object with two elements: xi (a numeric vector) and Omega (a matrix).
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with two elements: xi (a numeric vector) and Omega
##'   (a numeric matrix).
##'   
##' @name multivariateNormal
NULL

##' @rdname multivariateNormal

devMN = function(x, params, w = rep(1, nrow(x))){
    devMST(x = x, params = params, w = w, symmetr = TRUE, fixed.nu = 1e8)
}

##' @rdname multivariateNormal
gradDevMN = function(x, params){
    gradDevMST(x = x, params = params, w = w, symmetr = TRUE, fixed.nu = 1e8)
}

##' @rdname multivariateNormal
paramVec2ListMN = function(paramVec){
    ## d + d(d+1)/2 = len
    ## d^2 + 3d - 2*len = 0
    d = (-3 + sqrt(9 + 8 * length(paramVec))) / 2
    ## Add in alpha elements:
    paramVec = c(paramVec, rep(0, d))
    out = sn:::optpar2dplist(paramVec, p = 1, d = d)$dp
    names(out) = c("xi", "Omega", "alpha")
    out = out[c("xi", "Omega")]
    return(out)
}

##' @rdname multivariateNormal
paramList2VecMN = function(paramList){
    stopifnot(names(paramList) == c("xi", "Omega"))
    d = length(paramList[[1]])
    paramList$alpha = rep(0, d)
    paramVec = sn:::dplist2optpar(paramList)
    paramVec = paramVec[1:(length(paramVec) - d)]
    return(paramVec)
}
