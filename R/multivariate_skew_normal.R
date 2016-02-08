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
##' @param x A matrix of observations with one row per observation.
##' @param params A vector, typically as created by paramList2Vec called on a 
##'   list object with three elements: mu (a numeric vector) sigma (a matrix),
##'   and alpha (a numeric vector).
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with three elements: mean (a numeric vector), sigma
##'   (a numeric matrix) and alpha (a numeric vector).
##'   
##' @name multivariateSkewNormal
NULL

##' @rdname multivariateSkewNormal

devMSN = function(x, params){
    n = nrow(x)
    out = apply(x, 1, function(data){
        sn:::msn.pdev(y = matrix(data, nrow = 1), param = params,
                      x = matrix(1))
    })
    return(out)
}

##' @rdname multivariateSkewNormal
gradDevMSN = function(x, params){
    stop("There is no built in function for this!!!")
}

##' @rdname multivariateSkewNormal
paramVec2ListMSN = function(paramVec){
    ## d + d(d+1)/2 + d = len
    ## d^2 + 5d - 2*len = 0
    d = (-5 + sqrt(25 + 8 * length(paramVec))) / 2
    out = sn:::optpar2dplist(paramVec, p = 1, d = d)$dp
    names(out) = c("mu", "sigma", "alpha")
    return(out)
}

##' @rdname multivariateSkewNormal
paramList2VecMSN = function(paramList){
    d = length(paramList[[1]])
    paramVec = sn:::dplist2optpar(paramList)
    return(paramVec)
}
