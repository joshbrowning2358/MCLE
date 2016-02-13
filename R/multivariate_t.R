##' Multivariate t
##' 
##' Functions defining the multivariate t-distribution.
##' 
##' The deviance function returns the deviance of all the observations with the 
##' given parameters.  The gradDev function computes the gradient of this 
##' deviance.  The other two functions are useful for converting parameter 
##' vectors to lists and vice versa, and thus they return the list or vector.
##' 
##' @param x A numeric matrix of observations (one row per observation).
##' @param params A vector, typically as created by paramList2Vec called on a 
##'   list object with three elements: xi, Omega, nu.
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with three elements: xi (a vector prodiving the 
##'   center), Omega (a matrix describing the dispersion), and nu (a numeric
##'   providing the heaviness of the tails).
##'   
##' @name multivariateT
NULL

##' @rdname multivariateT
devMT = function(x, params, w = rep(1, NROW(x))){
    devMST(x = x, params = params, w = w, symmetr = TRUE)
}

gradDevMT = function(x, params, w = rep(1, NROW(x))){
    gradDevMST(x = x, params = params, w = w, symmetr = TRUE)
}

paramVec2ListMT = function(paramVec){
    # d + d*(d+1)/2 + 1 = length(paramVec)
    # d^2/2 + 3d/2 + 1 - length(paramVec) = 0
    # Quadratic formula to get:
    d = -3/2 + sqrt((3/2)^2 - 4 * 1/2 * (1-length(paramVec))) / 
        (2 * 1/2)
    # Add alpha onto paramVec
    paramVec = c(paramVec[1:(length(paramVec)-1)], rep(0, d),
                 paramVec[length(paramVec)])
    out = sn:::optpar2dplist(paramVec, p = 1, d = d)$dp
    out$alpha = NULL
    return(out)
}

paramList2VecMT = function(paramList){
    stopifnot(names(paramList) %in% c("xi", "Omega", "nu"))
    l = length(paramList$xi)
    paramList$alpha = rep(0, l)
    paramList = paramList[c("xi", "Omega", "alpha", "nu")]
    out = sn:::dplist2optpar(paramList)
    out = out[c(1:(l + l*(l+1)/2), length(out))]
    return(out)
}