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
##'   list object with three elements: beta, Omega, nu.
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with three elements: beta (a vector prodiving the 
##'   center), Omega (a matrix describing the dispersion), and nu (a numeric
##'   providing the heaviness of the tails).
##'   
##' @name multivariateT
NULL

##' @rdname multivariateT
devMT = function(x, params){
    l = length(params)
    ## Add in 0's for alpha:
    params = c(params[1:(l-1)], rep(0, ncol(x)), params[l])
    out = apply(x, 1, function(data){
        sn:::mst.pdev(y = matrix(data, nrow = 1), x = matrix(1), param = params)
    })
    return(out)
}

gradDevMT = function(x, params){
    l = length(params)
    ## Add in 0's for alpha:
    params = c(params[1:(l-1)], rep(0, ncol(x)), params[l])
    out = apply(x, 1, function(data){
        sn:::mst.pdev.grad(y = matrix(data, nrow = 1), x = matrix(1),
                           param = params)
    })
    return(out)
}

paramVec2ListMT = function(paramVec){
    # d + d*(d+1)/2 + 1 = length(paramVec)
    # d^2/2 + 3d/2 + 1 - length(paramVec) = 0
    # Quadratic formula to get:
    d = -3/2 + sqrt((3/2)^2 - 4 * 1/2 * (1-length(paramVec))) / 
        (2 * 1/2)
    out = sn:::optpar2dplist(paramVec, p = 1, d = d)$dp
    out$alpha = NULL
    return(out)
}

paramList2VecMT = function(paramList){
    l = length(paramList$mu)
    paramList$alpha = rep(0, l)
    paramList = paramList[c("mu", "sigma", "alpha", "nu")]
    out = sn:::dplist2optpar(paramList)
    out = out[c(1:(l + l*(l+1)/2), length(out))]
}