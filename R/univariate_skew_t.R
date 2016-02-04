##' Univariate Skew-t
##' 
##' Functions defining the univariate skew t-distribution.
##' 
##' The deviance function returns the deviance of all the observations with the 
##' given parameters.  The gradDev function computes the gradient of this 
##' deviance.  The other two functions are useful for converting parameter 
##' vectors to lists and vice versa, and thus they return the list or vector.
##' 
##' @param x A numeric vector of observations.
##' @param params A vector, typically as created by paramList2Vec called on a 
##'   list object with four elements: xi, omega, alpha, nu.
##' @param paramVec A vector of the parameters of the univariate distribution.
##' @param paramList A list with four elements: xi (a numeric value providing
##'   the center), omega (a numeric value describing the dispersion), alpha (a
##'   numeric value describing the skewness), and nu (a numeric providing the
##'   heaviness of the tails).
##'   
##' @name univariateSkewT
NULL

##' @rdname univariateSkewT
devUST = function(x, params){
    sn:::st.pdev(x = matrix(1, length(x)), y = as.matrix(x), dp = params,
                   w = rep(1, length(x)))
}

gradDevUST = function(x, params){
    sn:::st.pdev.gh(x = matrix(1, length(x)), y = as.matrix(x), dp = params,
                       w = rep(1, length(x)))
}

paramVec2ListUST = function(paramVec){
    out = list(xi = paramVec[[1]], omega = paramVec[[2]],
               alpha = paramVec[[3]], nu = paramVec[[4]])
    return(out)
}

paramList2VecUST = function(paramList){
    paramList = unlist(paramList)
}