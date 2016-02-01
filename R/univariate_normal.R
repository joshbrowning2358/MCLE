##' Univariate Normal
##' 
##' Functions defining the univariate normal distribution.  These should be 
##' used to construct a distribution class as defined in this package.
##' 
##' The deviance function returns the deviance of all the observations with the 
##' given parameters.  The gradDev function computes the gradient of this 
##' deviance.  The other two functions are useful for converting parameter
##' vectors to lists and vice versa, and thus they return the list or vector.
##' 
##' @param x A numeric vector of observations.
##' @param params A vector, typically as created by paramList2Vec called on a 
##'   list object with two elements: mu and sigma (both numerics).
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with two elements: mean (a numeric) and sigma
##'   (a numeric).
##'   
##' @name univariateNormal
NULL

##' @rdname multivariateNormal
devUN = function(x, params){
    params = paramVec2ListUN(params)
    -sum(dnorm(x, mean = params$mu, sd = params$sigma, log = TRUE))
}

gradDevUN = function(x, params){
    params = paramVec2ListUN(params)
    df.dmu = -(x-params$mu)/params$sigma^2
    df.dsigma = 1/params$sigma - (x-params$mu)^2/params$sigma^3
    return(c(sum(df.dmu), sum(df.dsigma)))
}

paramVec2ListUN = function(paramVec){
    n = (-1 + sqrt(1 + 4 * length(paramVec))) / 2
    return(list(
        mu = paramVec[1:n],
        sigma = matrix(paramVec[(n+1):length(paramVec)], nrow = n)
    ))
}

paramList2VecUN = function(paramList){
    return(c(paramList$mu, as.numeric(paramList$sigma)))
}