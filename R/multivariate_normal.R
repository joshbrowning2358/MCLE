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
##'   list object with two elements: mu (a numeric vector) and sigma (a matrix).
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with two elements: mean (a numeric vector) and sigma
##'   (a numeric matrix).
##'   
##' @name multivariateNormal
##' @importFrom mvtnorm dmvnorm
NULL

##' @rdname multivariateNormal

devMN = function(x, params){
    params = paramVec2ListMN(params)
    if(is.numeric(x)){
        if(length(params$mu) != 1)
            stop("One-dimensional data but not one-dimensional mu!")
        else
            x = matrix(x, ncol = 1)
    }
    -sum(mvtnorm::dmvnorm(x, mean = params$mu, sigma = params$sigma, log = TRUE))
}

##' @rdname multivariateNormal
gradDevMN = function(x, params){
    params = paramVec2List(params)
    df.dmu = -(x-params$mu)/params$sigma^2
    df.dsigma = 1/params$sigma - (x-params$mu)^2/params$sigma^3
    return(c(sum(df.dmu), sum(df.dsigma)))
}

##' @rdname multivariateNormal
paramVec2ListMN = function(paramVec){
    n = (-1 + sqrt(1 + 4 * length(paramVec))) / 2
    return(list(
        mu = paramVec[1:n],
        sigma = matrix(paramVec[(n+1):length(paramVec)], nrow = n)
    ))
}

##' @rdname multivariateNormal
paramList2VecMN = function(paramList){
    paramList$alpha = rep(0, length(paramList[[1]]))
    paramVec = sn:::msn.cp2cp(paramList)
    return(c(paramList$mu, as.numeric(paramList$sigma)))
}
