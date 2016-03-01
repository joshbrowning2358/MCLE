##' Univariate Poisson
##' 
##' Functions defining the univariate Poisson distribution.
##' 
##' The deviance function returns the deviance of all the observations with the 
##' given parameters.  The gradDev function computes the gradient of this 
##' deviance.  The other two functions are useful for converting parameter 
##' vectors to lists and vice versa, and thus they return the list or vector.
##' 
##' @param x A numeric vector of observations.  As the Poisson distribution is
##'   discrete, the x values should all be integers.
##' @param params A numeric value, providing the value of lambda.
##' @param paramVec A numeric vector (of length one) containing the value of lambda.
##' @param paramList A list with one element: lambda (a numeric value providing 
##'   the mean).
##'   
##' @name univariatePoisson
NULL

##' @rdname univariatePoisson
devPsn = function(x, params){
    dp = paramVec2ListPsn(params)
    -2*dpois(x, lambda = dp$lambda, log = TRUE)
}

##' @rdname univariatePoisson
gradDevPsn = function(x, params){
    dp = paramVec2ListPsn(params)
    -2*(x/dp$lambda - 1)
}

##' @rdname univariatePoisson
paramVec2ListPsn = function(paramVec){
    list(lambda = paramVec[[1]])
}

##' @rdname univariatePoisson
paramList2VecPsn = function(paramList){
    paramList = unlist(paramList)
}