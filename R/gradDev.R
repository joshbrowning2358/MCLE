##' Gradient of the Deviance
##' 
##' This function takes a vector or matrix of data and computes the gradient of
##' the deviance with respect to the parameters of a provided distribution and
##' parameter values.
##' 
##' @param x The object containing the data.  This should be an S3 object with 
##'   one element called "data".
##' @param params A numeric vector containing the parameters of the 
##'   distribution.
##'   
##' @result The gradient of the deviance for the passed data.
##'   

gradDev = function(x, params, ...){
    UseMethod("gradDev")
}

gradDev.normalUN = function(x, params){
    params = paramVec2List(params)
    df.dmu = -(x-params$mu)/params$sigma^2
    df.dsigma = 1/params$sigma - (x-params$mu)^2/params$sigma^3
    return(c(sum(df.dmu), sum(df.dsigma)))
}

params = list(mu = 0, sigma = 1, df = 10)
tst = structure(list(data = rnorm(10)), class="student")
dens(tst, params)