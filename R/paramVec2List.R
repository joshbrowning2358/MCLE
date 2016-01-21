##' Deviance of a dataset
##' 
##' This function takes a vector or matrix of data and computes the deviance 
##' given some assumed distribution and provided parameters.
##' 
##' 
##' @result The deviance of the passed data.
##' 

dev = function(x, params, ...){
    UseMethod("dev")
}

dev.normalUN = function(x, params, ...){
    sum(-dt((x$data - params$mu)/params$sigma, df = params$df))
}

params = list(mu = 0, sigma = 1, df = 10)
tst = structure(list(data = rnorm(10)), class="student")
dens(tst, params)