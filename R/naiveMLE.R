##' Naive MLE
##' 
##' This function computes a Maximum Likelihood Estimate by using optim() to 
##' optimize the joint log-likelihood.
##' 
##' @param data A matrix of data.  If univariate, data should have 1 column.  If
##'   multivariate, data should have one column for each dimension.
##' @param dist A distribution object (essentially a list with functions 
##'   defining the distribution, but see the distribution class defined in this 
##'   package).  See getNormalDistribution and the functions it loads for 
##'   examples of how to define this distribution object and the corresponding 
##'   functions.
##' @param initial Initial estimate for the parameters.  This should be a list 
##'   of parameters such that it could be passed to the paramList2Vec function 
##'   in the distribution object.
##' @param returnOptim Logical.  If FALSE (default) the function returns the 
##'   numeric estimates of the parameters, converted back into a list object. If
##'   TRUE, the function returns a list with the results of the call to optim (named "optim"
##'   in the output) as well as the parameters (named "solution")
##'   
##' @examples
##' \dontrun{
##' data = rnorm(100)
##' uniNorm = getNormalDistribution()
##' initial = list(mu = 0, sigma = 1)
##' naiveMLE(data, dist, initial)
##' naiveMLE(data, dist, initial, returnOptim = TRUE)
##' mean(data)
##' sd(data)
##' 
##' data = matrix(rnorm(200), nrow = 100)
##' distribution = normalMult
##' initial = list(mu = c(0, 0), sigma = diag(2))
##' naiveMLE(data, distribution, initial)
##' apply(data, mean)
##' cov(data)
##' }
##' 
##' @return See the returnOptim argument description.
##'   

naiveMLE = function(data, dist, initial, returnOptim = FALSE){
    ## Data Quality Checks
    if(!is(initial, "list"))
        stop("The 'initial' argument must be a list!")
    if(is.null(names(initial)))
        stop("The 'initial' argument must be a list with named values!")
    stopifnot(is.numeric(data))
    stopifnot(is(dist, "distribution"))
    
    initialVec = dist@paramList2Vec(initial)
    fn = function(par){
        dist@dev(x = data, params = par)
    }
    gr = function(par){
        dist$grad(x = data, params = par)
    }
    optimResult = optim(initialVec, fn = fn, gr = gr)
    soln = dist@paramVec2List(optimResult$par)
    if(!returnOptim){
        return(soln)
    } else {
        return(list(optim = optimResult, soln = soln))
    }
}