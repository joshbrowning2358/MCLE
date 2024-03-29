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
##' uniNorm = list(dev = devUN, grad = gradDevUN,
##'                paramList2Vec = paramList2VecUN,
##'                paramVec2List = paramVec2ListUN)
##' initial = list(mu = 0, sigma = 1)
##' naiveMLE(data, dist = uniNorm, initial)
##' naiveMLE(data, dist = uniNorm, initial, returnOptim = TRUE)
##' mean(data)
##' sd(data)
##' 
##' ust = list(dev = devUST, grad = gradDevUST,
##'            paramList2Vec = paramList2VecUST,
##'            paramVec2List = paramVec2ListUST)
##' initial = list(xi = 0, omega = 1, alpha = 0, nu = 100)
##' naiveMLE(data, dist = ust, initial)
##' 
##' data = matrix(rnorm(200), nrow = 100)
##' mst = list(dev = devMST, grad = gradDevMST,
##'            paramList2Vec = paramList2VecMST,
##'            paramVec2List = paramVec2ListMST)
##' initial = list(beta = c(0, 0), Omega = diag(c(1, 1)),
##'                alpha = c(0, 0), nu = 100)
##' naiveMLE(data, dist = mst, initial)
##' naiveMLE(data, dist = mst, initial, returnOptim = TRUE)
##' 
##' data = rpois(30, lambda = 4.7)
##' data = c(data, 100)
##' dist = list(dev = devPsn, grad = gradDevPsn,
##'            paramList2Vec = paramList2VecPsn,
##'            paramVec2List = paramVec2ListPsn)
##' initial = list(lambda = 1)
##' naiveMLE(data, dist = dist, initial)
##' mean(data)
##' }
##' 
##' @return See the returnOptim argument description.
##' 
##' @export
##'   

naiveMLE = function(data, dist, initial, returnOptim = FALSE){
    ## Data Quality Checks
    if(!is(initial, "list"))
        stop("The 'initial' argument must be a list!")
    if(is.null(names(initial)))
        stop("The 'initial' argument must be a list with named values!")
    stopifnot(is.numeric(data))
    stopifnot(is(dist, "list"))
    stopifnot(c("dev", "grad", "paramVec2List", "paramList2Vec") %in%
                  names(dist))
    
    initialVec = dist$paramList2Vec(initial)
    fn = function(par){
        sum(dist$dev(x = data, params = par))
    }
    gr = function(par){
        colSums(dist$grad(x = data, params = par))
    }
    optimResult = optim(initialVec, fn = fn, gr = gr)
    soln = dist$paramVec2List(optimResult$par)
    if(!returnOptim){
        return(soln)
    } else {
        return(list(optim = optimResult, soln = soln))
    }
}
