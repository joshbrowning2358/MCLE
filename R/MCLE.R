##' Maximum Constrained Likelihood Estimator
##' 
##' This function computes a "Maximum Constrained Likelihood Estimator" by using
##' optim() to optimize a constrained log-likelihood function.  The constraint 
##' is a function that is applied to the log-likelihood of each observation 
##' which, in the ideal case, prevents the log-likelihood from being too small 
##' and hence the observation too influential on the estimator.
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
##'   TRUE, the function returns a list with the results of the call to optim 
##'   (named "optim" in the output) as well as the parameters (named "solution")
##' @param constraintFunc A function taking a numeric vector of likelihood 
##'   values and returning the adjusted likelihood values.
##' @param constraintFuncDeriv A function providing the derivative of the
##'   constraint function.  This function should also accept a numeric vector
##'   and return a numeric vector of derivative values.
##'   
##' @examples
##' \dontrun{
##' data = rnorm(100)
##' uniNorm = list(dev = devUN, grad = gradDevUN,
##'                paramList2Vec = paramList2VecUN,
##'                paramVec2List = paramVec2ListUN)
##' initial = list(mu = 0, omega = 1)
##' MCLE(data, dist = uniNorm, initial)
##' MCLE(data, dist = uniNorm, initial, returnOptim = TRUE)
##' mean(data)
##' sd(data)
##' 
##' ust = list(dev = devUST, grad = gradDevUST,
##'            paramList2Vec = paramList2VecUST,
##'            paramVec2List = paramVec2ListUST)
##' initial = list(xi = 0, omega = 1, alpha = 0, nu = 100)
##' MCLE(data, dist = ust, initial)
##' 
##' data = matrix(rnorm(200), nrow = 100)
##' mst = list(dev = devMST, grad = gradDevMST,
##'            paramList2Vec = paramList2VecMST,
##'            paramVec2List = paramVec2ListMST)
##' initial = list(xi = c(0, 0), Omega = diag(c(1, 1)),
##'                alpha = c(0, 0), nu = 100)
##' MCLE(data, dist = mst, initial)
##' 
##' mt = list(dev = devMT, grad = gradDevMT,
##'           paramList2Vec = paramList2VecMT,
##'           paramVec2List = paramVec2ListMT)
##' initial = list(xi = c(0, 0), Omega = diag(c(1, 1)), nu = 100)
##' MCLE(data, dist = mt, initial)
##' 
##' msn = list(dev = devMSN, grad = gradDevMSN,
##'           paramList2Vec = paramList2VecMSN,
##'           paramVec2List = paramVec2ListMSN)
##' initial = list(xi = c(0, 0), Omega = diag(c(1, 1)), alpha = c(0, 0))
##' MCLE(data, dist = msn, initial)
##' 
##' mn = list(dev = devMN, grad = gradDevMN,
##'           paramList2Vec = paramList2VecMN,
##'           paramVec2List = paramVec2ListMN)
##' initial = list(xi = c(0, 0), Omega = diag(c(1, 1)), alpha = c(0, 0))
##' MCLE(data, dist = msn, initial)
##' 
##' data = rpois(30, lambda = 4.7)
##' data = c(data, 100)
##' dist = list(dev = devPsn, grad = gradDevPsn,
##'            paramList2Vec = paramList2VecPsn,
##'            paramVec2List = paramVec2ListPsn)
##' initial = list(lambda = 1)
##' MCLE(data, dist = dist, initial)
##' mean(data)
##' }
##' 
##' @return See the returnOptim argument description.
##' 
##' @export
##' 

MCLE = function(data, dist, initial, w = rep(1, NROW(data)), returnOptim = FALSE,
                constraintFunc = constraintET,
                constraintFuncDeriv = constraintDerivET){
    ## Data Quality Checks
    stopifnot(is.numeric(constraintFunc(c(1, 2, 3))))
    stopifnot(length(constraintFunc(c(1, 2, 3))) == 3)
    stopifnot(is.numeric(constraintFuncDeriv(c(1, 2, 3))))
    stopifnot(length(constraintFuncDeriv(c(1, 2, 3))) == 3)
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
        sum(constraintFunc(dist$dev(data, params = par, w = w)))
    }
    gr = function(par){
        # Chain rule: d/dtheta(f(L(x))) = f'(L(x)) * d/dtheta(L(x))
        constraintFuncDeriv(dist$dev(data, params = par, w = w)) %*%
            dist$grad(x = data, params = par, w = w)
    }
    optimResult = optim(initialVec, fn = fn, gr = gr)
    soln = dist$paramVec2List(optimResult$par)
    if(!returnOptim){
        return(soln)
    } else {
        return(list(optim = optimResult, soln = soln))
    }
}
