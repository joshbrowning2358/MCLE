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
NULL

##' @rdname multivariateNormal

devMN = function(x, params, w = rep(1, nrow(x))){
    n = nrow(x)
    ## Modified from sn:::msn.pdev
    d <- ncol(x)
    n <- sum(w)
    p <- 1
    dp. <- sn:::optpar2dplist(params, d = ncol(x), p = 1)
    logL <- w * sn:::dmsn(x, matrix(1, nrow = NROW(x)) %*% dp.$beta,
                          dp.$Omega, c(0, 0), log = TRUE)
    dev <- (-2) * logL
    return(dev)
}

##' @rdname multivariateNormal
gradDevMN = function(x, params){
    ## Taken from sn
    p = 1
    d = ncol(x)
    mu <- matrix(params[1:(p * d)], p, d)
    D <- exp(-2 * params[(p * d + 1):(p * d + d)])
    A <- diag(d)
    i0 <- p * d + d * (d + 1)/2
    A[!lower.tri(A, diag = TRUE)] <- params[(p * d + d + 1):i0]
    Oinv <- t(A) %*% diag(D, d, d) %*% A
    ## Formulas taken from
    ## http://stats.stackexchange.com/questions/27436/how-to-take-derivative-of-multivariate-normal-density
    muMat = matrix(mu, nrow = nrow(x), ncol = length(mu), byrow = TRUE)
    gradMu = Oinv %*% t(x - muMat)
    ## Gradient with respect to sigma, but what about with respect to
    ## optimization parameters?
    gradSigma = apply(x - muMat, 1, function(diffVec){
        -1/2*(Oinv - Oinv %*% matrix(diffVec, ncol = 1) %*% diffVec %*% Oinv)
    })
}

##' @rdname multivariateNormal
paramVec2ListMN = function(paramVec){
    ## d + d(d+1)/2 = len
    ## d^2 + 3d - 2*len = 0
    d = (-3 + sqrt(9 + 8 * length(paramVec))) / 2
    ## Add in alpha elements:
    paramVec = c(paramVec, rep(0, d))
    out = sn:::optpar2dplist(paramVec, p = 1, d = d)$dp
    names(out) = c("mu", "sigma", "alpha")
    out = out[c("mu", "sigma")]
    return(out)
}

##' @rdname multivariateNormal
paramList2VecMN = function(paramList){
    d = length(paramList[[1]])
    paramList$alpha = rep(0, d)
    paramVec = sn:::dplist2optpar(paramList)
    paramVec = paramVec[1:(length(paramVec) - d)]
    return(paramVec)
}
