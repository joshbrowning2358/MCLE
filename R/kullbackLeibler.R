##' Kullback-Leibler divergence
##' 
##' A simple function to compute the Kullback-Leibler divergence between two 
##' (univariate or multivariate) probability distribution functions.
##' 
##' @param pdf1 A function taking a numeric value (univariate) or vector
##'   (multivariate) specifying the coordinate for where the density is
##'   required.  This function should be the true or known density.
##' @param pdf2 Same as pdf1, but for the estimated density.
##' @param dimension Integer.  The dimension of the space.  Assumed to be 1.
##' @param lower 
##' @param upper
##' 
##' @return The Kullback-Leibler divergence.
##' 
##' @export
##' 
##' @import cubature
##' 
##' @example
##' # Univariate examples
##' pdf1 = dnorm
##' pdf2 = function(x){dnorm(x, mean = 0.3)}
##' pdf3 = function(x){dnorm(x, mean = 2)}
##' kullbackLeibler(pdf1, pdf2)
##' kullbackLeibler(pdf1, pdf3)
##' 
##' # Multivariate examples
##' pdf1 = function(x){sn::dmsn(x, Omega = diag(2), alpha = 1:2)}
##' pdf2 = function(x){sn::dmst(x, Omega = diag(2), alpha = 1:2)}
##' pdf3 = function(x){sn::dmst(x, Omega = diag(c(1,2)), alpha = 1:2)}
##' kullbackLeibler(pdf1, pdf2, dimension = 2)
##' kullbackLeibler(pdf1, pdf3, dimension = 2)

kullbackLeibler = function(pdf1, pdf2, dimension = 1,
                           lower = rep(-10, dimension),
                           upper = rep(10, dimension)){
    # Data Quality Checks
    stopifnot(length(lower) == dimension)
    stopifnot(length(upper) == dimension)
    
    f = function(x){
        pdf1(x) * log(pdf1(x) / pdf2(x))
    }
    adaptIntegrate(f, lowerLimit = lower, upperLimit = upper)$integral
}