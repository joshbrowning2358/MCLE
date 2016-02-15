##' Vectorized Outer Multiplication
##' 
##' When computing densities and gradients, we occassionally need to compute the
##' outer product between two vectors.  Often, these vectors are stored in
##' matrices and thus cannot be simply multiplied together.  This function
##' performs this multiplication for two matrices X and Y.
##' 
##' @param X The first matrix.  Each row contains a (row) vector x, and the
##'   final result is the outer product: t(x) times y.
##' @param Y The second matrix.  Each row contains a (row) vector y, and the
##'   final result is the outer product t(x) times y.
##' 
##' @return A list containing the results of the various matrix products.
##' 
##' @examples 
##' X = matrix(1:10, nrow = 5)
##' Y = matrix(1, nrow = 5, ncol = 2)
##' vecOuterMult(X, Y)
##' 

vecOuterMult = function(X, Y){
    Map(tcrossprod, split(X, row(X)), split(Y, row(Y)))
}