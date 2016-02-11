##' Multivariate Skew-t
##' 
##' Functions defining the multivariate skew t-distribution.
##' 
##' The deviance function returns the deviance of all the observations with the 
##' given parameters.  The gradDev function computes the gradient of this 
##' deviance.  The other two functions are useful for converting parameter 
##' vectors to lists and vice versa, and thus they return the list or vector.
##' 
##' @param x A numeric matrix of observations (one row per observation).
##' @param params A vector, typically as created by paramList2Vec called on a 
##'   list object with four elements: beta, Omega, alpha, nu.
##' @param paramVec A vector of the parameters of the multivariate distribution.
##' @param paramList A list with four elements: beta (a vector prodiving the 
##'   center), Omega (a matrix describing the dispersion), alpha (a vector 
##'   describing the skewness), and nu (a numeric providing the heaviness of the
##'   tails).
##'   
##' @name multivariateSkewT
NULL

##' @rdname multivariateSkewT
devMST = function(x, params, w = rep(1, nrow(x)), fixed.nu = NULL,
                  symmetr = FALSE){
    l = length(params)
    # Taken from sn:::mst.pdev:
    d <- ncol(x)
    p <- 1
    npar0 <- (p * d + d * (d + 1)/2)
    param1 <- c(params[1:npar0], if (symmetr) rep(0, d) else params[npar0 + 
        (1:d)], if (is.null(fixed.nu)) params[length(params)])
    dp <- paramVec2ListMST(param1)
    if (is.null(fixed.nu)){
        nu = dp$nu
    } else {
        nu = fixed.nu
    }
    logL <- w * sn:::dmst(x, matrix(1, nrow = nrow(x)) %*% dp$beta,
                          dp$Omega, dp$alpha, nu, log = TRUE)
    dev <- (-2) * logL
    return(dev)
}

gradDevMST = function(x, params, w = rep(1, nrow(x)), symmetr = FALSE, fixed.nu = FALSE){
    l = length(params)
    
    ## Taken from sn:::mst.pdev.grad
    d <- ncol(x)
    p <- 1
    beta <- matrix(params[1:(p * d)], p, d)
    D <- exp(-2 * params[(p * d + 1):(p * d + d)])
    A <- diag(d)
    i0 <- p * d + d * (d + 1)/2
    if (d > 1) 
        A[!lower.tri(A, diag = TRUE)] <- params[(p * d + d + 1):i0]
    if (symmetr){
        eta <- rep(0, d)
    } else {
        eta <- params[(i0 + 1):(i0 + d)]
    }
    if (is.null(fixed.nu)){
        nu <- exp(params[length(params)])
    } else {
        nu <- fixed.nu
    }
    Oinv <- t(A) %*% diag(D, d, d) %*% A
    u <- x - matrix(1, nrow = nrow(x)) %*% beta
    u.w <- u * w
    Q <- as.vector(rowSums((u %*% Oinv) * u.w))
    L <- as.vector(u.w %*% eta)
    if (nu < 10000){
        sf <- sqrt((nu + d)/(nu + Q))
    } else {
        sf <- sqrt((1 + d/nu)/(1 + Q/nu))
    }
    t. <- L * sf
    dlogft <- (-0.5) * sf^2
    dt.dL <- sf
    dt.dQ <- (-0.5) * L * sf/(Q + nu)
    logT. <- pt(t., nu + d, log.p = TRUE)
    dlogT. <- exp(dt(t., nu + d, log = TRUE) - logT.)
    stop("x isn't what you're calling x!  It's the matrix of 1's!!!")
#     Dbeta <- (-2 * t(x) %*% (u.w * dlogft) %*% Oinv - outer(as.vector(t(x) %*% 
#         (dlogT. * dt.dL * w)), eta) - 2 * t(x) %*% (dlogT. * 
#         dt.dQ * u.w) %*% Oinv)
    M1 = apply((u.w * dlogft), 1, function(vec){vec %*% Oinv})
    M2 = vecOuterMult(2 * x, (dlogT. * dt.dQ * u.w))
    M3 = vecOuterMult(x * (dlogT. * dt.dL * w), outer(rep(1, nrow(x)), eta))
    Dbeta <- mapply(function(m1, m2, m3){
        (m1 - m2) %*% Oinv - m3
    }, M1, M2, M3)
    Deta <- dlogT. * sf * u.w
    if (d > 1) {
        M <- lapply(vecOuterMult(u * dlogft + u * dlogT. * dt.dQ, u.w),
                    function(mat){
                        2 * (diag(D, d, d) %*% A %*% mat)
                    })
        DA <- lapply(M, function(x) x[!lower.tri(x, diag = TRUE)])
    } else {
        DA <- NULL
    }
    M <- lapply(vecOuterMult(u * dlogft + u * dlogT. * dt.dQ, u.w), function(mat)
        A %*% mat %*% t(A))
    if (d > 1){ 
        DD <- mapply(function(mat, weight){
            diag(mat) + 0.5 * weight/D
        }, mat = M, weight = w)
    } else {
        DD <- mapply(fun(mat, weight){
            mat + 0.5 * weight / D
        }, mat = M, weight = w)
    }
    stop("Keep adjusting here!!!")
    grad <- (-2) * c(Dbeta, DD * (-2 * D), DA, if (!symmetr) Deta)
    if (is.null(fixed.nu)) {
        df0 <- min(nu, 1e+08)
        if (df0 < 10000) {
            diff.digamma <- digamma((df0 + d)/2) - digamma(df0/2)
            log1Q <- log(1 + Q/df0)
        } else {
            diff.digamma <- log1p(d/df0)
            log1Q <- log1p(Q/df0)
        }
        dlogft.ddf <- 0.5 * (diff.digamma - d/df0 + (1 + d/df0) * 
            Q/((1 + Q/df0) * df0) - log1Q)
        eps <- 1e-04
        df1 <- df0 + eps
        if (df0 < 10000){
            sf1 <- sqrt((df1 + d)/(Q + df1))
        } else {
            sf1 <- sqrt((1 + d/df1)/(1 + Q/df1))
        }
        logT.eps <- pt(L * sf1, df1 + d, log.p = TRUE)
        dlogT.ddf <- (logT.eps - logT.)/eps
        Ddf <- sum((dlogft.ddf + dlogT.ddf) * w)
        grad <- c(grad, -2 * Ddf * df0)
    }
    return(grad)
}

paramVec2ListMST = function(paramVec){
    # d + d*(d+1)/2 + d + 1 = length(paramVec)
    # d^2/2 + 5d/2 + 1 - length(paramVec) = 0
    # Quadratic formula to get:
    d = -5/2 + sqrt((5/2)^2 - 4 * 1/2 * (1-length(paramVec))) / 
        (2 * 1/2)
    sn:::optpar2dplist(paramVec, p = 1, d = d)$dp
}

paramList2VecMST = function(paramList){
    sn:::dplist2optpar(paramList)
}