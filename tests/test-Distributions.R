context("Distributions")

test_that("Weights work appropriately:", {
})

test_that("Univariate Normal functions are ok:", {
    expect_equal(-dnorm(3, log = TRUE), devUN(3, params = c(0, 1)))
    expect_equal(-dnorm(-2, log = TRUE), devUN(-2, params = c(0, 1)))
    expect_equal(-dnorm(3, log = TRUE) - dnorm(2, log = TRUE),
                 devUN(c(2, 3), params = c(0, 1)))
    mu = 0
    sigma = 1
    nGrad = attr(numericDeriv(quote(dnorm(-2, mu, sigma, log = TRUE)),
                              c("mu", "sigma")), "gradient")
    expect_equal(as.numeric(round(-nGrad, 4)),
                 round(gradDevUN(x = -2, params = c(0, 1)), 4))
})

test_that("Multivariate Normal functions are ok:", {
})

test_that("Univariate Skew Normal functions are ok:", {
})

test_that("Multivariate Skew Normal functions are ok:", {
    params = list(xi = c(1, -1), Omega = diag(c(2.2, 1.3)), alpha = rnorm(2))

    # Define some objects to avoid typing params again
    myDmsn = function(x){
        -2 * sn::dmst(x, log = TRUE, xi = params$xi, Omega = params$Omega,
                 alpha = params$alpha, nu = 1e8)
    }
    snParams = paramList2VecMSN(params)

    # Density functions
    expect_equal(myDmsn(c(3, 1)), devMSN(matrix(c(3, 1), nrow = 1),
                                       params = snParams))
    expect_equal(c(myDmsn(c(2, -2)), myDmsn(c(-1, 3))),
                 devMSN(matrix(c(2, -2, -1, 3), nrow = 2, byrow = TRUE),
                        params = snParams))
    
    # Gradient functions
    snGrad = sn:::mst.pdev.grad(snParams, x = matrix(1), fixed.nu = 1e8,
                                y = matrix(c(2, 3), nrow = 1), w = 1)
    myGrad = gradDevMSN(x = matrix(c(2, 3), nrow = 1), snParams)
    # Just ensure they're close
    expect_less_than(max(abs(snGrad - myGrad)), 1e-8)
    y = matrix(rnorm(30), nrow = 10)
    snParams = paramList2VecMSN(list(xi = rnorm(3), Omega = diag(rnorm(3)^2),
                                     alpha = rnorm(3)))
    snGrad = sn:::mst.pdev.grad(snParams, x = matrix(1, nrow = 10),
                                y = y, w = rep(1, 10), fixed.nu = 1e8)
    myGrad = gradDevMSN(x = y, snParams)
    myGrad = colSums(myGrad)
    expect_less_than(max(abs(snGrad - myGrad)), 1e-8)
    # apply mst.pdev.grad and ensure it gives the same matrix
    yList = lapply(split(y, row(y)), matrix, nrow = 1)
    snVec = sapply(yList, sn:::mst.pdev.grad, param = snParams,
                   x = matrix(1, nrow = 1), w = 1, fixed.nu = 1e8)
    myGrad = gradDevMSN(x = y, snParams)
    expect_less_than(max(abs(t(snVec) - myGrad)), 1e-8)
    # Omega with off-diagonal
    Omega = diag(abs(rnorm(3)) + 1)
    Omega[upper.tri(Omega)] = runif(3, min = -1, max = 1)
    Omega[lower.tri(Omega)] = Omega[upper.tri(Omega)]
    snParams = paramList2VecMSN(list(xi = rnorm(3), Omega = Omega,
                                     alpha = rnorm(3)))
    y = matrix(rnorm(30), nrow = 10)
    snGrad = sn:::mst.pdev.grad(snParams, x = matrix(1, nrow = 10),
                                y = y, w = rep(1, 10), fixed.nu = 1e8)
    myGrad = gradDevMSN(x = y, snParams)
    myGrad = colSums(myGrad)
    expect_less_than(max(abs(snGrad - myGrad)), 1e-8)
})

test_that("Univariate t functions are ok:", {
})

test_that("Multivariate t functions are ok:", {
    params = list(xi = c(1, -1), Omega = diag(c(2.2, 1.3)), nu = 231)

    # Define some objects to avoid typing params again
    myDmt = function(x){
        -2 * sn::dmst(x, log = TRUE, xi = params$xi, Omega = params$Omega,
                 alpha = c(0, 0), nu = params$nu)
    }
    tParams = paramList2VecMT(params)

    # Density functions
    expect_equal(myDmt(c(3, 1)), devMT(matrix(c(3, 1), nrow = 1),
                                       params = tParams))
    expect_equal(c(myDmt(c(2, -2)), myDmt(c(-1, 3))),
                 devMT(matrix(c(2, -2, -1, 3), nrow = 2, byrow = TRUE),
                        params = tParams))
    
    # Gradient functions
    snGrad = sn:::mst.pdev.grad(tParams, x = matrix(1), symmetr = TRUE,
                                y = matrix(c(2, 3), nrow = 1), w = 1)
    myGrad = gradDevMT(x = matrix(c(2, 3), nrow = 1), tParams)
    # Just ensure they're close
    expect_less_than(max(abs(snGrad - myGrad)), 1e-8)
    y = matrix(rnorm(30), nrow = 10)
    tParams = paramList2VecMT(list(xi = rnorm(3), Omega = diag(rnorm(3)^2),
                                   nu = exp(rnorm(1))))
    snGrad = sn:::mst.pdev.grad(tParams, x = matrix(1, nrow = 10),
                                y = y, w = rep(1, 10), symmetr = TRUE)
    myGrad = gradDevMT(x = y, tParams)
    myGrad = colSums(myGrad)
    expect_less_than(max(abs(snGrad - myGrad)), 1e-8)
    # apply mst.pdev.grad and ensure it gives the same matrix
    yList = lapply(split(y, row(y)), matrix, nrow = 1)
    snVec = sapply(yList, sn:::mst.pdev.grad, param = tParams,
                   x = matrix(1, nrow = 1), w = 1, symmetr = TRUE)
    myGrad = gradDevMT(x = y, tParams)
    expect_less_than(max(abs(t(snVec) - myGrad)), 1e-8)
    # Omega with off-diagonal
    Omega = diag(abs(rnorm(3)) + 1)
    Omega[upper.tri(Omega)] = runif(3, min = -1, max = 1)
    Omega[lower.tri(Omega)] = Omega[upper.tri(Omega)]
    tParams = paramList2VecMT(list(xi = rnorm(3), Omega = Omega,
                                   nu = exp(rnorm(1))))
    y = matrix(rnorm(30), nrow = 10)
    snGrad = sn:::mst.pdev.grad(tParams, x = matrix(1, nrow = 10),
                                y = y, w = rep(1, 10), symmetr = TRUE)
    myGrad = gradDevMT(x = y, tParams)
    myGrad = colSums(myGrad)
    expect_less_than(max(abs(snGrad - myGrad)), 1e-8)
})

test_that("Univariate skew t functions are ok:", {
})

test_that("Multivariate skew t functions are ok:", {
    params = list(xi = c(1, -1), Omega = diag(c(2.2, 1.3)),
                  alpha = c(2, -1), nu = 231)

    # Define some objects to avoid typing params again
    myDmst = function(x){
        -2 * sn::dmst(x, log = TRUE, xi = params$xi, Omega = params$Omega,
                 alpha = params$alpha, nu = params$nu)
    }
    stParams = paramList2VecMST(params)

    # Density functions
    expect_equal(myDmst(c(3, 1)), devMST(matrix(c(3, 1), nrow = 1),
                                        params = stParams))
    expect_equal(c(myDmst(c(2, -2)), myDmst(c(-1, 3))),
                 devMST(matrix(c(2, -2, -1, 3), nrow = 2, byrow = TRUE),
                        params = stParams))
    
    # Gradient functions
    snGrad = sn:::mst.pdev.grad(stParams, x = matrix(1),
                                y = matrix(c(2, 3), nrow = 1), w = 1)
    myGrad = gradDevMST(x = matrix(c(2, 3), nrow = 1), stParams)
    # Just ensure they're close
    expect_less_than(max(abs(snGrad - myGrad)), 1e-8)
    y = matrix(rnorm(30), nrow = 10)
    stParams = paramList2VecMST(list(xi = rnorm(3), Omega = diag(rnorm(3)^2),
                                     alpha = rnorm(3), nu = exp(rnorm(1))))
    snGrad = sn:::mst.pdev.grad(stParams, x = matrix(1, nrow = 10),
                                y = y, w = rep(1, 10))
    myGrad = gradDevMST(x = y, stParams)
    myGrad = colSums(myGrad)
    expect_less_than(max(abs(snGrad - myGrad)), 1e-8)
    # apply mst.pdev.grad and ensure it gives the same matrix
    yList = lapply(split(y, row(y)), matrix, nrow = 1)
    snVec = sapply(yList, sn:::mst.pdev.grad, param = stParams,
                   x = matrix(1, nrow = 1), w = 1)
    myGrad = gradDevMST(x = y, stParams)
    expect_less_than(max(abs(t(snVec) - myGrad)), 1e-8)
    # Omega with off-diagonal
    Omega = diag(abs(rnorm(3)) + 1)
    Omega[upper.tri(Omega)] = runif(3, min = -1, max = 1)
    Omega[lower.tri(Omega)] = Omega[upper.tri(Omega)]
    stParams = paramList2VecMST(list(xi = rnorm(3), Omega = Omega,
                                     alpha = rnorm(3), nu = exp(rnorm(1))))
    y = matrix(rnorm(30), nrow = 10)
    snGrad = sn:::mst.pdev.grad(stParams, x = matrix(1, nrow = 10),
                                y = y, w = rep(1, 10))
    myGrad = gradDevMST(x = y, stParams)
    myGrad = colSums(myGrad)
    expect_less_than(max(abs(snGrad - myGrad)), 1e-8)
})