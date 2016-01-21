context("Distributions")

test_that("Univariate Normal functions are ok:", {
    expect_equal(-dnorm(3, log = TRUE), devNormalUni(3, params = c(0, 1)))
    expect_equal(-dnorm(-2, log = TRUE), devNormalUni(-2, params = c(0, 1)))
    expect_equal(-dnorm(3, log = TRUE) - dnorm(2, log = TRUE),
                 devNormalUni(c(2, 3), params = c(0, 1)))
    mu = 0
    sigma = 1
    nGrad = attr(numericDeriv(quote(dnorm(-2, mu, sigma, log = TRUE)),
                              c("mu", "sigma")), "gradient")
    expect_equal(as.numeric(round(-nGrad, 4)),
                 round(gradDevNormalUni(x = -2, params = c(0, 1)), 4))
})