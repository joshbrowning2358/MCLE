##' Simulate k vs MAD
##' 
##' This function simulates random variables from a (univariate) random number 
##' generator and compares the MAD absolute deviation with quantiles of the 
##' deviance.  Such a comparison can provide insight into the appropriate choice
##' of k, as k should be chosen based on which deviance values are "extreme" and
##' which are "reasonable".
##' 
##' In this simulation, nSize observations are drawn from the passed 
##' distribution.  From this data, the MAD and the quantile(s) of the deviance 
##' are computed.  This process is repeated nRuns times to generate the expected
##' relationship between the quantiles of the deviance and the MAD.
##' 
##' @param rand The random number generator function, accepting the number of 
##'   random observations to generate as well as any parameters passed in ....
##' @param dev The deviance function, accepting the value of an observation as 
##'   well as any parameters passed in ...
##' @param quantiles The quantiles of the deviance which are desired for 
##'   reporting.
##' @param nSize The number of random samples to draw within each run.
##' @param nRuns The number of total runs to perform.
##' @param ... Additional vectors of parameters.  Each passed argument must be 
##'   named and must be of type vector.  For each run, a sample will be drawn 
##'   from each vector and this value will be passed to both rand and dev.  See
##'   examples.
##' @param Logical.  Should the MAD be plotted on a log scale?
##'   
##' @return
##' 
##' @example
##' dev = function(x, ...){-2 * dnorm(x, ..., log = TRUE)}
##' simulateKvsMAD(rnorm, dev, sd = 10^seq(-3, 3, .5), mean = -2:2)
##' simulateKvsMAD(rnorm, dev, sd = 10^seq(-3, 3, .5), mean = -2:2,
##'                quantiles = c(0.8, 0.9, 0.95, 0.99))
##'                
##' dev = function(x, ...){-2 * sn::dst(x, ..., log = TRUE)}
##' simulateKvsMAD(sn::rst, dev, xi = -5:5, omega = 10^seq(-3, 3, .5),
##'                alpha = -5:5, nu = 10^seq(0, 4, 1/3),
##'                quantiles = c(0.8, 0.9, 0.95, 0.99))
##'   
##' @export
##' 

simulateKvsMAD = function(rand, dev, quantiles = c(0.95), nSample = 1000, nRuns = 100, ..., logMAD = TRUE){
    distArgs = list(...)
    results = lapply(1:nRuns, function(dummy){
        # Add additional arguments to call
        if(length(distArgs) > 0){
            additionalArgs = sapply(1:length(distArgs), function(i){
                paste(names(distArgs)[i], "=", sample(distArgs[[i]], size = 1))
            })
            additionalArgs = paste(",", paste(additionalArgs, collapse = ", "))
        } else {
            additionalArgs = ""
        }
        
        # Build function calls
        runSampleText = paste0("rand(nSample", additionalArgs, ")")
        computeDevText = paste0("dev(sampleData", additionalArgs, ")")
        
        sampleData = eval(parse(text = runSampleText))
        sampleDev = eval(parse(text = computeDevText))
        q = quantile(sampleDev, probs = quantiles)
        MAD = mad(sampleData)
        out = data.frame(MAD = MAD, q)
        colnames(out)[2] = "deviance"
        out$quantile = quantiles
        return(out)
    })
    results = do.call("rbind", results)
    
    out = ggplot(results, aes(x = MAD, y = deviance, color = quantile)) +
        geom_point()
    if(logMAD)
        out = out + scale_x_log10()
    out
}