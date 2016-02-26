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
##' dev = function(x, ...){-2 * sn::dmst(x, ..., log = TRUE)}
##' simulateKvsMAD(sn::rmst, dev, xi = list(c(-1, -1), c(0, 1), c(-2, 2)),
##'                Omega = list(diag(2) * .1, diag(2), diag(2)*10),
##'                alpha = list(c(0, 1)), nu = 10^seq(0, 4, 1/3),
##'                quantiles = c(0.8, 0.9, 0.95, 0.99))
##' 
##' @export
##' 
##' @import ggplot2
##' 

simulateKvsMAD = function(rand, dev, quantiles = c(0.95), nSample = 1000, nRuns = 100, ..., logMAD = TRUE){
    distArgs = list(...)
    results = lapply(1:nRuns, function(dummy){
        # Add additional arguments to call
        if(length(distArgs) > 0){
            sampledArgs = sapply(distArgs, sample, size = 1)
            # Assign values to the local environment.  Assignment is done this
            # way to allow sampling of complex objects like matrices.  Simply
            # pasting them into arguments is problematic.
            for(i in 1:length(distArgs)){
                assign(x = names(sampledArgs)[i], value = sampledArgs[[i]])
            }
            additionalArgs = lapply(names(sampledArgs), function(name){
                paste(",", name, "=", name)
            })
            additionalArgs = do.call(paste0, additionalArgs)
        } else {
            additionalArgs = ""
        }
        
        # Build function calls
        runSampleText = paste0("rand(nSample", additionalArgs, ")")
        computeDevText = paste0("dev(sampleData", additionalArgs, ")")
        
        # Sample data, compute deviance
        sampleData = eval(parse(text = runSampleText))
        sampleDev = eval(parse(text = computeDevText))
        dimension = NCOL(sampleData)
        
        # Calculate relevant statistics
        q = quantile(sampleDev, probs = quantiles)
        if(dimension > 1){
            MAD = apply(sampleData, 2, mad)
            # Collapse to a single number via multiplication
            MAD = prod(MAD)
        } else {
            MAD = mad(sampleData)
        }
        
        # Create output data
        out = data.frame(MAD = MAD, q)
        colnames(out)[2] = "deviance"
        out$quantile = quantiles
        # Add arguments to output data, but only if all arguments are simple
        # numerics
        if(length(sampledArgs) > 0 & all(lapply(sampledArgs, length) == 1)){
            toMerge = data.frame(matrix(sampledArgs, nrow = 1))
            colnames(toMerge) = names(distArgs)
            out = merge(toMerge, out)
        }
        return(out)
    })
    results = do.call("rbind", results)
    
    toPlot = reshape2::melt(results, id.vars = c("quantile", "deviance"))
    out = ggplot(toPlot, aes(x = value, y = deviance, color = quantile)) +
        geom_point() + facet_wrap( ~ variable)
    if(logMAD)
        out = out + scale_x_log10()
    out
}