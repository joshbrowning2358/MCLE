#Function should take a distribution, parameters, outlier contamination counts, and some measure of outlier scale (as well as number of simulations).
#Then, a simulation should be performed comparing the MLE and MCLE(alpha=.05,.01,.005,.001)

#n specifies size of simulated dataset
#nSim specifies number of simulations to run
sim.MCLE <- function( distr, dp, n
    ,contamPerc=.1
    ,contamScale=10
    ,nSim=1000
    ,robAlpha=c(0.05, 0.01, 0.005, 0.001)){
  allowed.distr = c("mvn")
  stopifnot(distr %in% allowed.distr)
  stopifnot(contamPerc<=1 & contamPerc>=0)
  if(distr=="mvn"){
    stopifnot(is(dp,"list"))
    stopifnot(any(names(dp)==c("mu", "Rinv")))
  }
  
  #Generate seeds for simulations
  seeds = round(runif(nSim)*1000000)
  
  simResults = lapply(seeds, function(seed){
    set.seed(seed)
    data = do.call( paste0(distr,".sim"), args=list(n=n, dp=dp))
    out = list()
    for(rAlpha in robAlpha){
      out[[length(out)+1]] = MCLE( data, distr=distr, robAlpha=rAlpha, robFuncType="bounded")
      names(out)[length(out)] = paste("bounded:", rAlpha)
    }
    out[[length(out)+1]] = MCLE( data, distr=distr, robAlpha=0.01, robFuncType="none")
    names(out)[length(out)] = "MLE"
    #concatenate estimates into vectors
    out = lapply(out, function(subList){
      do.call("c", lapply( subList, function(x){as.numeric(x)} ) )
    } )
    out = data.frame( do.call("rbind", out) )
    out$MCLE = rownames(out)
    out$seed = seed
    return(out)
  } )
  simResults = do.call("rbind", simResults)
  
  if("Rinv" %in% names(dp)){
    dp$Sigma = solve(dp$Rinv%*%t(dp$Rinv))
    dp$Rinv = NULL
  }
  
  goal = do.call("c", lapply( dp, function(x){as.numeric(x)} ) )
  errResults = simResults[,1:length(goal)] -
    matrix(goal,nr=nrow(simResults), nc=length(goal), byrow=T)
  errResults = cbind(errResults, simResults[,-(1:length(goal))])
  MSE = ddply( errResults, "MCLE", function(df){
    df = df[,1:length(goal)]
    out = data.frame(mean = colMeans(df)
               ,median = apply(df, 2, median)
               ,max = apply(df, 2, max)
               ,min = apply(df, 2, min)
               ,q97.5 = apply(df, 2, quantile, probs=.975)
               ,q02.5 = apply(df, 2, quantile, probs=.025) )
    out = melt(out, id.vars=NULL)
    return(out)
  } )
  return(MSE)
}