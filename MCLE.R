#Function should take data, a distribution, and alpha and return the fitted parameters.
MCLE <- function(data, distr, robAlpha, robFuncType){
#  allowed.distr = c("mvst", "mvsn", "mvt", "mvn") 
  allowed.distr = c("mvn") 
  if(!distr %in% allowed.distr)
    stop("distr must be in '", paste(allowed.distr, collapse="', '"), "'")
  if(robAlpha>=1 | robAlpha<=0)
    stop("robAlpha must be in (0,1)")
  
  devFunc = eval( parse( text=paste0(distr, ".dev") ) )
  start = do.call( eval( parse( text=paste0(distr, ".start") ) ), args=list(data=data) )
  k = do.call( eval( parse( text=paste0("get.k.",distr) ) ), args=list(x=data, robAlpha=robAlpha) )
  robustFunc = get.robust.func(k=k, robFuncType=robFuncType)[[1]]
  
  soln = nlminb(start=start, function(params){
    obsNLL = devFunc(data, params=params)
    robNLL = robustFunc(obsNLL)
    return(sum(robNLL))
  })
  
  convFunc = eval(parse( text=paste0(distr,".params2dp") ) )
  out = convFunc(soln$par)
  if("Rinv" %in% names(out)){
    out$Sigma = t(out$Rinv)%*%out$Rinv
    out$Rinv = NULL
  }
  return( out )
}