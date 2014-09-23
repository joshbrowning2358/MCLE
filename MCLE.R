setwd("~/GitHub/MCLE/")
source("mvn.dev.R")
source("mvt.dev.R")
source("mvst.dev.R")
source("mvsn.dev.R")
source("get.robust.func.R")
library(robustbase)

#Function should take data, a distribution, and alpha and return the fitted parameters.
MCLE <- function(data, distr, robAlpha, robFuncType){
  allowed.distr = c("mvst", "mvsn", "mvt", "mvn") 
  if(!distr %in% allowed.distr)
    stop("distr must be in '", paste(allowed.distr, collapse="', '"), "'")
  if(robAlpha>=1 | robAlpha<=0)
    stop("robAlpha must be in (0,1)")
  
  devFunc = eval( parse( text=paste0("robust.",distr, ".dev") ) )
  start = do.call( eval( parse( text=paste0(distr, ".start") ) ), args=list(data=data) )
  
  soln = nlminb(start=start, function(params){
    sum(apply(data, 1, devFunc, params=params, robAlpha=robAlpha, robFuncType=robFuncType))
  }, control=list(trace=1)
  )
  params6 = c(14.5543,23.0313,-6.77108,379.037,90.6066,316.534,33.1582,-18.5518,326.433)
  devFunc(data[1,], params6, .01, "bounded")
  devFunc(data[1,], mvn.start(data), .01, "bounded")
  
  convFunc = eval(parse( text=paste0(distr,".params2dp") ) )
  out = convFunc(soln$par)
  if("Rinv" %in% names(out)){
    out$Sigma = t(out$Rinv)%*%out$Rinv
    out$Rinv = NULL
  }
  return( out )
}

data = matrix(rnorm(300), nrow=100)
MCLE( data, "mvn", robAlpha=.01, robFuncType="none" )
MCLE( data, "mvsn", robAlpha=.01, robFuncType="none" )
MCLE( data, "mvt", robAlpha=.01, robFuncType="none" )
MCLE( data, "mvst", robAlpha=.01, robFuncType="none" )
MCLE( data, "mvn", robAlpha=.01, robFuncType="bounded" )
MCLE( data, "mvsn", robAlpha=.01, robFuncType="bounded" )
MCLE( data, "mvt", robAlpha=.01, robFuncType="bounded" )
MCLE( data, "mvst", robAlpha=.01, robFuncType="bounded" )