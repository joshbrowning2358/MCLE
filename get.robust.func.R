get.robust.func <- function(k, robFuncType){
  #Data quality checks
  if(k<=0)
    stop("k must be positive")
  if(!robFuncType %in% c("bounded", "none"))
    stop("Invalid robFuncType!")
  
  if(robFuncType=="bounded")
    return( function(nll){
      if(nll>k)
        return(2*k-k*exp(-nll/k+1))
      else
        return(nll)
    } )
  if(robFuncType=="none")
    return( function(nll){nll})
}