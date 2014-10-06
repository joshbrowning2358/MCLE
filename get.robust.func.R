get.robust.func <- function(k, robFuncType){
  #Data quality checks
  if(k<=0)
    stop("k must be positive")
  if(!robFuncType %in% c("bounded", "none"))
    stop("Invalid robFuncType!")
  
  if(robFuncType=="bounded"){
    robFunc = function(nll){
      ifelse(nll<k, nll, k+1-exp(-(nll-k)) )
    }
    robFuncGrad = function(nll){
      ifelse(nll>k,exp(-(nll-k)),1)
    }
  }
  if(robFuncType=="none"){
    robFunc = function(nll){nll}
    robFuncGrad = function(nll){1}
  }
  
  return( list(robFunc=robFunc, robFuncGrad=robFuncGrad) )
}