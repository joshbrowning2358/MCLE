library(robustbase)

#Left quantile weight
#Breakdown point is p/2
LQW = function(data, p=0.25){
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(!is.null(dim(data)))
    stop("data must be a numeric vector!")
  qs = quantile(data, probs=c((1-p)/2, p/2, 0.25) )
  return( - (qs[1] + qs[2] - 2*qs[3])/(qs[1]-qs[2]) )
}

#Right quantile weight
#Breakdown point is (1-q)/2
RQW = function(data, q=0.75){
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(!is.null(dim(data)))
    stop("data must be a numeric vector!")
  qs = quantile(data, probs=c((1+q)/2, 1-q/2, 0.75) )
  return( (qs[1] + qs[2] - 2*qs[3])/(qs[1]-qs[2]) )
}

#Left medcouple
#Breakdown point is 0.125
LMC = function(data){
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(!is.null(dim(data)))
    stop("data must be a numeric vector!")
  -mc( data[data<median(data)] )
}

#Left medcouple
#Breakdown point is 0.125
RMC = function(data){
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(!is.null(dim(data)))
    stop("data must be a numeric vector!")
  mc( data[data>median(data)] )
}

#Peakedness
#Breakdown point is 0.125
P = function(data){
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(!is.null(dim(data)))
    stop("data must be a numeric vector!")
  qs = quantile(data, probs=c(0.875, 0.125, 0.75, 0.25))
  return( (qs[1] - qs[2])/(qs[3]-qs[4]) ) 
}