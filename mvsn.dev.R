robust.mvsn.dev <- function(x, params, robAlpha, robFuncType){
  p = mvsn.params2p(params=params)
  params = c(params, Inf) #Append on alpha and nu
  robust.mvst.dev(x=x, params=params, robAlpha=robAlpha, robFuncType=robFuncType)
}

mvsn.start <- function(data){
  mu = colMeans(data)
  Sigma = var(data)
  R = chol(Sigma)
  Rinv = solve(R)
  return( mvsn.dp2params(mu, Rinv, rep(0,length(mu))) )
}

mvsn.params2p <- function(params){
  #length(params)=p(p+1)/2+2*p, use quadratic formula to find p
  p = (-5/2+sqrt(25/4+4/2*length(params)))
  if(p != floor(p))
    stop("p is not an integer!  params is not of the right length")
  return(p)
}

mvsn.params2dp <- function(params){
  p = mvsn.params2p(params)
  
  Rinv = matrix(0, nr=p, nc=p)
  Rinv[upper.tri(Rinv,diag=T)] = params[1:(p*(p+1)/2) + p]
  return(list( mu=params[1:p]
              ,Rinv=Rinv
              ,alpha=params[1:p+p+p*(p+1)/2]) )
}

mvsn.dp2params <- function(mu, Rinv, alpha){
  p = length(mu)
  if(any(dim(Rinv)!=p))
    stop("Rinv must be pxp!")
  if(length(alpha)!=p)
    stop("alpha must be of length p!")
  
  return( c(mu, Rinv[upper.tri(Rinv,diag=T)], alpha ) )
}