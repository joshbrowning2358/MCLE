robust.mvn.dev <- function(x, params, robAlpha, robFuncType){
  p = mvn.params2p(params=params)
  params = c(params, rep(0,p), Inf) #Append on alpha and nu
  robust.mvst.dev(x=x, params=params, robAlpha=robAlpha, robFuncType=robFuncType)
}

robust.mvn.dev.grad <- function(x, params, robAlpha, robFuncType){
  p = mvn.params2p(params=params)
  params = c(params, rep(0,p), Inf) #Append on alpha and nu
  grad = robust.mvst.dev.grad(x=x, params=params, robAlpha=robAlpha, robFuncType=robFuncType)
  return( grad[1:(p+p*(p+1)/2)] )
}

mvn.start <- function(data){
  mu=colMeans(data)
  Sigma=var(data)
  R = chol(Sigma)
  Rinv = solve(R)
  return( mvn.dp2params(mu, Rinv ) )
}

mvn.params2p <- function(params){
  #length(params)=p(p+1)/2+p, use quadratic formula to find p
  p = (-3/2+sqrt(9/4+4*1/2*length(params)))
  if(p != floor(p))
    stop("p is not an integer!  params is not of the right length")
  return(p)
}

mvn.params2dp <- function(params){
  p = mvn.params2p(params)
  
  Rinv = matrix(0, nr=p, nc=p)
  Rinv[upper.tri(Rinv,diag=T)] = params[1:(p*(p+1)/2) + p]
  return(list( mu=params[1:p]
              ,Rinv=Rinv ) )
}

mvn.dp2params <- function(mu, Rinv){
  p = length(mu)
  if(any(dim(Rinv)!=p))
    stop("Rinv must be pxp!")
  
  return( c(mu, Rinv[upper.tri(Rinv,diag=T)] ) )
}