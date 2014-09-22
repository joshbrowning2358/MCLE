robust.mvt.dev <- function(x, params, robAlpha, robFuncType){
  p = mvt.params2p(params=params)
  n = length(params)
  params = c(params[-n], rep(0,p), params[n]) #Append on alpha
  robust.mvst.dev(x=x, params=params, robAlpha=robAlpha, robFuncType=robFuncType)
}

mvt.start <- function(data){
  mu=colMeans(data)
  Sigma=var(data)
  R = chol(Sigma)
  Rinv = solve(R)
  return( mvt.dp2params(mu, Rinv, 10) )
}

mvt.params2p <- function(params){
  #length(params)=p(p+1)/2+p+1, use quadratic formula to find p
  p = (-3/2+sqrt(9/4+4/2*(length(params)-1)))
  if(p != floor(p))
    stop("p is not an integer!  params is not of the right length")
  return(p)
}

mvt.params2dp <- function(params){
  p = mvt.params2p(params)
  
  Rinv = matrix(0, nr=p, nc=p)
  Rinv[upper.tri(Rinv,diag=T)] = params[1:(p*(p+1)/2) + p]  
  return(list( mu=params[1:p]
              ,Rinv=Rinv
              ,nu=exp(params[1+p*(p+1)/2+p])) )
}

mvt.dp2params <- function(mu, Rinv, nu){
  p = length(mu)
  if(any(dim(Rinv)!=p))
    stop("Rinv must be pxp!")
  
  return( c(mu, Rinv[upper.tri(Rinv,diag=T)], log(nu) ) )
}