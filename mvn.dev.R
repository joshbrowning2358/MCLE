library(sn)
library(robustbase)
library(mnormt)

get.k.mvn <- function(x, robAlpha){
  p = mvn.params2p(params=params)
  Sigma = covMcd(x)$cov
  k = p/2*log(2*pi)+log(det(Sigma))+qchisq(1-robAlpha,df=p)
  return(k)
}

mvn.dev <- function(x, params){
  p = mvn.params2p(params=params)
  dp = mvn.params2dp(params)
  detSigma = 1/prod(diag(dp$Rinv))^2
  y = (x-matrix(dp$mu,nr=nrow(x),nc=ncol(x),byrow=T))%*%dp$Rinv
  return(p/2*log(2*pi) + log(detSigma) + apply(y*y,1,sum))
}

# mvn.dev.grad <- function(x, params, robAlpha, robFuncType){
#   p = mvn.params2p(params=params)
#   params = c(params, rep(0,p), Inf) #Append on alpha and nu
#   grad = robust.mvst.dev.grad(x=x, params=params, robAlpha=robAlpha, robFuncType=robFuncType)
#   return( grad[1:(p+p*(p+1)/2)] )
# }

mvn.start <- function(data){
  if(is.null(dim(data)))
    data = matrix(data,ncol=1)
  mu=colMeans(data)
  Sigma=covMcd(data)$cov
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

mvn.sim <- function(n, dp){
  mu = dp$mu
  Rinv = dp$Rinv
  Sigma = solve(Rinv%*%t(Rinv))
  rmnorm(n=n, mean=mu, varcov=Sigma)
}