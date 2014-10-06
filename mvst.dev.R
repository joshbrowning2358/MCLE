robust.mvst.dev <- function(x, params, robAlpha, robFuncType){
  p = mvst.params2p(params)
  dp = mvst.params2dp(params)
  xi = dp[[1]]
  Rinv = dp[[2]]
  alpha = dp[[3]]
  nu = dp[[4]]

  #Data quality checks
  if(min(diag(Rinv))<=0)
      stop("Rinv must be positive definite")
  if(nu<=0)
    stop("nu must be positive")

  if(robAlpha>=1 | robAlpha<=0 )
    stop("robAlpha must be in (0,1)")

  #omega = sqrt(diag(Omega)).  But, Omega=R^TR, so omega=diag(R)=1/diag(Rinv)
  omega = 1/diag(Rinv)
  
  #We want (x-xi)^T Omega^{-1} (x-xi)
  #If y=R^{-T}(x-xi), then y^T y=(x-xi)^T R^{-1} R^{-T} (x-xi), which is what we want.
  y = t( t(Rinv)%*%t(x-xi) )
  #Rinv is diagonal, so determinant is product of diagonal.
  #Omega=Rinv^{-T}%*%Rinv^{-1}, so det(Omega)=1/det(Rinv)^2
  detOmega = 1/prod(diag(Rinv))^2
  

  if(nu>100000){ #Approximately multivariate skew-normal, use that density
    #Begin penalizing when observations are larger than F(1-alpha)
    #For ease of use, ignore skewness
    #Then, quadratic form has p d.o.f. and numerator has nu d.o.f.
    detEst = det( covMcd(x)$cov )
    k = -log(2) + p/2*log(2*pi) + 1/2*detEst + qf(1-robAlpha, df1=p, df2=nu)
    robFunc = get.robust.func( k=k, robFuncType )[[1]]

    return( robFunc( -log(2) + p/2*log(2*pi) + 1/2*log(detOmega) + 1/2*rowSums(y^2)
      - t( log( pnorm( (alpha/omega)%*%t(x-xi) ) ) ) ) )
  } else {
    #Begin penalizing when observations are larger than F(1-alpha)
    #For ease of use, ignore skewness
    #Then, quadratic form has p d.o.f. and numerator has nu d.o.f.
    detEst = det( covMcd(x)$cov )
    k = -log(2) + p/2*log(2*pi) + 1/2*detEst + qf(1-robAlpha, df1=p, df2=nu)
    robFunc = get.robust.func( k=k, robFuncType )[[1]]
    warning("Appropriate choice of k has not been implemented!")
    warning("Need to vectorize this part of mvst.dev!")
    
    Q = t(y)%*%y
    return( robFunc( -log(2) - lgamma((nu+p)/2) + 1/2*log(detOmega) + p/2*log(pi*nu) +
        lgamma(nu/2) + (nu+p)/2*log(1+Q/nu) -
        log( pt( (t(alpha)/omega)%*%matrix(x-xi)*sqrt((nu+p)/(Q+nu)), nu+p ) ) ) )
  }
}

# robust.mvst.dev.grad <- function(x, params, robAlpha, robFuncType){
#   p = mvst.params2p(params)
#   dp = mvst.params2dp(params)
#   xi = dp[[1]]
#   Rinv = dp[[2]]
#   alpha = dp[[3]]
#   nu = dp[[4]]
# 
#   #Data quality checks
#   if(min(diag(Rinv))<=0)
#       stop("Rinv must be positive definite")
#   if(nu<=0)
#     stop("nu must be positive")
# 
#   if(robAlpha>=1 | robAlpha<=0 )
#     stop("robAlpha must be in (0,1)")
# 
#   #omega = sqrt(diag(Omega)).  But, Omega=R^TR, so omega=diag(R)=1/diag(Rinv)
#   omega = 1/diag(Rinv)
#   
#   y = Rinv%*%(x-xi)
#   #Rinv is diagonal, so determinant is product of diagonal.
#   #Omega=Rinv^{-T}%*%Rinv^{-1}, so det(Omega)=1/det(Rinv)^2
#   detOmega = 1/prod(diag(Rinv))^2
#   
#   #Begin penalizing when observations are larger than F(1-alpha)
#   #For ease of use, ignore skewness (i.e. assume mv normal or mv t for these purposes)
#   #Then, quadratic form has p d.o.f. and numerator has nu d.o.f.
#   robFunc = get.robust.func( k=qf(1-robAlpha, df1=p, df2=nu), robFuncType )[[1]]
#   robFuncGrad = get.robust.func( k=qf(1-robAlpha, df1=p, df2=nu), robFuncType )[[2]]
#     
#   if(nu>100000){ #Approximately multivariate skew-normal, use that density
#     if(p==1){
#       #Function is -log(2) + p/2*log(2*pi) + 1/2*log(detOmega) + robFunc( 1/2*y^2 - log( pnorm( alpha/omega*(x-xi) ) ) )
#       dXi = robFuncGrad( 1/2*y^2 - log( pnorm( alpha/omega*(x-xi) ) ) ) *
#         (y*Rinv - dnorm(alpha/omega*(x-xi))/pnorm(alpha/omega*(x-xi))*(-alpha/omega) )
#     } else {
#       #Function is -log(2) + p/2*log(2*pi) + 1/2*log(detOmega) + robFunc( 1/2*t(y)%*%y - log( pnorm( (t(alpha)/omega)%*%matrix(x-xi) ) ) ) )
#       dXi = robFuncGrad( 1/2*t(y)%*%y - log( pnorm( (t(alpha)/omega)%*%matrix(x-xi) ) ) ) *
#         (t(Rinv)%*%y - dnorm((t(alpha)/omega)%*%matrix(x-xi))/pnorm((t(alpha)/omega)%*%matrix(x-xi))*(-alpha/omega) )
#     }
#   } else {
#     Q = t(y)%*%y
#     if(p==1)
#       return( -log(2) + p/2*log(2*pi) + 1/2*log(detOmega) +
#           robFunc( 1/2*(x-xi)*OmegaInv*(x-xi) - log( pnorm( alpha/omega*(x-xi) ) ) ) )
#     return( -log(2) - lgamma((nu+p)/2) + 1/2*log(detOmega) + p/2*log(pi*nu) +
#         lgamma(nu/2) + robFunc( (nu+p)/2*log(1+Q/nu) -
#         log( pt( (t(alpha)/omega)%*%matrix(x-xi)*sqrt((nu+p)/(Q+nu)), nu+p ) ) ) )
#   }
# }

mvst.start <- function(data){
  mu=colMeans(data)
  Sigma=var(data)
  R = chol(Sigma)
  Rinv = solve(R)
  return( mvst.dp2params(mu, Rinv, rep(0,length(mu)), log(10)) )
}

mvst.params2p <- function(params){
  #length(params)=p(p+1)/2+2p+1, use quadratic formula to find p
  p = (-5/2+sqrt(25/4+4/2*(length(params)-1)))
  if(p!=floor(p))
    stop("p is not an integer!  params is not of the right length")
  return(p)
}

mvst.params2dp <- function(params){
  p = mvst.params2p(params)
  
  Rinv = matrix(0, nr=p, nc=p)
  Rinv[upper.tri(Rinv,diag=T)] = params[1:(p*(p+1)/2) + p]  
  return(list( mu=params[1:p]
              ,Rinv=Rinv
              ,alpha=params[1:p+p*(p+1)/2+p]
              ,nu=exp(params[1+p*(p+1)/2+2*p])) )
}

mvst.dp2params <- function(mu, Rinv, alpha, nu){
  p = length(mu)
  if(any(dim(Rinv)!=p))
    stop("Sigma must be pxp!")
  if(length(alpha)!=p)
    stop("alpha must be of length p!")
  if(length(nu)!=1)
    stop("nu must be of length 1!")
  if(nu<=0)
    stop("nu must be positive!")
  
  return( c(mu, Rinv[upper.tri(Rinv,diag=T)], alpha, log(nu) ) )
}