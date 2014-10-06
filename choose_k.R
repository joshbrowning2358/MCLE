choose_k <- function(data, robAlpha, alpha=rep(0,ncol(data)), nu=Inf){
  omega_hat = covMcd(data)
  xi = apply(data, 2, median)
  p = ncol(data)
  if(nu==Inf){
    #Constants dependent only on parameters:
    baseDev = -log(2) + p/2*log(2*pi) + 1/2*log(det(omega_hat))
    #Contribution to deviance based on observation:
    varDev = qf(1-robAlpha, df1=p, df2=Inf)
    #  - t( log( pnorm( (alpha/omega)%*%t(x-xi) ) ) ) # Ignore contribution to skewness, for now
    k = baseDev + varDev
  } else {
    k = 
  }
  return(k)
}