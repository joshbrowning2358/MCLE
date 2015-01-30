library(sn)
library(foreach)
library(ggplot2)
library(reshape)

setwd("~/Professional Files/Mines/Research/Robust Estimators/")

#First, examine for multivariate t, i.e. alpha=0
params = data.frame(nu=c(10^seq(0,4,.1),Inf))
params = merge(params, data.frame(alpha=seq(-5,5,.25)) )
p = .25
params$LQW = sapply(1:nrow(params), function(i){
  qs = qst(c((1-p)/2, p/2, 0.25), alpha=params$alpha[i], nu=params$nu[i] )
  return( - (qs[1] + qs[2] - 2*qs[3])/(qs[1]-qs[2]) )
} )
q = .75
params$RQW = sapply(1:nrow(params), function(i){
  qs = qst(c((1+q)/2, 1-q/2, 0.75), alpha=params$alpha[i], nu=params$nu[i] )
  return( (qs[1] + qs[2] - 2*qs[3])/(qs[1]-qs[2]) )
} )
params$SK2 = sapply(1:nrow(params), function(i){
  qs = qst(c(0.25,0.5,0.75), alpha=params$alpha[i], nu=params$nu[i] )
  return( (qs[1] + qs[3] - 2*qs[2])/(qs[3]-qs[1]) )
} )
toPlot = melt(params, id.vars=c("nu","alpha") )
save(toPlot, file="Results/quantile_estimator_to_nu.RData")
ggplot( toPlot, aes(x=value, y=alpha, color=nu, group=nu) ) +
  scale_color_continuous(trans="log10") +
  facet_wrap( ~ variable, scale="free" ) +
  geom_line()

ggplot( toPlot, aes(x=value, y=nu, color=alpha, group=alpha) ) +
  scale_y_log10() +
  facet_wrap( ~ variable, scale="free" ) +
  geom_line()
