##' Check Distribution
##'
##' @param object An S4 object.
##' 
##' @return TRUE if there are no errors, otherwise the error messages.
##' 

checkDistribution = function(object){
#     errors = character()
#     if(!c("x", "params") %in% args(object@dev)){
#         msg = "Function dev() must have 'x' and 'params' as arguments!"
#         errors = c(errors, msg)
#     }
#     if(!c("x", "params") %in% args(object@gradDev)){
#         msg = "Function gradDev() must have 'x' and 'params' as arguments!"
#         errors = c(errors, msg)
#     }
#     if("paramVec" != args(object@paramVec2List)){
#         msg = "Function paramVec2List() must have 'paramVec' as the only argument!"
#         errors = c(errors, msg)
#     }
#     if("paramList" != args(object@paramList2Vec)){
#         msg = "Function paramList2Vec() must have 'paramList' as the only argument!"
#         errors = c(errors, msg)
#     }
}

##' Distribution Class
##'
##' \describe{
##'     \item{\code{dev}:}{The function computing the deviance for the
##'     distribution.}
##'     \item{\code{gradDev}:}{The function computing the gradient of the
##'     deviance.}
##'     \item{\code{paramList2Vec}:}{A function converting a named list of the
##'     parameters to a vector for optimization.}
##'     \item{\code{paramVec2List}:}{A function converting a vector of the
##'     parameters to a named list for easier interpretation/output.}
##' }
##' 
##' @export distribution
##' 

distribution = setClass(Class = "distribution",
    slots = c(dev = "function", gradDev = "function",
              paramList2Vec = "function", paramVec2List = "function"),
    validity = checkDistribution,
    package = "MCLE")