##' Get Normal Distribution
##' 
##' @return A distribution object for the normal distribution.
##' 
##' @export
##' 

getNormalDistribution = function(){
    distribution(dev = devUN, gradDev = gradDevUN,
                 paramList2Vec = paramList2VecUN,
                 paramVec2List = paramVec2ListUN)
}