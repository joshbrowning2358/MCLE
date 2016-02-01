##' Define functions
##' 
##' @param distribution A character value, currently allowed to be one of 
##'   "normal", "student t", "skew normal" or "skew t".
##' @param dimension A numeric value specifying the dimensionality of the
##'   distribution.
##' 
##' @return An object containing the relevant functions.
##' 
##' @export
##' 

defineFunctions = function(distribution, dimension = 1){
    # Data Quality Checks
    stopifnot(distribution %in% c("normal", "student-t", "skew normal", "skew-t"))
    stopifnot(is.numeric(dimension) & dimension >= 1)
    
    # Extract relevant functions
    if(distribution == "normal" & dimension == 1){
        deviance = devUN
        gradDeviance = gradDevUN
        paramVec2List = paramVec2ListUN
        paramList2Vec = paramList2VecUN
    } else {
        stop("Not yet implemented!")
    }
    return(list(deviance = deviance, gradDeviance = gradDeviance,
                paramVec2List = paramVec2List, paramList2Vec = paramList2Vec))
}