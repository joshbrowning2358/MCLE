##' Define functions
##' 
##' @param distribution A character value, currently allowed to be one of 
##'   "normal", "student t", "skew normal" or "skew t".
##' 
##' @return An object containing the relevant functions.
##' 
##' @export
##' 

defineFunctions = function(distribution){
    # Data Quality Checks
    stopifnot(distribution %in% c("normal", "student-t", "skew normal", "skew-t"))
    
    # Extract relevant functions
    if(distribution == "normal"){
        deviance = devMN
        gradDeviance = gradDevMN
        paramVec2List = paramVec2ListMN
        paramList2Vec = paramList2VecMN
    } else if(distribution == "skew normal"){
        deviance = devMSN
        gradDeviance = gradDevMSN
        paramVec2List = paramVec2ListMSN
        paramList2Vec = paramList2VecMSN
    } else if(distribution == "student-t"){
        deviance = devMT
        gradDeviance = gradDevMT
        paramVec2List = paramVec2ListMT
        paramList2Vec = paramList2VecMT
    } else if(distribution == "skew-t"){
        deviance = devMST
        gradDeviance = gradDevMST
        paramVec2List = paramVec2ListMST
        paramList2Vec = paramList2VecMST
    } else {
        stop("Not yet implemented!")
    }
    return(list(deviance = deviance, gradDeviance = gradDeviance,
                paramVec2List = paramVec2List, paramList2Vec = paramList2Vec))
}