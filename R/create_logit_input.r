create_logit_input <- function(choice) {
    #' create_logit_input
    #'
    #' Creates a data matrix that is consistent with the built in model forms
    #'
    #' @param choice A (number of observations) x 1 vector of chosen locations
    #' @return dataCompile: a data matrix
    #' @export
    #' @examples
    #'
    
    dataCompile <- matrix(rep(diag(max(choice)), each = dim(choice)[1]),
        nrow = dim(choice)[1])
    
    return(dataCompile)
    
}
