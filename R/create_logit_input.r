create_logit_input <- function(modelInputData) {
#creates a data matrix in model form that is consistent with the built in model forms
#inputs:
#modelInputData: data structure with information on variable choices and the alternative matrix (distances)
#outputs:
#dataCompile: a data matrix

#' create_logit_input
#'
#' Creates a data matrix that is consistent with the built in model forms
#'
#' @param modelInputData Data structure with information on variable choices and the alternative matrix (distances)
#' @keywords fish
#' @export
#' @examples
#' outputs:
#' dataCompile - a data matrix

x9 <- diag(modelInputData[["alts"]]); # makes matrix of choice possibilites
x8 <- matrix(x9, 1, modelInputData[["alts"]]*modelInputData[["alts"]])
x7 <- matrix(rep(x9,each=modelInputData[["obs"]]),nrow=modelInputData[["obs"]])
miles <- modelInputData[["zonalChoices"]]/modelInputData[["scaleszonal"]]
x <- cbind(x7, miles)
 
dataCompile <- x

return(dataCompile)

}