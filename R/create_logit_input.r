create_logit_input <- function(choice) {
#' create_logit_input
#'
#' Creates a data matrix that is consistent with the built in model forms
#'
#' @param choice A (number of observations) x 1 vector of chosen locations
#' @keywords fish
#' @export
#' @examples
#' outputs:
#' dataCompile - a data matrix

x9 <- diag(max(choice)); # makes matrix of choice possibilites
x8 <- matrix(x9, 1, max(choice)*max(choice))
x7 <- matrix(rep(x9,each=dim(choice)[1]),nrow=dim(choice)[1])

dataCompile <- x7

return(dataCompile)

}