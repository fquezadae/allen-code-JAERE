clogit <- function(starts3,dat,alts,otherdat) {
#' clogit
#'
#' Conditional logit likelihood
#'
#' @param starts3 Starting values
#' @param d Data matrix, see output from shiftSortX
#' @param alts Number of alternative choices in model
#' @keywords fish
#' @export
#' @examples
#' outputs:
#' ld - negative log likelihood

ld1 <- list()
starts3 <- as.matrix(starts3)
otherdat <- as.matrix(otherdat)

for(i in 1:dim(dat)[1])
{

betas1 <- c((starts3[1,]*otherdat[i,]), as.matrix(starts3[2,]))
betas <- t(as.matrix(betas1))

djz <- t(dat[i,3:dim(dat)[2]])

dj <- matrix(djz, nrow = alts, ncol = dim(betas)[2])

xb <- dj%*%t(betas)
xb <- xb - xb[1]
exb <- exp(xb)
ld1[[i]] <- (-log(t(exb)%*%(rep(1, alts))))

}

ld <- (-do.call("sum", ld1))

return(ld)

}