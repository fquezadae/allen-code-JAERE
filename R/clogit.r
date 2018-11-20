clogit <- function(starts3,dat,otherdat,alts) {
#' clogit
#'
#' Conditional logit likelihood
#'
#' @param starts3 Starting values
#' @param dat Data matrix, see output from shiftSortX, alternatives with distance by column bind
#' @param otherdat Other data used in model (as list)
#' @param alts Number of alternative choices in model
#' @return ld - negative log likelihood
#' @export
#' @examples
#'

ld1 <- list()
starts3 <- as.matrix(starts3)
pred_catch <- as.matrix(otherdat$predicted_catch)

for(i in 1:dim(dat)[1])
{

betas1 <- c((starts3[1,]*pred_catch[i,]), as.matrix(starts3[2,]))
betas <- t(as.matrix(betas1))
betas[,dim(betas)[2]] <- betas[,dim(betas)[2]]*otherdat$zi[i,]

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