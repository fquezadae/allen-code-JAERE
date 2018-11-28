logit_zone_avgcat <- function(starts3,dat,otherdat,alts) {
#' logit_zone_avgcat
#'
#' Average catch conditional logit procedure
#'
#' @param starts3 Starting values with dimensions equal to number of alternatives + 1
#' @param dat Data matrix, see output from shiftSortX, alternatives with distance by column bind
#' @param otherdat Other data used in model (as list). Any number of interaction variables (e.g. vessel characteristics that affect 
#' how much disutility is suffered by traveling a greater distance) are allowed. However, the user must place these in `otherdat` 
#' as list objects named `intdat` respectively. Note the variables within `intdat` have no naming restrictions. 
#' Also note that `intdat` variables are dimension *(number of observations) x 1*, to be interacted with the distance to each alternative.
#' If there are no other data, the user can set `intdat` variables as ones with dimension *(number of observations) x 1*.
#' @param alts Number of alternative choices in model
#' @return ld - negative log likelihood
#' @export
#' @examples
#'

ld1 <- list()
intdat <- (otherdat$intdat)

starts3 <- as.matrix(starts3)
gridcoef <- as.matrix(starts3[1:(alts-1),])
intcoef <- as.matrix(starts3[alts:(alts-1+length(intdat)),])

for(i in 1:dim(dat)[1])
{

betas1 <- c(as.matrix(gridcoef), 
			t(as.matrix(do.call(rbind,lapply(intdat,`[`,i,))))%*%as.matrix(intcoef))
betas <- t(as.matrix(betas1))

djz <- t(dat[i,(alts+3):dim(dat)[2]])

dj <- matrix(djz, nrow = alts, ncol = dim(betas)[2])

xb <- dj%*%t(betas)
xb <- xb - xb[1]
exb <- exp(xb)
ld1[[i]] <- (-log(t(exb)%*%(rep(1, alts))))

}

ld <- (-do.call("sum", ld1))

return(ld)

}