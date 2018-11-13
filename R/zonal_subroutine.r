zonal_subroutine <- function(modelInputData,dataCompile,betain,altcoeff,optimOpt,func,otherdat) {
#' zonal_subroutine
#'
#' Subroutine to run chosen discrete choice model
#'
#' @param modelInputData Information on model preference and parameters as LIST
#' @param dataCompile Data matrix from make_logit_input function as MATRIX
#' @param betain Initial parameter estimates for cost/distance
#' @param altcoeff Initial parameter estimates for revenue/location-specific covariates
#' @param optimOpt Optimization options [max function evaluations, max iterations, (reltol) tolerance of x]
#' @param func Name of likelihood function
#' @keywords fish
#' @export
#' @examples
#' outputs:
#' OutLogit - [outmat1 se1 tEPM2] (coefs, ses, tstats)
#' clogitoutput - optimization information
#' seoumat2 - ses
#' MCM - Model Comparison metrics
#' H1 - inverse hessian

errorExplain <- NULL
OutLogit <- NULL
clogitoutput <- NULL
seoutmat2 <- NULL
MCM <- NULL
H1 <- NULL
fr <- func #clogit

betain <- as.matrix(betain) #may be better to demand inputs a certain way? or change in here?

ab <- modelInputData[["alts"]] + max(dim(betain))
d <- shiftSortX(dataCompile,modelInputData[["choice"]],modelInputData[["alts"]],ab,modelInputData[["catch"]]/modelInputData[["scalescatch"]])

MCR <- 1
starts2 <- c(altcoeff,betain)

# dat <- d[,(modelInputData[["alts"]]+3):dim(d)[2]] #skips first zone choice

LL_start <- fr(starts2,d,modelInputData[["alts"]],otherdat)

if (is.null(LL_start) || is.nan(LL_start) || is.infinite(LL_start)) { #haven't checked what happens when error yet
   errorExplain <- "Initial function results bad (Nan, Inf, or undefined)"
   return("Initial function results bad (Nan, Inf, or undefined)") 
}

#############################################################################
mIter <- optimOpt[2]
MaxFunEvals <- optimOpt[1]
TolX  <-  optimOpt[3]

controlin <- list(maxit=mIter, reltol=TolX)

res <- tryCatch({

optim(starts2, fr, dat=d, alts=modelInputData[["alts"]], otherdat=otherdat, control = controlin, hessian = TRUE)

}, error = function(e) {

    errorExplain <- "Optimization error"
    return("Optimization error") 
	
})

q2 <- res[["par"]]
LL <- res[["value"]]
output <- list(counts = res[["counts"]], convergence = res[["convergence"]], mesage=res[["message"]])
H <- res[["hessian"]]
	
#Model Comparison metrocs (MCM)
param <- max(dim(as.matrix(starts2)))
obs <- dim(dataCompile)[1]
MCMAIC <- 2*param - 2*LL 

MCMAICc <- MCMAIC + (2*param*(param+1))/(obs - param - 1)

MCMBIC <- -2*LL + param*log(obs)

MCMPseudoR2 <- (LL_start-LL)/LL_start

MCM <- list(MCMAIC=MCMAIC, MCMAICc=MCMAICc, MCMBIC=MCMBIC, MCMPseudoR2=MCMPseudoR2)

if (is.null(H) == FALSE){
H1 <- solve(H)
diag2 <- diag(H1)
se2 <- sqrt(diag2)

outmat2 <- t(q2)
seoutmat2 <- t(se2)
clogitoutput <- output
tLogit <- t(outmat2/se2)
OutLogit <- cbind(t(outmat2), as.matrix(se2), (tLogit))
}

#############################################################################

return(list(errorExplain=errorExplain,OutLogit=OutLogit,clogitoutput=clogitoutput,seoutmat2=seoutmat2,MCM=MCM,H1=H1))

}