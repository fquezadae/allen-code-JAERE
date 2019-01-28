discretefish_subroutine <- function(catch,choice,distance,otherdat,initparams,optimOpt,func,methodname) {
#' discretefish_subroutine
#'
#' Subroutine to run chosen discrete choice model
#'
#' @param catch Data corresponding to actual zonal choice
#' @param choice Data corresponding to actual catch
#' @param distance Data corresponding to distance
#' @param otherdat Other data (as a list)
#' @param initparams Initial parameter estimates for revenue/location-specific covariates then cost/distance
#' @param optimOpt Optimization options [max function evaluations, max iterations, (reltol) tolerance of x]
#' @param func Name of likelihood function
#' @return
#' OutLogit - [outmat1 se1 tEPM2] (coefs, ses, tstats) \cr 
#' optoutput - optimization information \cr 
#' seoumat2 - ses \cr 
#' MCM - Model Comparison metrics \cr 
#' H1 - inverse hessian \cr 
#' @export
#' @examples
#'

errorExplain <- NULL
OutLogit <- NULL
optoutput <- NULL
seoutmat2 <- NULL
MCM <- NULL
H1 <- NULL
fr <- func #e.g. clogit

ab <- max(choice) + 1 #no interactions in create_logit_input - interact distances in likelihood function instead
dataCompile <- create_logit_input(choice)

d <- shiftSortX(dataCompile,choice,catch,distance,max(choice),ab)

starts2 <- initparams

LL_start <- fr(starts2,d,otherdat,max(choice))

if (is.null(LL_start) || is.nan(LL_start) || is.infinite(LL_start)) { #haven't checked what happens when error yet
   errorExplain <- "Initial function results bad (Nan, Inf, or undefined), check 'ldglobalcheck'"
   return("Initial function results bad (Nan, Inf, or undefined), check 'ldglobalcheck'") 
}

#############################################################################
mIter <- optimOpt[1] #should add something to default options here if not specified
relTolX <- optimOpt[2]
reportfreq <- optimOpt[3]
detailreport <- optimOpt[4]


controlin <- list(trace=detailreport,maxit=mIter,reltol=relTolX,REPORT=reportfreq)

res <- tryCatch({

optim(starts2, fr, dat=d, otherdat=otherdat, alts=max(choice), method = methodname, control = controlin, hessian = TRUE)

}, error = function(e) {

    return("Optimization error, check 'ldglobalcheck'") 
	
})

if (res[[1]][1] == "Optimization error, check 'ldglobalcheck'") {

	return(list(errorExplain=res,OutLogit=OutLogit,optoutput=optoutput,seoutmat2=seoutmat2,MCM=MCM,H1=H1))
	
}

q2 <- res[["par"]]
LL <- res[["value"]]
output <- list(counts = res[["counts"]], convergence = res[["convergence"]], optim_message=res[["message"]])
H <- res[["hessian"]]
	
#Model comparison metrics (MCM)
param <- max(dim(as.matrix(starts2)))
obs <- dim(dataCompile)[1]
AIC <- 2*param - 2*LL 

AICc <- AIC + (2*param*(param+1))/(obs - param - 1)

BIC <- -2*LL + param*log(obs)

PseudoR2 <- (LL_start-LL)/LL_start

MCM <- list(AIC=AIC, AICc=AICc, BIC=BIC, PseudoR2=PseudoR2)

if (is.null(H) == FALSE){

H1 <- tryCatch({

solve(H)

}, error = function(e) {

    return("Error, singular, check 'ldglobalcheck'") 
	
})

diag2 <- tryCatch({

diag(H1)

}, error = function(e) {

    return("Error, NAs, check 'ldglobalcheck'") 
	
})

se2 <- tryCatch({ 

sqrt(diag2) 

}, warning = function(war) {

	print("Check 'ldglobalcheck'")
	sqrt(diag2) 
})

outmat2 <- t(q2)
seoutmat2 <- t(se2)
optoutput <- output
tLogit <- t(outmat2/se2)
OutLogit <- cbind(t(outmat2), as.matrix(se2), (tLogit))
}

#############################################################################

return(list(errorExplain=errorExplain,OutLogit=OutLogit,optoutput=optoutput,seoutmat2=seoutmat2,MCM=MCM,H1=H1))

}