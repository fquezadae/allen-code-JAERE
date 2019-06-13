logit_avgcat_v <- function(starts3, dat, otherdat, alts) {
    #' logit_avgcat_v
    #'
    #' Average catch conditional logit procedure
    #'
    #' @param starts3 Starting values with dimensions equal to (number of alternatives - 1)*(number of grid-varying variables) + (number of interactions).
    #' Recall the user must drop one alternative for identification. e.g. c([grid-varying variables -1], [interaction variables]).
    #' @param dat Data matrix, see output from shift_sort_x, alternatives with distance by column bind
    #' @param otherdat Other data used in model (as list). Any number of grid-varying variables (e.g. variables that affect catches across locations) or 
    #' interaction variables (e.g. vessel characteristics that affect how much disutility is suffered by traveling a greater distance) are allowed. \cr \cr
    #' However, the user must place these in `otherdat` as list objects named `griddat` and `intdat` respectively. Note the variables
    #' within `griddat` and `intdat` have no naming restrictions. \cr \cr
    #' Also note that `griddat` variables are dimension *(number of observations) x 1*, #' to be interacted with each alternative, 
    #' while `intdat` variables are dimension *(number of observations) x 1*, to be interacted with the distance to each alternative. \cr \cr
    #' If there are no other data, the user can set `griddat` as ones with dimension *(number of observations) x 
    #' (number of alternatives)* and `intdat` variables as ones with dimension *(number of observations) x 1*.
    #' @param alts Number of alternative choices in model
    #' @param project Name of project
    #' @param expname Expected catch table
    #' @param mod.name Name of model run for model result output table
    #' @return ld - negative log likelihood
    #' @export
    #' @examples
    #'
	
	griddat <- as.matrix(do.call(cbind, otherdat$griddat))
    intdat <- as.matrix(do.call(cbind, otherdat$intdat))
	
	gridnum <- dim(griddat)[2]
	intnum <- dim(intdat)[2]
	#get number of variables
	
	obsnum <- dim(griddat)[1]

    starts3 <- as.matrix(starts3)
    gridcoef <- as.matrix(starts3[1:(gridnum * (alts - 1)), ])
    intcoef <- as.matrix(starts3[((gridnum * (alts - 1)) + 1):(((gridnum * 
        (alts - 1))) + intnum), ])
    
   	#############################################

	gridbetas <- (matrix(gridcoef,obsnum,(alts-1)*gridnum,byrow=TRUE)*griddat[, rep(1:gridnum, each=(alts-1))])
	dim(gridbetas) <- c(nrow(gridbetas), (alts-1), gridnum)
	gridbetas <- rowSums(gridbetas,dim=2)
	
	intbetas <- .rowSums(intdat*matrix(intcoef,obsnum,intnum,byrow=TRUE),obsnum,intnum)
	
	betas <- matrix(c(gridbetas, intbetas),obsnum,(alts-1+1))

	djztemp <- betas[1:obsnum,rep(1:ncol(betas), each = (alts))]*dat[, (alts+3):(dim(dat)[2])]
	dim(djztemp) <- c(nrow(djztemp), ncol(djztemp)/((alts-1)+1), (alts-1)+1)

	prof <- rowSums(djztemp,dim=2)
	profx <- prof - prof[,1]

	exb <- exp(profx)

	ldchoice <- (-log(rowSums(exb)))

	#############################################
	
	ld <- -sum(ldchoice)
	
    if (is.nan(ld) == TRUE) {
        ld <- .Machine$double.xmax
    }
    
    ldsumglobalcheck <- ld
    assign('ldsumglobalcheck', value = ldsumglobalcheck, pos = 1)
    paramsglobalcheck <- starts3
    assign('paramsglobalcheck', value = paramsglobalcheck, pos = 1)
    ldglobalcheck <- unlist(as.matrix(ldchoice))
    assign('ldglobalcheck', value = ldglobalcheck, pos = 1)
	
    return(ld)
    
}
