logit_c <- function(starts3, dat, otherdat, alts) {
    #' Conditional logit likelihood
    #'
    #' Conditional logit likelihood
    #'
    #' @param starts3 Starting values as a vector (num). For this likelihood,
    #'     the order takes: c([grid-specific parameters], [cost (distance)
    #'     parameters]). \cr \cr
	#'     The grid-specific parameters function and cost parameters are of
	#'     length (# of grid-specific variables) and (# of cost parameters)
	#'     respectively.
    #' @param dat Data matrix, see output from shift_sort_x, alternatives with
	#'     distance.
    #' @param otherdat Other data used in model (as list containing objects
	#'     griddat and intdat). \cr \cr
    #'     For grid-specific variables `griddat` and cost variables to be
	#'     interacted with distance `intdat`, any number of variables are
	#'     allowed, as a list of matrices. Note the variables (each as a matrix)
	#'     within `griddat` and `intdat` have no naming restrictions.
    #'     Grid-specific variables may correspond to catches that vary by
    #'     location, or interaction variables may be vessel characteristics that
	#'     affect how much disutility is suffered by traveling a greater
    #'     distance. Note in this likliehood the grid-specific variables vary
    #'     across alternatives, where each variable may have been estimated in a 
    #'     previous procedure. Also note that `griddat` variables are dimension
    #'     *(number of observations) x (number of alternatives)*, while `intdat`
    #'     variables are dimension *(number of observations) x 1*, to be
    #'     interacted with the distance to each alternative. If there are no
	#'     other data, the user can set `griddat` as ones with dimension
	#'     *(number of observations) x (number of alternatives)* and `intdat`
	#'     variables as ones with dimension *(number of observations) x 1*.
    #' @param alts Number of alternative choices in model as length 1 vector
	#'     (num).
    #' @return ld: negative log likelihood
    #' @export
    #' @examples
    #' data(zi)
    #' data(catch)
    #' data(choice)
    #' data(distance)
    #' data(si)
    #'
    #' optimOpt <- c(1000,1.00000000000000e-08,1,0)
    #'
    #' methodname <- 'BFGS'
    #'
    #' kk <- 4
    #'
    #' si2 <- matrix(sample(1:5,dim(si)[1]*kk,replace=TRUE),dim(si)[1],kk)
    #' zi2 <- sample(1:10,dim(zi)[1],replace=TRUE)
    #'
    #' otherdat <- list(griddat=list(predicted_catch=as.matrix(predicted_catch),
	#'     si2=as.matrix(si2)), intdat=list(zi=as.matrix(zi),
	#'     zi2=as.matrix(zi2)))
    #'
    #' initparams <- c(2.5, 2, -1, -2)
    #'
    #' func <- logit_c
    #'
    #' results <- discretefish_subroutine(catch,choice,distance,otherdat,
	#'     initparams,optimOpt,func,methodname)
    #'
    
    griddat <- as.matrix(do.call(cbind, otherdat$griddat))
    intdat <- as.matrix(do.call(cbind, otherdat$intdat))
    
    gridnum <- dim(griddat)[2]/alts
    intnum <- dim(intdat)[2]
    # get number of variables
    
    obsnum <- dim(griddat)[1]
    
    starts3 <- as.matrix(starts3)
    gridcoef <- as.matrix(starts3[1:gridnum, ])
    intcoef <- as.matrix(starts3[((gridnum + intnum) - intnum + 1):
	    (gridnum + intnum), ])
    # split parameters for grid and interactions
    
    gridbetas <- (matrix(rep(gridcoef, each = alts), obsnum, alts * gridnum,
	    byrow = TRUE) * griddat)
    dim(gridbetas) <- c(nrow(gridbetas), alts, gridnum)
    gridbetas <- rowSums(gridbetas, dim = 2)
    
    intbetas <- .rowSums(intdat * matrix(intcoef, obsnum, intnum, byrow = TRUE),
        obsnum, intnum)
    
    betas <- matrix(c(gridbetas, intbetas), obsnum, (alts + 1))
    
    djztemp <- betas[1:obsnum, rep(1:ncol(betas), each = alts)] *
	    dat[, 3:(dim(dat)[2])]
    dim(djztemp) <- c(nrow(djztemp), ncol(djztemp)/(alts + 1), alts + 1)
    
    prof <- rowSums(djztemp, dim = 2)
    profx <- prof - prof[, 1]
    
    exb <- exp(profx)
    
    ldchoice <- (-log(rowSums(exb)))
    
    ld <- -sum(ldchoice)
    
    if (is.nan(ld) == TRUE) {
        ld <- .Machine$double.xmax
    }
    
    ldsumglobalcheck <- ld
    assign("ldsumglobalcheck", value = ldsumglobalcheck, pos = 1)
    paramsglobalcheck <- starts3
    assign("paramsglobalcheck", value = paramsglobalcheck, pos = 1)
    ldglobalcheck <- unlist(as.matrix(ldchoice))
    assign("ldglobalcheck", value = ldglobalcheck, pos = 1)
    
    return(ld)
    
}
