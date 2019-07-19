epm_weibull <- function(starts3, dat, otherdat, alts) {
    #' Expected profit model weibull catch function
    #'
    #' Expected profit model weibull catch function
    #'
    #' @param starts3 Starting values as a vector (num). For this likelihood,
    #'     the order takes: c([catch function parameters], [cost (distance)
    #'     parameters], [catch sigma(s)], [scale parameter]). \cr \cr
	#'     The catch function and cost parameters are of length (# of catch
    #'     variables)*kk and (# of cost variables) respectively, where kk equals
    #'     the number of alternatives. The catch sigma(s) are either of length 1
    #'     or length kk (if the analyst is estimating a location-specific catch
	#'     parameter). The scale parameter is of length 1.
    #' @param dat Data matrix, see output from shift_sort_x, alternatives with
	#'     distance.
    #' @param otherdat Other data used in model (as list containing objects
	#'     griddat, intdat, and prices). \cr \cr
    #'     For grid-specific variables `griddat` and cost variables to be
	#'     interacted with distance `intdat`, any number of variables are
    #'     allowed, as a list of matrices. Note the variables (each as a matrix)
    #'     within `griddat` and `intdat` have no naming restrictions. Also note
    #'     that `griddat` variables are dimension *(number of observations) x
	#'     (number of alternatives)*, while `intdat` variables are dimension
    #'     *(number of observations) x 1*, to be interacted with the distance to
    #'     each alternative. Grid-specific variables may correspond to catches
	#'     that vary by location, or interaction variables may be vessel
	#'     characteristics that affect how much disutility is suffered by
    #'     traveling a greater distance. Note in this likelihood the
	#'     grid-specific variables are the variables in the catch equation, and
    #'     each variable varies across observations but not for each location:
    #'     they are grid-specific due to the location-specific coefficients. If
	#'     there are no other data, the user can set `griddat` as ones with
    #'     dimension *(number of observations) x (number of alternatives)* and
	#'     `intdat` variables as ones with dimension *(number of observations)
	#'     x 1*. \cr \cr
    #'     The variable prices is a matrix of dimension *(number of
	#'     observations) x 1*, corresponding to prices.
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
    #' data(prices)
    #'
    #' catch[catch<0] <- 0.00001
    #' #Note weibull catch distribution.
    #'
    #' optimOpt <- c(1000,1.00000000000000e-08,1,0)
    #'
    #' methodname <- 'BFGS'
    #'
    #' si2 <- sample(1:5,dim(si)[1],replace=TRUE)
    #' zi2 <- sample(1:10,dim(zi)[1],replace=TRUE)
    #'
    #' otherdat <- list(griddat=list(si=as.matrix(cbind(si,si,si,si)),
	#'     si2=as.matrix(cbind(si2,si2,si2,si2))),
    #'     intdat=list(zi=as.matrix(zi),zi2=as.matrix(zi2)),
	#'     pricedat=list(prices=as.matrix(prices)))
    #'
    #' initparams <- c(2.5, 2.0, 1.5, 1.0, 1.1, 1.05, 0.9, 0.8, -0.8, -0.4, 3,
	#'     2, 3.5, 2.5, 1)
    #'
    #' func <- epm_weibull
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
    
    pricedat <- as.matrix(unlist(otherdat$pricedat))
    
    starts3 <- as.matrix(starts3)
    gridcoef <- as.matrix(starts3[1:(gridnum * alts), ])
    
    intcoef <- as.matrix(starts3[(((gridnum * alts) + intnum) - intnum + 1):
	    ((gridnum * alts) + intnum), ])
    
    if ((dim(starts3)[1] - ((gridnum * alts) + intnum + 1)) == alts) {
	
        k <- as.matrix(starts3[((gridnum * alts) + intnum + 1):
		    ((gridnum * alts) + intnum + alts), ])
        signum <- alts
	
    } else {
	
        k <- as.matrix(starts3[((gridnum * alts) + intnum + 1), ])
        signum <- 1
	
    }
    # if number of parameters before scale parameter is not equal to number of
	    # alts, use first parameter as sigma catch
    
    k <- sqrt(k^2)
    
    sigmac <- as.matrix(starts3[((gridnum * alts) + intnum + 1 + signum), ])
    # end of vector
    
    gridbetas <- (matrix(gridcoef, obsnum, alts * gridnum, byrow = TRUE) *
	    griddat)
    dim(gridbetas) <- c(nrow(gridbetas), alts, gridnum)
    gridbetas <- rowSums(gridbetas, dim = 2)
    
    gridmu <- sqrt(gridbetas^2)
    
    expgridcoef <- gridmu * matrix(gamma(1 + (1/k)), obsnum, alts)
    
    intbetas <- .rowSums(intdat * matrix(intcoef, obsnum, intnum, byrow = TRUE),
        obsnum, intnum)
    
    betas <- matrix(c((expgridcoef * matrix(pricedat, obsnum, alts)), intbetas), 
        obsnum, (alts + 1))
    
    djztemp <- betas[1:obsnum, rep(1:ncol(betas), each = alts)] *
	    dat[, 3:(dim(dat)[2])]
    dim(djztemp) <- c(nrow(djztemp), ncol(djztemp)/(alts + 1), alts + 1)
    
    prof <- rowSums(djztemp, dim = 2)
    profx <- prof - prof[, 1]
    
    exb <- exp(profx/matrix(sigmac, dim(prof)[1], dim(prof)[2]))
    
    ldchoice <- (-log(rowSums(exb)))
    
    yj <- dat[, 1]
    cj <- dat[, 2]
    
    if (signum == 1) {
        empk <- k
    } else {
        empk <- k[cj]
    }
    
    empgridbetas <- t(gridcoef)
    dim(empgridbetas) <- c(nrow(empgridbetas), alts, gridnum)
    
    empgriddat <- griddat
    dim(empgriddat) <- c(nrow(empgriddat), alts, gridnum)
    
    empgridmu <- .rowSums(empgridbetas[, cj, ] * empgriddat[, 1, ], obsnum,
	    gridnum)
    # note grid data same across all alternatives
    
    empgridmu <- sqrt(empgridmu^2)
    
    ldcatch <- (matrix((log(empk)), obsnum)) + (matrix((-(empk)), obsnum) *
	    log(empgridmu)) + (matrix((empk - 1), obsnum) * log(yj)) +
		(-((yj/empgridmu)^(matrix(empk, obsnum))))
    
    ld1 <- ldcatch + ldchoice
    
    ld <- -sum(ld1)
    
    if (is.nan(ld) == TRUE) {
        ld <- .Machine$double.xmax
    }
    
    ldsumglobalcheck <- ld
    assign("ldsumglobalcheck", value = ldsumglobalcheck, pos = 1)
    paramsglobalcheck <- starts3
    assign("paramsglobalcheck", value = paramsglobalcheck, pos = 1)
    ldglobalcheck <- unlist(as.matrix(ld1))
    assign("ldglobalcheck", value = ldglobalcheck, pos = 1)
    
    return(ld)
    
}
