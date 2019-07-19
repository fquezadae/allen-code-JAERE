logit_correction_estscale <- function(starts3, dat, otherdat, alts) {
    #' Full information model with Dahl's correction function that estimates
	#'     scale instead of cost (for personal MC use)
    #'
    #' Full information model with Dahl's correction function that estimate
	#'     scale instead of cost (for personal MC use)
    #'
    #' @param starts3 Starting values as a vector (num). For this likelihood,
	#'     the order takes: c([marginal utility from catch], [catch function
	#'     parameters], [polynomial starting parameters], [catch sigma], [scale
	#'     parameter]). \cr \cr
    #'     The number of polynomial interaction terms is currently set to 2, so
	#'     given the chosen degree 'polyn' there should be
	#'     (((polyn+1)*2) + 2)*kk polynomial starting parameters, where kk
	#'     equals the number of alternatives. The marginal utility from catch,
	#'     scale parameter, and catch sigma are of length 1 respectively. The
	#'     catch function parameters are of length (# of catch variables)*kk.
    #' @param dat Data matrix, see output from shift_sort_x, alternatives with
	#'     distance.
    #' @param otherdat Other data used in model (as list containing objects
	#'     griddat, intdat, startloc, polyn, and distance). \cr \cr
    #'     Any number of grid-specific variables are allowed, as a list of
	#'     matrices. Note for this likelihood only one interaction variable
    #'     is allowed (as we are estimating the scale parameter given a cost
	#'     parameter set to (-1)), as a list of a single matrix. Note the
	#'     variables (each as a matrix) within `griddat` and `intdat` have no
	#'     naming restrictions. Also note that `griddat` variables are dimension
	#'     *(number of observations) x (number of alternatives)*, while `intdat`
	#'     variables are dimension *(number of observations) x 1*, to be
	#'     interacted with the distance to each alternative. Grid-specific
	#'     variables may correspond to catches that vary by location, or
	#'     interaction variables may be vessel characteristics that affect how
	#'     much disutility is suffered by traveling a greater distance. Note in
	#'     this likelihood the grid-specific variables are the variables in the
	#'     catch equation, and each variable varies across observations but not
	#'     for each location: they are grid-specific due to the
	#'     location-specific coefficients. If there are no other data, the user
	#'     can set `griddat` as ones with dimension *(number of observations) x
	#'     (number of alternatives)* and `intdat` variables as ones with
	#'     dimension *(number of observations) x 1*. \cr \cr
    #'     The variable startloc is a matrix of dimension
	#'     *(number of observations) x 1*, that corresponds to the starting
	#'     location when the agent decides between alternatives. \cr \cr
    #'     The variable polyn is a vector of length 1 corresponding to the
	#'     chosen polynomial degree. \cr \cr
    #'     The variable distance is a matrix of dimension
	#'     *(number of observations) x (number of alternatives)* corresponding
	#'     to the distance to each alternative.
    #' @param alts Number of alternative choices in model as length 1 vector
	#'     (num).
    #' @return ld: negative log likelihood
    #' @section Details:
    #'     The computation of the Nth Fibonacci number in this function uses the
    #'     formula obtained from solving the recurrence relation
    #'     \deqn{F_{n} = F_{n-1} + F_{n-2},} which translates to the equation,
    #'     \eqn{x^{n} = x^{n-1} + x^{n-2}}, and then to \eqn{x^2=x+1}. The two
    #'     solutions to the latter equation  are \eqn{\phi} and \eqn{\psi},
    #'     where \deqn{\phi=\frac{1+\sqrt{5}}{2}, \qquad
	#'     \psi=\frac{1-\sqrt{5}}{2}.}
    #'     This eventually leads to the formula used in this function,
    #'     \deqn{F_{n} =\frac{\phi^{n} - \psi^{n}}{\sqrt{5}}.}
    #' @export
    #' @examples
    #' data(zi)
    #' data(catch)
    #' data(choice)
    #' data(distance)
    #' data(si)
    #' data(startloc)
    #'
    #' optimOpt <- c(1000,1.00000000000000e-08,1,0)
    #'
    #' methodname <- 'BFGS'
    #'
    #' polyn <- 3
    #' kk <- 4
    #'
    #' si2 <- sample(1:5,dim(si)[1],replace=TRUE)
    #'
    #' otherdat <- list(griddat=list(si=as.matrix(cbind(si,si,si,si)),
	#'     si2=as.matrix(cbind(si2,si2,si2,si2))),intdat=list(zi=as.matrix(zi)),
	#'     startloc=as.matrix(startloc),polyn=polyn,
	#'     distance=as.matrix(distance))
    #'
    #' initparams <- c(3, 0.5, 0.4, 0.3, 0.2, 0.55, 0.45, 0.35, 0.25,
	#'     rep(0, (((polyn+1)*2) + 2)*kk), 1, 3)
    #'
    #' func <- logit_correction_estscale
    #'
    #' results <- discretefish_subroutine(catch,choice,distance,otherdat,
	#'     initparams,optimOpt,func,methodname)
    #'
    
    griddat <- as.matrix(do.call(cbind, otherdat$griddat))
    intdat <- as.matrix(do.call(cbind, otherdat$intdat))
    
    gridnum <- dim(griddat)[2]/alts
    intnum <- dim(intdat)[2]
    # get number of variables
    
    startloc <- (otherdat$startloc)
    distance <- otherdat$distance
    
    polyn <- otherdat$polyn
    
    starts3 <- as.matrix(starts3)
    
    revcoef <- as.matrix(starts3[1:1, ])
    
    gridlength <- (gridnum * alts) + ((((polyn + 1) * 2) + 2) * alts)
    
    gridcoef <- as.matrix(starts3[2:(1 + gridlength), ])
    
    intcoef <- (-1)
    
    sigmaa <- as.matrix(starts3[((1 + gridlength) + 1), ])
    signum <- 1
    
    sigmac <- as.matrix(starts3[((1 + gridlength) + 1 + signum), ])
    # end of vector
    
    obsnum <- dim(griddat)[1]
    
    gridbetas <- (matrix(gridcoef[1:(alts * gridnum), ], obsnum, alts * gridnum,
        byrow = TRUE) * griddat)
    dim(gridbetas) <- c(nrow(gridbetas), alts, gridnum)
    gridbetas <- rowSums(gridbetas, dim = 2)
    
    intbetas <- .rowSums(intdat * matrix(intcoef, obsnum, intnum, byrow = TRUE), 
        obsnum, intnum)
    
    betas <- matrix(c((gridbetas * matrix(revcoef, obsnum, alts)), intbetas),
	    obsnum, (alts + 1))
    
    djztemp <- betas[1:obsnum, rep(1:ncol(betas), each = alts)] *
	    dat[, 3:(dim(dat)[2])]
    dim(djztemp) <- c(nrow(djztemp), ncol(djztemp)/(alts + 1), alts + 1)
    
    prof <- rowSums(djztemp, dim = 2)
    profx <- prof - prof[, 1]
    
    exb <- exp(profx/matrix(sigmac, dim(prof)[1], dim(prof)[2]))
    
    ldchoice <- (-log(rowSums(exb)))
    
    revside <- gridbetas * matrix(revcoef, obsnum, alts)
    costside <- distance * intbetas
    
    probprof <- revside + costside
    
    probprofx <- probprof - probprof[, 1]
    
    probexb <- exp(probprofx/matrix(sigmac, dim(probprof)[1], dim(probprof)[2]))
    
    probs <- probexb/matrix(rowSums(probexb), obsnum, alts)
    
    yj <- dat[, 1]
    cj <- dat[, 2]
    
    locstay <- model.matrix(~as.factor(startloc) - 1)
    locmove <- model.matrix(~as.factor(cj) - 1)
    
    probstay <- probs * locstay
    probmove <- probs * locmove
    
    intpoly <- 2
    
    movemat <- matrix(c(locmove, (matrix(probmove, obsnum, alts * polyn)^
	    matrix(rep(1:polyn, each = alts), obsnum, alts * polyn, byrow = TRUE)),
		(matrix(probmove, obsnum, alts * intpoly) * matrix(rowSums(probstay),
		obsnum, alts * intpoly))^matrix(rep(1:intpoly, each = alts), obsnum,
		alts * intpoly, byrow = TRUE)), obsnum, alts * (polyn + 1 + intpoly)) *
		matrix(!(startloc == cj), obsnum, alts * (polyn + 1 + intpoly))
	# 1 is for constant
    
    staymat <- matrix(c(locstay, (matrix(probstay, obsnum, alts * polyn)^
	    matrix(rep(1:polyn, each = alts), obsnum, alts * polyn, byrow = TRUE))),
		obsnum, alts * (polyn + 1)) * matrix((startloc == cj), obsnum,
		alts * (polyn + 1))
	# 1 is for constant
    
    Xvar <- matrix(c(griddat * matrix(locmove, obsnum, gridnum * alts), staymat, 
        movemat), obsnum, dim(gridcoef)[1])
    
    empcatches <- Xvar %*% gridcoef
    
    ldcatch <- matrix((-(0.5) * log(2 * pi)), obsnum) + (-(0.5) *
	    log(matrix(sigmaa, obsnum)^2)) + (-(0.5) * (((yj - empcatches)/
		(matrix(sigmaa, obsnum)))^2))
    
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
