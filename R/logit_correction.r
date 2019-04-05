logit_correction <- function(starts3, dat, otherdat, alts) {
    #' logit_correction
    #'
    #' Full information model with Dahl's correction function
    #'
    #' @param starts3 Starting values. e.g. c([grid-varying variables], [interaction variables], [catch variance], [scale parameter]).
    #' @param dat Data matrix, see output from shift_sort_x, alternatives with distance by column bind
    #' @param otherdat Other data used in model (as list). Any number of grid-varying variables (e.g. expected catch that varies by location) or 
    #' interaction variables (e.g. vessel characteristics that affect how much disutility is suffered by traveling a greater distance) are allowed. \cr \cr
    #' However, the user must place these in `otherdat` as list objects named `griddat` and `intdat` respectively. Note the variables #' within `griddat` 
    #' and `intdat` have no naming restrictions. \cr \cr
    #' Also note that `griddat` variables are  dimension *(number of observations) x #' (number of alternatives)*, while `intdat` variables are 
    #' dimension *(number of observations) x 1*, to be interacted with the distance to each #' alternative. \cr \cr
    #' If there are no other data, the user can set `griddat` as ones with dimension *(number of observations) x (number of alternatives)*
    #' and `intdat` variables as ones with dimension *(number of observations) x 1*.
    #' @param alts Number of alternative choices in model
    #' @return ld - negative log likelihood
    #' @export
    #' @examples
    #'
    
    ld1 <- matrix(ncol=1, nrow=dim(dat)[1])

    griddat <- (otherdat$griddat)
    intdat <- (otherdat$intdat)
    startloc <- (otherdat$startloc)
    
	polyn <- otherdat$polyn
	gridnum <- otherdat$gridnum
	intnum <- otherdat$intnum
	
    starts3 <- as.matrix(starts3)
	
	revcoef <- as.matrix(starts3[1:1, ])
	
	gridlength <- (gridnum * alts) + ((((polyn+1)*2)+2)*alts)
	
    gridcoef <- as.matrix(starts3[2:(1 + gridlength), ])
    
    # intcoef <- as.matrix(starts3[(1 + gridlength + 1):((1 + gridlength) + 
		# length(intdat)),
		# ])
    
	intcoef <- (-1)
	
    # if ((dim(starts3)[1] - (1 + gridlength + length(intdat) + 
        # 1)) == alts) {
        # sigmaa <- as.matrix(starts3[((1 + gridlength) + length(intdat) + 
            # 1):((1 + gridlength) + length(intdat) + alts), ])
        # signum <- alts
    # } else {
        # sigmaa <- as.matrix(starts3[((1 + gridlength) + length(intdat) + 
            # 1), ])
        # signum <- 1
    # }
    
    # sigmac <- as.matrix(starts3[((1 + gridlength) + length(intdat) + 
        # 1 + signum), ])  #end of vector
		
	# if ((dim(starts3)[1] - (1 + gridlength + 
        # 1)) == alts) {
        # sigmaa <- as.matrix(starts3[((1 + gridlength) + 
            # 1):((1 + gridlength) + alts), ])
        # signum <- alts
    # } else {
        # sigmaa <- as.matrix(starts3[((1 + gridlength) + 
            # 1), ])
        # signum <- 1
    # }
	
	sigmaa <- as.matrix(starts3[((1 + gridlength) + 
            1), ])
	signum <- 1
    
    sigmac <- as.matrix(starts3[((1 + gridlength) + 
        1 + signum), ])  #end of vector
		
    for (i in 1:dim(dat)[1]) {
        
        betas1 <- c(t((griddat[i,]) * t(matrix(gridcoef[1:alts, ], alts, gridnum))) %*% rep(as.matrix(
            revcoef), gridnum), t((intdat[i,])) %*% as.matrix(intcoef))
        betas <- t(as.matrix(betas1))
        
        djz <- t(dat[i, 3:dim(dat)[2]])
        
        dj <- matrix(djz, nrow = alts, ncol = dim(betas)[2])
        
        xb <- dj %*% t(betas)
        xbm <- xb - xb[1]
        exb <- exp(xbm/matrix(sigmac, dim(xbm)[1], 1))
        
        ldchoice <- (-log(t(exb) %*% (rep(1, alts))))
        
		probs <- exb/sum(exb)
				
        yj <- dat[i, 1]
        cj <- dat[i, 2]
		oriloc <- as.numeric((startloc[i,]))
		
		locstay <- rep(0, alts)
		locstay[oriloc] <- 1
		
		locmove <- rep(0, alts)
		locmove[cj] <- 1
		
		if ((alts - (cj - 1) + 1) > alts) {
		probsorder <- probs[1:(alts - (cj - 1))]
		} else {
		probsorder <- c(probs[(alts - (cj - 1) + 1):alts], probs[1:(alts - (cj - 1))])
		}
		
		probstay <- probsorder[oriloc]*locstay
		
		probmove <- probsorder[cj]*locmove
		
		probmoveint <- probsorder[oriloc]*probsorder[cj]*locmove
		
		if (oriloc == cj) {
			Xvar <- c(((griddat[i,])*locmove), locstay)
				for (pp in 1:polyn) {
					Xvar <- c(Xvar, probstay^pp)
				}
			Xvar <- c(Xvar, rep(0, alts*(polyn+1+2))) #constant plus int
		} else {
			Xvar <- c(((griddat[i,])*locmove), rep(0, alts*(polyn+1)), locmove)
				for (pp in 1:polyn) {
					Xvar <- c(Xvar, probmove^pp)
				}
			Xvar <- c(Xvar, probmoveint, probmoveint^2)
		}
						
        empcatches <- t(as.matrix(Xvar))%*%gridcoef
        
        if (signum == 1) {
            empsigmaa <- sigmaa
        } else {
            empsigmaa <- sigmaa[cj]
        }
        
        ldcatch1 <- (-(0.5) * log(2 * pi))
        ldcatch2 <- (-(0.5) * log(empsigmaa^2))
        ldcatch3 <- (-(0.5) * (((yj - empcatches)/(empsigmaa))^2))
        ldcatch <- ldcatch1 + ldcatch2 + ldcatch3
        
        ld1[i,] <- ldcatch + ldchoice
        
    }
    
	ld <- -sum(ld1)
    
    if (is.nan(ld) == TRUE) {
        ld <- .Machine$double.xmax
    }
    
    # ldsumglobalcheck <<- ld
    # paramsglobalcheck <<- starts3
    # ldglobalcheck <<- unlist(as.matrix(ld1))
    
    return(ld)
    
}
