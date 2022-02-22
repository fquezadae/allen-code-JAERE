mc_func_true <- function(mcnum, constp=c(0,0,0,0), cdev, iinum, bwnum=1,
    sertype = "gev", chodev=1, alpha=3, betac = -1, regconstant = 0) {
    #' Monte Carlo function
    #'
    #' MC function for Chen et al. 2022
    #'
    #' @param mcnum MC iteration
    #' @param constp Catch regression constant values
    #' @param cdev Catch devation standard deviation
    #' @param iinum Number of observations at each location
    #' @param bwnum Bandwidth value for kernel
    #' @param sertype Choice error type
    #' @param chodev Choice error parameter
    #' @param alpha Marginal utility from catch
    #' @param betac Marginal disutility from distance
    #' @param regconstant Catch regression constant exists 
    #' @return Output from all models
    #' @export
    
# for reproducibility, pass different mcnum for each MC iteration
set.seed(mcnum)

bw <- bwnum

# number of locations
kk <- 4
ii <- iinum 

# constant for catch equation
yconst <- t((constp)*matrix(1,kk,sum(ii)))

# catch parameters
betavar <- as.matrix(c(1.50, 1.25, 1.00, 0.75))


# distance matrix
distance <- rbind(t(as.matrix(c(0.0, 0.5, 0.5, 0.707))), 
    t(as.matrix(c(0.5, 0.0, 0.707, 0.5))), 
    t(as.matrix(c(0.5, 0.707, 0, 0.5))), 
    t(as.matrix(c(0.707, 0.5, 0.5, 0))))*3

# make vessel characteristics
si <- list()
zi <- list()
for (l in 1:kk) {
    si[[l]] <- matrix(sample(1:5,ii[l]*1,replace=TRUE),ii[l],1)
    zi[[l]] <- matrix(sample(1:10,ii[l]*1,replace=TRUE),ii[l],1)
}

# make catch deviation
bik <- list()	
for (l in 1:kk) { 

    biksigma <- c(cdev,cdev,cdev,cdev)
    bik[[l]] <- matrix(0,ii[l],kk)

    for (m in 1:kk) {

    bik[[l]][,m] <- matrix(rnorm(ii[l],0,biksigma[m]),ii[l],1)

    }

}

# make choice error
wijk <- list()
for (l in 1:kk) {
    if (sertype == "gev") {
    # when chodev = 1 mu = 0 and beta = 1
    wijk[[l]] <- matrix(((1 - chodev)*(-digamma(1)))-(chodev)*
        log(-log(runif(4000, min = 0, max = 1))), ii[l],kk) 
    } else if (sertype == "norm") {
    wijk[[l]] <- matrix(rnorm(ii[l]*kk,0,chodev),ii[l],kk)
    } else {
    wijk[[l]] <- matrix((-log(rexp(ii[l]*kk,1))),ii[l],kk) 
    # -log(exp(1)) is standard type 1 extreme value i.e. gumbel beta=1 mu=0
    }
}

# generate catches, have vessels make choices, and track them
choice <- list()
yikchosen <- list()
siout <- list()
siout2 <- list()
ziout <- list()
ziout2 <- list()
startlocout <- list()
distanceout <- list()
predyik <- list()
predyik2 <- list()
# for vessels starting at location j
for (j in 1:kk) {

Vijk <- list()
yik <- list()

# choose between k locations
for (k in 1:kk) {

tijk <- betac*distance[j,k]*zi[[j]] + wijk[[j]][,k]

yik[[k]] <- yconst[1,k] + betavar[k,]*si[[j]] + bik[[j]][,k]


Vijk[[k]] <- alpha*yik[[k]] + tijk

}

# choose max utility
choice[[j]] <- as.matrix(which(t(apply(matrix(unlist(Vijk),ii[j],kk),1,max) ==
    matrix(unlist(Vijk),ii[j],kk)))-(((1:ii[j])-1)*kk))
yikchosen[[j]] <- as.matrix(diag(matrix(unlist(yik),ii[j],kk)[,choice[[j]]]))

# track vessel characteristics
siout[[j]] <- si[[j]]
ziout[[j]] <- zi[[j]]
startlocout[[j]] <- rep(j,ii[j])
distanceout[[j]] <- t(matrix(rep(distance[j,],ii[j]),kk,ii[j]))

}

# unlist tracked variables for regression
dummymatrix <- as.matrix(model.matrix(~as.factor(unlist(choice))-1))

zifin <- data.frame(V1 = as.numeric(unlist(ziout)))
startlocfin <- data.frame(V1 = as.numeric(unlist(startlocout)))
sifin <- data.frame(V1 = as.numeric(unlist(siout)),V2 = 
    as.numeric(unlist(siout)),V3 = as.numeric(unlist(siout)),V4 = 
    as.numeric(unlist(siout)))
sifin <- matrix(as.numeric(unlist(siout)),sum(ii),kk)

sireg <- data.frame(sifin*cbind(dummymatrix))
colnames(sireg) <- c("V1","V2","V3","V4")

sireg$yireg <- unlist(yikchosen)

# uncorrected catch regression
catchreg <- (lm(data=sireg, yireg~V1+V2+V3+V4-1))

# make predicted catch data
catchfin <- data.frame(V1 = unlist(yikchosen))
choicefin <- data.frame(V1 = unlist(choice))
predicted_catchfin <- data.frame(sifin*t(matrix(coef(catchreg),kk,sum(ii))))
colnames(predicted_catchfin) <- c("V1","V2","V3","V4")

distancefin <- data.frame(do.call(rbind,distanceout))
colnames(distancefin) <- c("V1","V2","V3","V4")

real_catchfin <- data.frame(sifin[,1:4]*t(matrix(betavar,kk,sum(ii))))
colnames(real_catchfin) <- c("V1","V2","V3","V4")

# data to pass to RUM model
intdatfin <- list(zi=zifin)
griddatfin <- list(si=sifin)
startlocdatfin <- list(startloc=startlocfin)

# Optimization options for the 
# maximum number of function evaluations, maximum iterations, and the 
# relative tolerance of x. Then, how often to report output, and whether to 
# report output.
optimOpt <- c(100000,1.00000000000000e-8,1,0)
methodname = "BFGS"

otherdatfin <- list(griddat=list(as.matrix(real_catchfin)),
    intdat=list(as.matrix(zifin)))

initparams <- c(3, 6) 

func <- logit_c_estscale
methodname = "BFGS"

# hypothetical "true" model
results_save_true <- discretefish_subroutine(catchfin,choicefin,distancefin,
    otherdatfin,initparams,optimOpt,func,methodname)

# output
par_res_base <- list(
    true = results_save_true,
    uncorrectcoef=coef(catchreg),
    choicetab=table(choicefin)
    )

return(par_res_base)
    
}
