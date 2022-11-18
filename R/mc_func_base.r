mc_func_base <- function(mcnum, constp=c(0,0,0,0), cdev, iinum, bwnum=1,
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
bikN <- list()	
for (l in 1:kk) { 

    biksigma <- c(cdev,cdev,cdev,cdev)
    bik[[l]] <- matrix(0,ii[l],kk)
    bikN[[l]] <- matrix(0,ii[l],kk)

    for (m in 1:kk) {

    bik[[l]][,m] <- matrix(rnorm(ii[l],0,biksigma[m]),ii[l],1)
    bikN[[l]][,m] <- matrix(rnorm(ii[l],0,0),ii[l],1)

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
yikT <- list()

# choose between k locations
for (k in 1:kk) {

tijk <- betac*distance[j,k]*zi[[j]] + wijk[[j]][,k]

yik[[k]] <- yconst[1,k] + betavar[k,]*si[[j]] + bik[[j]][,k]
yikT[[k]] <- yik[[k]] + bikN[[j]][,k]

Vijk[[k]] <- alpha*yik[[k]] + tijk

}

# choose max utility
choice[[j]] <- as.matrix(which(t(apply(matrix(unlist(Vijk),ii[j],kk),1,max) ==
    matrix(unlist(Vijk),ii[j],kk)))-(((1:ii[j])-1)*kk))
yikchosen[[j]] <- as.matrix(diag(matrix(unlist(yikT),ii[j],kk)[,choice[[j]]]))

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

# make probabilities based on "cells"
makeprobs <- data.frame(cbind(si=unlist(siout),zi=unlist(zifin),
    startl=unlist(startlocfin),choice=unlist(choicefin)))

makeprobs <- makeprobs %>%
    group_by(si, zi, startl) %>%
    mutate(tot = n(), staynum = sum(startl == choice)) %>%
    mutate(probstay = staynum/tot) %>%
    data.frame()

makeprobs <- makeprobs %>%
    group_by(si, zi, startl, choice) %>%
    mutate(chonum = n()) %>%
    mutate(probmove = chonum/tot) %>%
    data.frame()

makeprobs <- makeprobs %>%
    mutate(probstayz = case_when(
        (startl != choice) ~ 0,
        T ~ probstay)) %>%
    mutate(probmovez = case_when(
        (startl == choice) ~ 0,
        T ~ probmove)) %>%
    data.frame()
    
makeprobs <- makeprobs %>%
    mutate(stayv = as.numeric(startl == choice), 
        movev = as.numeric(startl !=choice)) %>%
    data.frame()

mm <- model.matrix(~as.factor(makeprobs$choice)-1)

# model matrix for corrected catch regression
siregcor <- cbind(sireg, makeprobs$stayv*mm, mm*makeprobs$probstayz, 
    mm*makeprobs$probstayz^2, mm*makeprobs$probstayz^3, mm*makeprobs$movev, 
    mm*makeprobs$probmovez, mm*makeprobs$probmovez^2, mm*makeprobs$probmovez^3, 
    mm*makeprobs$probmovez*makeprobs$probstay, 
    mm*(makeprobs$probmovez*makeprobs$probstay)^2)
colnames(siregcor) <- c("si1", "si2", "si3", "si4", "yi", "stayc1", "stayc2",
    "stayc3", "stayc4", "stay11", "stay12", "stay13", "stay14", "stay21",
    "stay22", "stay23", "stay24", "stay31", "stay32", "stay33", "stay34",
    "movec1", "movec2", "movec3", "movec4", "move11", "move12", "move13",
    "move14", "move21", "move22", "move23", "move24", "move31", "move32",
    "move33", "move34", "movei11", "movei12", "movei13", "movei14", "movei21",
    "movei22", "movei23", "movei24")

# corrected catch regression 
catchcorreg <- (lm(data=siregcor, yi~.-1))

predicted_corcatchfin <- data.frame(sifin*t(matrix(coef(catchcorreg)[1:4],
    kk,sum(ii))))
colnames(predicted_corcatchfin) <- c("V1","V2","V3","V4")

# data to pass to RUM model
intdatfin <- list(zi=zifin)
griddatfin <- list(si=sifin)
startlocdatfin <- list(startloc=startlocfin)

polyn <- 3
polyintnum <- 1
polyconstant <- 1
singlecor <- 0

otherdatfin <- list(griddat=list(as.matrix(sifin)), noCgriddat = NA,
    intdat=list(as.matrix(zifin)), startloc=as.matrix(startlocfin), 
    polyn=polyn, polyintnum = polyintnum, regconstant = regconstant, 
    polyconstant = polyconstant, singlecor = singlecor)

initparams <- unname(c(1, coef(catchreg), rep(0, 
        (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
        rep(1,dim(zifin)[2]), 1))

# Optimization options for the 
# maximum number of function evaluations, maximum iterations, and the 
# relative tolerance of x. Then, how often to report output, and whether to 
# report output.
optimOpt <- c(100000,1.00000000000000e-8,1,0)
methodname = "BFGS"

otherdatfin$distance <- as.matrix(distancefin)
otherdatfin$bw <- bw

func <- logit_correction_polyint_estscale

initcount <- 0

# naive (shotgun) search for starting parameters for minimization
searchspace <- 10000

changevec <- unname(c(1, rep(0, length(coef(catchreg))),
    rep(1, (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
    rep(1,dim(zifin)[2]), 1)) 
    
results <- explore_startparams(searchspace, initparams, dev = 2, 
    logit_correction_polyint_estscale, catchfin, choicefin, distancefin, otherdatfin,
    changevec)

LLmat <- data.frame(cbind(1:searchspace,unlist(results$saveLLstarts)))
LLmatorder <- LLmat[order(LLmat$X2),]

initparamssave <- results$savestarts[LLmatorder$X1[1:100]]
 
results_save_correction <- list()
results_save_correction$OutLogit <- matrix(NA, 2, 2)

# solve minimization, making sure convergence or maximum attempts are reached    
while (any(is.na(as.numeric(results_save_correction$OutLogit[,2]))) == TRUE &
    initcount < 12) {
    
initcount <- initcount + 1

initparams <- initparamssave[[initcount]]

results_save_correction <- discretefish_subroutine(catchfin,choicefin,
    distancefin,otherdatfin,initparams,optimOpt,func,methodname)

}

otherdatfin <- list(griddat=list(as.matrix(predicted_catchfin)),
    intdat=list(as.matrix(zifin)))

# Initial paramters for revenue then cost.
initparams <- c(3, 6) 

func <- logit_c_estscale
methodname = "BFGS"

# uncorrected model
results_save_uncorrected <- discretefish_subroutine(catchfin,choicefin,
    distancefin,otherdatfin,initparams,optimOpt,func,methodname)

otherdatfin <- list(griddat=list(as.matrix(predicted_corcatchfin)),
    intdat=list(as.matrix(zifin)))

initparams <- c(3, 6)

func <- logit_c_estscale
methodname = "BFGS"

# two stage model
results_save_2st <- discretefish_subroutine(catchfin,choicefin,
    distancefin,otherdatfin,initparams,optimOpt,func,methodname)

otherdatfin <- list(griddat=list(as.matrix(real_catchfin)),
    intdat=list(as.matrix(zifin)))

initparams <- c(3, 6) 

func <- logit_c_estscale
methodname = "BFGS"

# hypothetical "true" model
results_save_true <- discretefish_subroutine(catchfin,choicefin,distancefin,
    otherdatfin,initparams,optimOpt,func,methodname)

# test statistics from catch regression
SSEr <- sum(catchreg$residuals^2)
SSEu <- sum(catchcorreg$residuals^2)
J <- length(coef(catchcorreg))-length(coef(catchreg))
NK <- sum(ii)-length(coef(catchcorreg))

Ft <- ((SSEr-SSEu)/J)/(SSEu/NK)

# output
par_res_base <- list(
    correction=results_save_correction, 
    uncorrected=results_save_uncorrected, twostage = results_save_2st,
    true = results_save_true,
    uncorrectcoef=coef(catchreg), corcatchregcoef = catchcorreg$coef, 
    uncorcatchregvcov = vcov(catchreg),
    corcatchregvcov = vcov(catchcorreg),
    Ft = Ft,
    choicetab=table(choicefin),
    probs=makeprobs, init = initcount
    )

return(par_res_base)
    
}
