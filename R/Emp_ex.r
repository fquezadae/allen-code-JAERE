library(plyr)
library(dplyr)
options(dplyr.width = Inf)
library(maps)
library(ggplot2)
library(barebones.FishSET)
library(doParallel)
library(foreach)
options(scipen = 99)

finaldata <- read.csv(paste0(
    #CV haul data here, confidential data is not included))

# first hauls don't have a previous location (not using ports)
finaldata <- subset(finaldata, is.na(finaldata$lag.adfgstat6) == FALSE)

# looked up two vessels by hand  
finaldata$Grosstons[finaldata$ADFGnumber == 56153] <- 1318  
finaldata$Grosstons[finaldata$ADFGnumber == 56154] <- 1318  

# year 2015
ismall <- 2014
ibig <- 2016
obsnum <- 20
finaldata <- subset(finaldata, finaldata$Year > ismall)
finaldata <- subset(finaldata, finaldata$Year < ibig)

# B season
finaldata <- subset(finaldata, finaldata$Week >= 23)
finaldata <- subset(finaldata, finaldata$Week <= 44)

outname <- paste("CVhaul",ismall,ibig,obsnum,sep="_")
dirout <- paste0(# your directory here)

# get spatial data to make sure it's BSAI
adfgmatch <- read.csv(paste0(dirout, "Groundfish_Statistical_Areas_2001.csv"))
adfgmatch <- adfgmatch[, c("STAT_AREA", "FMP_AREA_CODE")]
adfgmatch$ADFGstat6 <- adfgmatch$STAT_AREA
adfgmatch$STAT_AREA <- NULL

#585430 probably a combination of 585431 and 585432, both which don't exist
#in your polygons but do exist on the map (as a square split in two)
adfgmatch <- rbind(adfgmatch, c("GOA", 585430))

#655401 probably a typo from your polygons, there's no 655410 in that file,
#but 655410 exists on the map
adfgmatch <- rbind(adfgmatch, c("BSAI", 655401))

adfgmatch$ADFGstat6 <- as.numeric(adfgmatch$ADFGstat6)
finaldata <- merge(finaldata, adfgmatch, by = c("ADFGstat6"), all.x = TRUE, 
    all.y=FALSE)
finaldata <- finaldata[finaldata$FMP_AREA_CODE == "BSAI", ]

###############################################################################
###############################################################################

subsetttold <- c(0)
subsettt <- sort(as.numeric(unique(c(levels(as.factor(finaldata$lag.adfgstat6)),
    levels(as.factor(finaldata$ADFGstat6))))))

# iteratively make sure each location has at least three unique vessels that
# moved there (for confidentiality) and 20 observations (need number of
# covariates in catch equation plus number of polynomial terms including 
# intercepts and interaction terms for just identified)
while (length(subsetttold) != length(subsettt) || 
    any(subsetttold != subsettt)) {

subsetttold <- subsettt

finaldata$movedum <- finaldata$ADFGstat6 != finaldata$lag.adfgstat6

uniquecount <- finaldata %>%
    group_by(ADFGstat6,movedum) %>%
    summarise(n_distinct(ADFGnumber)) #number of different vessels
    
uniquecount <- data.frame(uniquecount)
names(uniquecount)[3] <- "countend"

uniquecount2 <- uniquecount$ADFGstat6[(uniquecount$countend >= 3)]

subsettt1 <- uniquecount2[duplicated(uniquecount2)]

###############################################################################

uniquecountobs <- finaldata %>%
    group_by(ADFGstat6) %>%
    tally() #number of obs

names(uniquecountobs)[2] <- "countend"

uniquecountobs2 <- uniquecountobs[(uniquecountobs$countend >= obsnum),] %>% 
    group_by(ADFGstat6) %>%
    filter(n() > 0)

uniquecountobs3 <- uniquecountobs2$ADFGstat6

subsettt2 <- intersect(uniquecountobs3,subsettt1)

subsettt <- subsettt2
    
finaldata <- finaldata[(finaldata$ADFGstat6 %in% subsettt) & 
    (finaldata$lag.adfgstat6 %in% subsettt), ]

subsettt <- as.numeric(subsettt)
    
}
 
finaldata <- finaldata[with(finaldata, order(Id)), ]
finaldata2 <- finaldata
finaldata2$Portname <- droplevels(finaldata2$Portname)

# get data for mapping
adfgstat6 = read.table(paste0(dirout, "pollockcvhaul_adfgstat6_map.csv"), 
    header=TRUE, 
    sep=",")
names(adfgstat6)[2] <- "ADFGstat6"
names(adfgstat6)[3] <- "centroidlonin"
names(adfgstat6)[4] <- "centroidlatin"
adfgstat6 <- adfgstat6[c("ADFGstat6", "centroidlonin", "centroidlatin")]
adfgstat6 <- adfgstat6[!duplicated(adfgstat6$ADFGstat6), ]
adfgin2 <- adfgstat6[adfgstat6$ADFGstat6 %in% subsettt, ]

deg2rad <- function(deg) return(deg*pi/180)
# Haversine formula (hf) Mario Pineda-Krch
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(pmin(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

# get centroids in radians for distance
adfgin2$centroidloninrad <- deg2rad(adfgin2$centroidlonin)
adfgin2$centroidlatinrad <- deg2rad(adfgin2$centroidlatin)
adfgin2 <- adfgin2[with(adfgin2, order(ADFGstat6)), ]
finaldata2 <- finaldata2[with(finaldata2, order(ADFGstat6)), ]
adfgin2$locjj <- as.factor(as.numeric(as.factor(adfgin2$ADFGstat6)))

# normalize vessel characteristics
finaldata2$Horsepower <- finaldata2$Horsepower/mean(finaldata2$Horsepower)
finaldata2$Length <- finaldata2$Length/mean(finaldata2$Length)
finaldata2$Age <- finaldata2$Age/mean(finaldata2$Age)
finaldata2$Grosstons <- finaldata2$Grosstons/mean(finaldata2$Grosstons)

# make distance matrix
M <- finaldata2[,c("Horsepower","lag.lonrad","lag.latrad")]
distanceRnoHP <- t(apply(M, 1, function(x) gcd.hf(adfgin2$centroidloninrad,
    adfgin2$centroidlatinrad, x["lag.lonrad"], x["lag.latrad"])))

# make factors for starting and chosen locations
finaldata2$startlocjj <- as.factor(as.numeric(
    as.factor(finaldata2$lag.adfgstat6)))
finaldata2$chosenlocjj <- as.factor(as.numeric(as.factor(finaldata2$ADFGstat6)))

dat <- finaldata2

# now run model
scaledata <- 100
distancescale <- 100
regconstant <- 0
weightscale <- 1/100

#make data matrices. no scaling here on vessel characteristics (besides
#normalization above). (finding right scaling can help convergence). scaling on
#distance and weight as noted in the paper
dataout = finaldata2[,c("chosenlocjj", "Weight", "Age", "Horsepower",
    "Grosstons", "Duration")] #weight is in kg
dataout$Weight <- (dataout$Weight/1000*weightscale) #weight now in MT/100
dataout$Age <- (dataout$Age/scaledata*100)
dataout$Horsepower <- (dataout$Horsepower/scaledata*100)
dataout$Grosstons <- (dataout$Grosstons/scaledata*100)
finaldata2[,c("Length")] <- (finaldata2[,c("Length")]/scaledata*100)
finaldata2[,c("Grosstons")] <- (finaldata2[,c("Grosstons")]/scaledata*100)
catchcovariate <- as.numeric(dataout$Age)*
    as.numeric(dataout$Horsepower)
Ycol <- finaldata2[,c("Length","startlocjj","chosenlocjj","Weight")]
Ycol$Weight <- (Ycol$Weight/1000*weightscale) #weight now in MT/100
Ycol$Length <- catchcovariate

#uncorrected catch estimates 
Xmat = as.matrix(
    model.matrix(~ chosenlocjj - 1, data=Ycol)*catchcovariate)
betaout = coef(lm(Ycol$Weight~Xmat-1))

#sample counts
count = ddply(finaldata2, .(chosenlocjj), summarise, rating.mean=length(Week))
count$perc = count$rating.mean/dim(finaldata2)[1]

ii <- as.matrix(summary(dataout$chosenlocjj))
kk <- max(as.numeric(dataout$chosenlocjj))

#make variable matrix for costs
zifin <- data.frame(V1 = rep(1,length(dataout$Age)),
    V2 = as.numeric(dataout$Grosstons),
    V3 = as.numeric(dataout$Horsepower),
    V4 = as.numeric(dataout$Age))

startlocfin <- data.frame(V1 = as.numeric(Ycol$startlocjj))

#vessel catch covariate
sifin <- data.frame(cbind(
    matrix(catchcovariate,sum(ii),kk)))
colnames(sifin) <- sub("X", "V", colnames(sifin))

# make data for predicted catches
betaoutmat <- matrix(betaout,kk,1)
sifinpred <- data.frame(matrix(catchcovariate,sum(ii),kk))

#predicted catches
predicted_catchfin <- data.frame(sifinpred[,1:kk]*
    t(matrix((betaoutmat[,1]),kk,sum(ii))))

catchfin <- data.frame(V1 = dataout$Weight)
choicefin <- data.frame(V1 = as.numeric(dataout$chosenlocjj))

# data for distance 
distancefin <- data.frame(distanceRnoHP/distancescale)
colnames(distancefin) <- sub("X", "V", colnames(distancefin))

intdatfin <- list(zi=zifin)
griddatfin <- list(si=sifin)
startlocdatfin <- list(startloc=startlocfin)

# parameters for model, including number of polynomial terms, bandwidth, etc.
polyintnum <- 0
polyconstant <- 1
polyn <- 2
bw <- 0.90
singlecor <- 1
reltolin <- 9
reltol <- 1*10^(-reltolin)

optimOpt <- c(10000,reltol,1,0) #Optimization options for the 
    #maximum number of function evaluations, maximum iterations, and the 
    #relative tolerance of x. Then, how often to report output, and whether to 
    #report output.
methodname <- "BFGS"

# first uncorrected discrete choice estimates    
otherdatfin <- list(griddat=list(predcatch = as.matrix(predicted_catchfin)),
    intdat=list(as.matrix(zifin)))
initparams <- c(1, rep(-1,dim(zifin)[2]))

func <- logit_c
methodname = "BFGS"

# print(Sys.time())
results_save_uncorrected <- discretefish_subroutine(catchfin,choicefin,
    distancefin,otherdatfin,initparams,optimOpt,func,methodname)
# print(Sys.time())

###############################################################################

# then corrected discrete choice estimates    
otherdatfin <- list(griddat = list(as.matrix(sifin)),
    noCgriddat = NA,
    intdat = list(as.matrix(zifin)), startloc = as.matrix(startlocfin), 
    polyn = polyn, polyintnum = polyintnum, regconstant = regconstant, 
    polyconstant = polyconstant, singlecor = singlecor)

initparams <- unname(c(1, betaout, rep(0, 
    (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
    rep(-1,dim(zifin)[2]), 1))

otherdatfin$bw <- bw

func <- logit_correction_polyint
otherdatfin$distance <- as.matrix(distancefin)

###############################################################################

# model doesn't always converge, make sure hessian is invertible
initcount <- 0
    
results_save_correction <- discretefish_subroutine(catchfin,choicefin,
    distancefin,otherdatfin,initparams,optimOpt,func,methodname)

searchspace <- 10000

# this is a "shotgun" vector of initial starting parameters to try
changevec <- unname(c(1, rep(0, length(betaoutmat[,1])),
    rep(1, (((polyn+polyconstant)*(1+(1-singlecor))) + polyintnum)*kk), 
    rep(1,dim(zifin)[2]), 1)) 
results <- explore_startparams(searchspace, initparams, dev = 2, 
    logit_correction_polyint, catchfin, choicefin, distancefin, otherdatfin,
    changevec)

LLmat <- data.frame(cbind(1:searchspace,unlist(results$saveLLstarts)))
LLmatorder <- LLmat[order(LLmat$X2),]

initparamssave <- results$savestarts[LLmatorder$X1[1:100]]

#try up to 20 times as long as hessian isn't invertible
while ((any(is.na(as.numeric(results_save_correction$OutLogit[,2]))) == TRUE) & 
    initcount < 20) {

initcount <- initcount + 1
    
initparams <- initparamssave[[initcount]]

results_save_correction <- discretefish_subroutine(catchfin,choicefin,
    distancefin,otherdatfin,initparams,optimOpt,func,methodname)

}

# get likelihood and other model estimates
alts <- dim(otherdatfin$distance)[2]
ab <- max(choicefin) + 1
    # no interactions in create_logit_input - interact distances in likelihood
        # function instead
dataCompile <- create_logit_input(choicefin)
d <- shift_sort_x(dataCompile, choicefin, catchfin, distancefin, 
    max(choicefin), 
    ab)
LLout <- logit_correction_polyint_res(
    results_save_correction$OutLogit[,1], dat=d, otherdat=otherdatfin, alts)

# try nested model with no correction
initparamscorr <- initparams
initparams <- c(initparamscorr[1:(1+kk)], 
initparamscorr[(length(initparamscorr) -
    dim(zifin)[2]):length(initparamscorr)])

func <- logit_nocorrection

results_save_nocorrection <- discretefish_subroutine(catchfin,choicefin,
    distancefin,otherdatfin,initparams,optimOpt,func,methodname)

#test nested model
pchisq(2*(results_save_correction$MCM$LL-results_save_nocorrection$MCM$LL), 
    df = length(results_save_correction$OutLogit[,1]) - 
    length(results_save_nocorrection$OutLogit[,1]), lower.tail = FALSE)

outdat <- (list(correction=results_save_correction,
    nocorrection=results_save_nocorrection,
    uncorrected=results_save_uncorrected, uncorrectcoef=betaout, 
    choicetab=table(choicefin), probs = probmovesave #note probmovesave global
    initparams = initparamscorr, LLout = LLout))
    
saveRDS(outdat, paste0(dirout, outname, "_",
    initcount, ".rds"))

results <- outdat

# make figures and results now
##############################################################################
FIML <- as.data.frame(results$correction$OutLogit)
uncorrect <- as.data.frame(results$uncorrected$OutLogit)
colnames(FIML) <- c("coef","se","tstat")
colnames(uncorrect) <- c("coef","se","tstat")

uncorrectLL = as.data.frame((2*dim(uncorrect)[1] - 
    results$uncorrected$MCM$AIC)/2)
colnames(uncorrectLL) <- c("LL")

# print model statistics
print("Uncorrect LL")
print(uncorrectLL)
print("FIML LL")
print(LLout)
print("No correct MCM")
print(results$nocorrection$MCM)
print("FIML MCM")
print(results$correction$MCM)

# get catch coefficients
uncorrectcatch = as.data.frame(results$uncorrectcoef)
colnames(uncorrectcatch) <- c("coef")
FIMLcatch = as.data.frame(FIML[2:(1+dim(uncorrectcatch)[1]),1])
colnames(FIMLcatch) <- c("coef")

# get catch covariate data
matFIML <- sifin
coefFIML <- (as.matrix(FIMLcatch$coef))

# get full information revenue matrix
FIMLrev <- FIML[1,1]*((t(matrix(rep(t(coefFIML),dim(matFIML)[1]),
    ncol=dim(matFIML)[1])))*matFIML)
covarnum <- dim(FIMLrev)[2]/alts

print("FIML catch")
rowSums(matrix(as.matrix(FIMLcatch),alts,covarnum))
print("uncorrected catch")
rowSums(matrix(as.matrix(uncorrectcatch),alts,covarnum))

# make data for results (just pulling from above)
matuncorrect <- sifin
coefuncorrect <- (as.matrix(uncorrectcatch$coef))
uncorrectrev <- uncorrect[1,1]*((t(matrix(rep(t(coefuncorrect),
    dim(matuncorrect)[1]),ncol=dim(matuncorrect)[1])))*matuncorrect)
costmat <- zifin
distanceR <- distancefin 
cost <- 4
numparams <- dim(FIML)[1]
noClen <- 0

FIMLcost <- matrix(rowSums(t(matrix(FIML[(numparams-cost):(numparams-1),1], 
    dim(costmat)[2], dim(costmat)[1]))*(costmat)), dim(distanceR)[1], alts)*
    distanceR
uncorrectcost <- matrix(rowSums(t(matrix(uncorrect[
    ((noClen + 2): 
    ((noClen + 2) + cost-1)),1], dim(costmat)[2], 
    dim(costmat)[1]))*(costmat)), dim(distanceR)[1], alts)*distanceR

# make map of data
adfgstat6 <- read.table(paste0(dirout,
    "adfgstat6_map_113017.csv"), header=FALSE, sep=",")
names(adfgstat6)[1] <- "path_lon"
names(adfgstat6)[2] <- "path_lat"
names(adfgstat6)[3] <- "centroid_lon"
names(adfgstat6)[4] <- "centroid_lat"
names(adfgstat6)[5] <- "ADFGstat6"
names(adfgstat6)[6] <- "Piece"
subadfgstat6 <- adfgstat6[adfgstat6$ADFGstat6 %in% subsettt, ]
world <- map_data("world")

#extent of map
minlon <- min(subadfgstat6$path_lon)*1.001 #lon negative
maxlon <- max(subadfgstat6$path_lon)*0.999
minlat <- min(subadfgstat6$path_lat)*0.999
maxlat <- max(subadfgstat6$path_lat)*1.001

# get summary statistics
subadfgstat6$ADFGstat6fac <- as.factor(subadfgstat6$ADFGstat6)
avgcatchdat <- ddply(dat, .(ADFGstat6),summarise, Averagecatch =
    mean(Weight/1000))
obscountdat <- ddply(dat, .(ADFGstat6),summarise, observations =
    length(Weight))  
startobscountdat <- ddply(dat, .(lag.adfgstat6),summarise, observations =
    length(Weight))  
names(startobscountdat)[1] <- "ADFGstat6"    
subadfgstat6 <- merge(subadfgstat6, avgcatchdat, by = c("ADFGstat6"),
    all.x=TRUE,all.y=TRUE)
subadfgstat61 <- merge(subadfgstat6, obscountdat, by = c("ADFGstat6"),
    all.x=TRUE,all.y=TRUE)
subadfgstat6 <- rbind(subadfgstat61)

#this reorders pieces 1 and 2 (otherwise maps polygons wrong)
subadfgstat6$Piece[subadfgstat6$Piece==1] <- 3
subadfgstat6$Piece[subadfgstat6$Piece==2] <- 1
subadfgstat6$Piece[subadfgstat6$Piece==3] <- 2
subadfgstat6$Piece <- as.factor(subadfgstat6$Piece)

mappedsub <- ggplot(data=subadfgstat6) +
geom_path(aes(x=path_lon,y=path_lat, group = ADFGstat6), color="black", 
    size=0.375) + 
geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region), fill="grey", 
    color="black", size=0.375) +
geom_polygon(data=subadfgstat61,aes(x=path_lon,y=path_lat,group = ADFGstat6fac, 
    fill = Averagecatch), color="black",size=0.375) +
xlim(minlon,maxlon) + ylim(minlat,maxlat) + 
scale_fill_gradient(name="Average catches\n(metric tons)", low = "#e8e8e8", 
    high = "#000000") +
geom_text(aes(x=centroid_lon,y=centroid_lat,
    label = observations, color=factor(Piece)), size=3, show.legend=F, 
    color = "white") +
geom_point(aes(x=0,y=0,color=factor(Piece)),size = 0) +
    guides(color = "none") +
    theme_bw() +
    scale_color_grey(start = 0.2, end = 0.8) +
theme(text = element_text(size=10), legend.title=element_text(lineheight=.6), 
    axis.title.y=element_text(vjust=1.5), legend.position=c(.15,.25)) + 
xlab("Longitude") + ylab("Latitude")  

ggsave(mappedsub, file = paste0(dirout,
    "Figure_4",".eps"), width = 7.5, height = 4.21875, dpi=800)

# make polynomial plots    
FIML <- as.data.frame(results$correction$OutLogit)
uncorrect <- as.data.frame(results$uncorrected$OutLogit)
colnames(FIML) <- c("coef","se","tstat")
colnames(uncorrect) <- c("coef","se","tstat")

z1 <- rep(1,101)
z2 <- seq(0,1,0.01)
z3 <- seq(0,1,0.01)^2
z4 <- seq(0,1,0.01)
z5 <- seq(0,1,0.01)

#find standard error using delta method at different points
polylistSE <- list()
for (i in 1:kk) {

seout <- list()

for (a in 1:length(z1)) {

zz1 <- z1[a]
zz2 <- z2[a]
zz3 <- z3[a]
zz4 <- z4[a]
zz5 <- z5[a]

form <-
sprintf("~((x%.0f*(%f) + x%.0f*(%f) + x%.0f*(%f))*(1-(1-exp(-(%f)/(1-(%f))))))", 
    (1+kk+i),zz1,(1+2*kk+i),zz2,(1+3*kk+i),zz3,zz4,zz5)
       
seout[[a]] <- deltamethod(as.formula(form), 
    results$correction$OutLogit[1:numparams,1], 
    results$correction$H1)
    
}

polylistSE[[i]] <- unlist(seout)

}

# make polynomials from results
polynomials = matrix(FIML[(1+kk+1):(1+4*kk),1], kk)
polylist <- list()
for (i in 1:kk) {
polylist[[i]] <- (polynomials[i,1]*rep(1,101) + polynomials[i,2]*seq(0,1,0.01) + 
    polynomials[i,3]*seq(0,1,0.01)^2)*
    (1-(1-exp(-seq(0,1,0.01)/(1-seq(0,1,0.01)))))
}

forplot <- data.frame(cbind(unlist(polylist), unlist(polylistSE), 
    rep(seq(0,1,0.01), kk), rep(1:kk, each=101)))
forplot$chosenlocjj <- forplot$X4

forplot <- merge(forplot, unique(finaldata2[,c("chosenlocjj", "ADFGstat6")]), 
    by = c("chosenlocjj"), all.x = TRUE, all.y = TRUE)

forplot$issig <- abs(forplot$X1)/(forplot$X2) > 1.96
 
polyline = ggplot(data=forplot) + 
    geom_line(aes(y=X1,x=X3)) + 
    geom_point(aes(y=X1,x=X3), size = forplot$issig*1.5) +
    facet_wrap(~ADFGstat6, scales = "free", ncol = 2) +
    theme(text = element_text(size=8), axis.title.y=element_text(vjust=1.5)) + 
    xlab("Probability of location choice") + 
    ylab("Correction function")

testplot <- data.frame(cbind(finaldata2$chosenlocjj, results$probs))
testplot$chosenlocjj <- testplot$X1

testplot <- merge(testplot, unique(finaldata2[,c("chosenlocjj", "ADFGstat6")]), 
    by = c("chosenlocjj"), all.x = TRUE, all.y = TRUE)

histprobs = ggplot(data=testplot) + 
    geom_histogram(aes(x=X2), bins=10) +
    facet_wrap(~ADFGstat6, scales = "free_y", ncol = 2) +
    theme(text = element_text(size=8), axis.title.y=element_text(vjust=1.5)) + 
    xlab("Probability of location choice") + 
    ylab("Number of observations")

library(ggpubr)

figure2 <- ggarrange(polyline, histprobs,
    labels = c("A", "B"),
    ncol = 2, nrow = 1)

ggsave(paste0(dirout,
    "Figure_B1.eps"), figure2, 
    width = 7.5, height = 4.21875, dpi=800, device=cairo_ps)
    
# make welfare plot
# these adfg areas are the chinoook savings area
pol = as.numeric(as.character(unique(finaldata2$chosenlocjj[
    finaldata2$ADFGstat6 %in% 
    c("645530","645501","655500","665500","665430","655430")])))

# get utility matrix
offsetu = 0
FIMLu = exp(FIMLrev+FIMLcost+offsetu)
uncorrectu = exp(uncorrectrev+uncorrectcost+offsetu)

# utilities with and without policy area
# prep for plot
ggdata = as.data.frame(cbind(rbind(as.matrix(log(rowSums((FIMLu))) - 
    log(rowSums((FIMLu[-pol])))),
	as.matrix(log(rowSums((uncorrectu))) - log(rowSums((uncorrectu[-pol]))))),
	rbind(as.matrix(rep("FIML",dim(FIMLu)[1])), as.matrix(rep("Uncorrected",
    dim(FIMLu)[1])))))
colnames(ggdata) <- c("Udiff","Method")
ggdata$Udiff = as.numeric(as.character(ggdata$Udiff))
ggdata$Uperc = ggdata$Udiff/
    as.numeric(abs(rbind(as.matrix(log(rowSums((FIMLu)))),
	as.matrix(log(rowSums((uncorrectu)))))))

finaldata2$ID = as.factor(paste(finaldata2$Age, finaldata2$Horsepower, 
    finaldata2$Grosstons, sep="_"))
    
ggdata = cbind(ggdata,rbind(finaldata2[,c("Age","Horsepower","Grosstons","ID")],
    finaldata2[,c("Age","Horsepower","Grosstons","ID")]))

# function to get quantiles
qfun <- function(x, q = 4) {
    quantile <- cut(x, breaks = quantile((x), probs = 0:q/q), 
        include.lowest = TRUE, labels = 1:q)
    quantile
}

# subset data by vessel characteristic quantiles
ggdata = ddply(ggdata, .(Method), mutate, Agequant = qfun(Age))
ggdata = ddply(ggdata, .(Method), mutate, Horsepowerquant = qfun(Horsepower))
ggdata = ddply(ggdata, .(Method), mutate, Grosstonsquant = qfun(Grosstons))
ggdata$Horsepowerquant = mapvalues(ggdata$Horsepowerquant, from = 
    c("1", "2", "3", "4"), 
to = c("0-25%", "25-50%", "50-75%", "75-100%"))

ggdata$startlocjj =  rbind(as.matrix(Ycol$startlocjj),
    as.matrix(Ycol$startlocjj))
    
ggdata = ggdata[!(ggdata$startlocjj %in% pol),]

mu <- ddply(ggdata, c("Method", "Horsepowerquant"), summarise, 
    grp.mean=mean(Uperc))
    
welfareviolin = ggplot(data=ggdata, aes(fill = Horsepowerquant)) + 
    geom_boxplot(aes(y=Uperc*100,x=Method), 
        coef = 0, 
        outlier.shape = NA, 
        notch=TRUE) + 
    coord_cartesian(ylim = quantile(ggdata$Uperc*100, c(0.075, 0.925))) +
    theme_bw() +
    scale_fill_grey(start = 0.2, end = 0.8) +
    theme(text = element_text(size=10), axis.title.y=element_text(vjust=1.5),
        axis.title.x=element_blank()) + 
    theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
    ) + 
    labs(fill='Horsepower') +
    ylab("Percentage Welfare Loss")

ggsave(welfareviolin, file = paste0(dirout,
    "Figure_6",".eps"), width = 7.5, height = 4.21875, dpi=800)

# make maps with predicted catches
subadfgstat6$ADFGstat6fac = paste(subadfgstat6$ADFGstat6, subadfgstat6$Piece, 
    sep="_")

#make into metric tons (undo scaling from estimates)
mincatch <- 1/100
# make predicted weights with both methods    
FIMLweights = cbind(ddply(finaldata2, .(ADFGstat6,chosenlocjj), summarise, 
    rating.mean=length(chosenlocjj)),
    ((colMeans(FIMLrev)/FIML[1,1])/
    mincatch))
Uncorrectweights = cbind(ddply(finaldata2, .(ADFGstat6,chosenlocjj), summarise, 
    rating.mean=length(chosenlocjj)),
    ((colMeans(uncorrectrev)/uncorrect[1,1])/
    mincatch))    
names(FIMLweights)[3] <- "numsamples"
names(FIMLweights)[4] <- "RelativeCatch"
FIMLweights$Method = "FIML"
names(Uncorrectweights)[3] <- "numsamples"
names(Uncorrectweights)[4] <- "RelativeCatch"
Uncorrectweights$Method = "Uncorrected"

# merge predicted catches into spatial polygon data
subadfgstat6$orderobs <- 1:dim(subadfgstat6)[1]
tempadfg1 <- merge(subadfgstat6,FIMLweights,by=c("ADFGstat6"))
tempadfg2 <- merge(subadfgstat6,Uncorrectweights,by=c("ADFGstat6"))
tempadfg1 =  tempadfg1[with(tempadfg1, order(orderobs)), ]
tempadfg2 =  tempadfg2[with(tempadfg2, order(orderobs)), ]
subadfgstat6catch = rbind(tempadfg1, tempadfg2) 

#double up map data for facet
world2 = world
world2$Method = "Uncorrected"
world$Method = "FIML"
world = rbind(world,world2)

#chinook policy area
polids <- c("645530","645501","655500","665500","665430","655430")
polshapeori = ((adfgstat6[adfgstat6$ADFGstat6 %in% 
    polids[(polids %in% unique(subadfgstat6$ADFGstat6))],]))

# make polygon of policy area
hulls <- concaveman(as.matrix(polshapeori[,1:2]), concavity = 2, 
    length_threshold = 0)
hulls <- as.data.frame(hulls)
hulls$path_lon <- hulls$V1
hulls$path_lat <- hulls$V2
hulls <- hulls[-36,]
polshape = hulls
polshape2 = hulls
polshape2$Method = "Uncorrected"
polshape$Method = "FIML"
polshape = rbind(polshape,polshape2)
polshape$Policy = "CSSA"

mappedsub <- ggplot() +
    geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region), 
        color="black", size=0.375) +
    geom_path(data=subadfgstat6catch,aes(x=path_lon,y=path_lat,
        group = ADFGstat6fac), color="black", size=0.375) +
    geom_polygon(data=subadfgstat6catch,aes(x=path_lon,y=path_lat,
        group = ADFGstat6fac, fill = RelativeCatch), color="black",size=0.375) +
    geom_polygon(data=polshape,aes(x=path_lon,y=path_lat,colour=Policy), 
        alpha = 0,
        size=1.75) +
    geom_text(data=subadfgstat6catch,aes(x=centroid_lon,y=centroid_lat,
        label = round(RelativeCatch,2)),size=2.0, 
        color="#cccccc",
        position=position_dodge(width = 5)) +
    scale_color_manual(values = "black") +
    facet_wrap(~Method) +
    labs(fill='Metric Tons') +
    scale_fill_gradient(low = "white", high = "black") +
    xlim(minlon,maxlon) + ylim(minlat,maxlat) + 
    theme_bw() +
    theme(text = element_text(size=10), axis.title.y=element_text(vjust=1.5),
    legend.position = c(.90, .79),
    legend.text=element_text(size=7.5),
    legend.title=element_text(size=10),
    legend.key.size = unit(0.25, "cm"),
    legend.spacing.y = unit(0.10, 'cm')) + 
    xlab("Longitude") + ylab("Latitude")

ggsave(mappedsub, file = paste0(dirout,
    "Figure_5",".eps"), width = 7.5, height = 4.21875, dpi=800)
    