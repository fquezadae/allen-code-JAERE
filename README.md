Barebones FishSET
=========
---

This is a barebones version of the FishSET package, primarily used for directly 
calling some of the lower-level functions.

# Installation #
---

To install with vignette:

    > install.packages("devtools")
	> library(devtools)
	
Set the directory you've downloaded the package into:

    > setwd("J:/Fishperson/Directory_containing_barebones.FishSET")

Install with vignettes:

    > install("barebones.FishSET", build_vignettes = TRUE)
	
Check out documentation and vignette:

    > help(package="barebones.FishSET")
    > vignette("barebones.FishSET.vignette", package="barebones.FishSET")

Please see vignette for usage examples.

You can create a stand-alone .html vignette by then running:

	> setwd("barebones.FishSET")
	> devtools::build_vignettes()
    
Reproducible example
-----------------

The code used to reproduce the work in Chen et al. can be found below. This 
also serves as an use example of the package.

``` r
library(barebones.FishSET)
library(dplyr)
library(doParallel)
library(foreach)
 
# check how many workers your machine has
registerDoParallel(cores = 8) 
    
# number of observations per location
iinum <- rep(1000,4)
# number of MC iterations
mcnum <- 100
# no constant model no kernel
constp <- rep(0,4)
bwnum <- -1

# GEV choice error
sertype = "gev"
chodev = 1

# utility and catch parameters
alpha <- 3
betac <- -1
regconstant <- 0

# set a directory here
dirout <- paste0("U:\\AFSC_data_code\\",
    "fiml_correction_runs\\dahlMCs\\")

for (it in seq(1,4,by=0.5)) {

# run a model for each catch devation for figure 3
cdev = it

print(Sys.time()) #about 3 hours each
par_res_base <- foreach(nn = 1:(mcnum), .packages = c("barebones.FishSET", 
    "dplyr")) %dopar% mc_func_base(nn, constp, cdev, iinum, bwnum, sertype, 
    chodev, alpha, betac, regconstant)
print(Sys.time())

saveRDS(par_res_base, paste0(dirout, "testMC2\\MC_",cdev, ".rds"))

}

cdev = 0

print(Sys.time())
par_res_base <- foreach(nn = 1:(mcnum), .packages = c("barebones.FishSET", 
    "dplyr")) %dopar% mc_func_true(nn, constp, cdev, iinum, bwnum, sertype, 
    chodev, alpha, betac, regconstant)
print(Sys.time())

saveRDS(par_res_base, paste0(dirout, "testMC\\MC_",cdev, ".rds"))

```