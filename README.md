Barebones FishSET
=========
---

This is a barebones version of the FishSET package, primarily used for directly 
calling some of the lower-level functions. Please note an actively maintained 
repository can be found at: 
https://github.com/allen-chen-noaa-gov/barebones.FishSET.

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

The code used to reproduce the work in "Chen, Y. Allen, Alan C. Haynie, and 
Christopher M. Anderson. 2022. Journal of the Association of Environmental and 
Resource Economists." can be found below. This also serves as an use example of 
the package. Please note the package functions should be generalizable to any 
discrete choice problem - please see vignette for examples.

``` r
library(barebones.FishSET)
library(dplyr)
library(doParallel)
library(foreach)
 
# check how many workers your machine has
registerDoParallel(cores = #) 
    
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
dirout <- paste0(#)

for (it in seq(1,4,by=0.5)) {

# run a model for each catch devation for figure 3
cdev = it

print(Sys.time()) #about 3 hours each with 8 workers
par_res_base <- foreach(nn = 1:(mcnum), .packages = c("barebones.FishSET", 
    "dplyr")) %dopar% mc_func_base(nn, constp, cdev, iinum, bwnum, sertype, 
    chodev, alpha, betac, regconstant)
print(Sys.time())

saveRDS(par_res_base, paste0(dirout, "MC_", cdev, ".rds"))

}

# this model is for the base case with no private information
cdev = 0

print(Sys.time())
par_res_base <- foreach(nn = 1:(mcnum), .packages = c("barebones.FishSET", 
    "dplyr")) %dopar% mc_func_true(nn, constp, cdev, iinum, bwnum, sertype, 
    chodev, alpha, betac, regconstant)
print(Sys.time())

saveRDS(par_res_base, paste0(dirout, "MC_", cdev, ".rds"))

```

The code to reproduce the simulation figures is included as a script. These 
produce figures 1 through 3. Please see MC_make_figs.r in the exec folder.
The code to run the empirical example is also included as a script, however, 
the confidential data used to run the empirical example is not included. Please 
see Emp_ex.r in the exec folder.
