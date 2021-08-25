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