#Require these packages
require(cubature)
require(Rcpp)  
require(matrixStats)
require(tripack)
require(vegan)
require(geoR)
require(LDcorSV)
require(sna)
require(igraph)
require(compiler)
require(parallel)
require(strataG)



#Source in these files!
source("holoStats.R")
source("spatial_sumstats_functions_alvaradoserrano.R")
source("additional_stats_functions.R")
source("helper-functions.R")
source("segmentGLM.R")
source("runFSC_step.R") 
sourceCpp("helpers.cpp") 
source("integrate-mig-mat.R") 
source("plothist.R") 
source("pophist-cells2.R") 

#The line below still doesn't work...
#sourceCpp("helpers.cpp", cacheDir = "/mnt/research/TIMBER/spatsimsEW/cache/")

