###
### script to run simulations for evaluation of spatial summary stats
### GRID not individuals 
### STEP change growth model
###

library(holoSimCell)

#Run these in $TMPDIR on the cluster so they don't step on one another
codedir = "/mnt/research/TIMBER/spatsimsEW/code"
simdir = system("echo $TMPDIR", intern = TRUE)
outdir = "/mnt/research/TIMBER/spatsimsEW/OUT"

#Set the row of the scenarios table you're working with 
#and the number of replicates using command line arguments
#Usage: Rscript sim_spatstat.R 1 1 100
args <- commandArgs(TRUE)
i <- as.character(args[1])
nreps <- as.numeric(args[2])   

if(length(args) == 0){
    i <- 1
    repl <- 1
    nreps <- 10
}

#Read in the scenarios file (alternative parameter combinations for spatial stat evaluation)
setwd(codedir)

parmsets <- read.csv("scenarios.csv", header = TRUE)
parms <- parmsets[i,]

loc_parms <- data.frame(marker = "snp",
                        nloci = parms$nloci,           
                        seq_length = parms$locus_length,
                        mu = parms$mu)

loc_parms2 <- loc_parms     #Will use this object for simulations, increase # of loci until we get loc_parms$nloci SNPs
loc_parms2$nloci <- round(loc_parms2$nloci*1.25) 	#Simulate more loci than needed to account for monomorphic sites

preLGMparms <- data.frame(preLGM_t = parms$preLGM_t/parms$G,
                    preLGM_Ne = parms$preLGM_Ne,
                    ref_Ne = parms$ref_Ne)

#Set the bounds of the landscape
nrow = 15
ncol = 15

sample_pops <- seq(1, nrow*ncol, 11)
sample_n <- rep(25,length(sample_pops))

#Loop over the replicates 
for(repl in 1:nreps) {

    setwd(simdir)

    pops <- getpophist.cells(h=nrow*ncol,  ##demography, num habitats - product of xdim & ydim
                         xdim=ncol,        ##num cols - fixed at 15
                         ydim=nrow,        ##num rows - fixed at 15
                         maxtime=parms$texp/parms$G,  ##num time clicks to simulate (yrs, decades cents?) (call em decades here)  #!# reduced this a bit for testing, from 12000
                         lambda=1.25, ##intrisic rate of growth (discrete-time...right? can be modified with deltLambda)  #!# increased this a bit to avoid extinctions late in the sim
                                         #deltLambda: see function definition
                         K=parms$Ne,        ##carry capacity (can be modified with deltK)
                                         #deltK: see function definition
                         refs=eval(parse(text=as.character(parms$refs))),    ##pop ids of refugia
                         refsz=parms$ref_Ne, ##sizes of refugia indicated in refs   #!# increased this from 10 to avoid early extinction due to dem. stoch.
                         sz=1,   ##number of spatial units per grid cell
                                        #dispersal traits
                         distance.fun=distancePDF, #takes a length and 5 parameters and creates mig matrix
                         shortscale=parms$shortscale,
                         longmean=parms$longmean,
                         shortshape=1,    #shape of short-distance
                         mix=parms$mix,         #proportion of LDD
                         popDispInfl=function(x){log(x+1)},
                         samptime = 1,  #!# sample every X generations
                         CVn = NULL,       #!# coefficient of variation X% 
                         pois.var = FALSE, #!# Poisson distribution for demographic stochasticity
                         extFUN = NULL,
                         hab_suit = NULL
                         )

    parms_out <- c(parms, loc_parms2, preLGMparms)

    #Run fastsimcoal in a while loop
    #Only stops when loc_parms$nloci SNPs or more are in the output
    fscout <- NULL
    while(is.null(fscout)) {
        tmp <- runFSC_step(popsobj = pops,  		#Pass the 'pops' object from getpophist.cells()
					preLGMparms = preLGMparms,	    #This has parms for the refuge, preLGM size and timing
					sample_pops = sample_pops,	    #The IDs of sampled populations
					sample_n = sample_n,		    #The sample size per population
					label = paste0("simscen_",i,round(runif(1),5)),	#Label for FSC simulation files
					delete_files = TRUE,		    #Logical - clear out .par, .arp, and other FSC outputs?
					num_cores = 1,				    #Number of processors to use for FSC
					exec = "fsc26",				    #Executable for FSC (needs to be in $PATH variable)
					loc_parms = loc_parms2,		    #Vector of locus parameters
					found_Ne = parms$foundNe)		#Founding size 

        SNPnum <- get.gSum(tmp)
        message(paste("Coalescent simulation complete:", SNPnum$nvar, "variable sites"))

        #This is a check to see if there are enough variable sites in the dataset...
        if(SNPnum$nvar > loc_parms$nloci) {
            variableloci <- which(SNPnum$nall == 2)
            keeploci <- sample(variableloci,loc_parms$nloci,replace = FALSE)
            fscout <- tmp[,keeploci]
            rm(tmp)
        } else if(SNPnum$nvar == loc_parms$nloci) {
            fscout <- tmp
            rm(tmp)
        } else {
            message("too few SNPs, doubling number of simulated loci and trying again")
            loc_parms2$nloci <- loc_parms2$nloci*2
            rm(tmp)
        }
    }

    #Naming the strata - uses grid cell # with left-padding to make all the same length
    #Treated as a factor through, do not index with this padded value
    fscout <- nameStrat(fscout = fscout, pops = pops, sample_pops = sample_pops, sample_n = sample_n)

    #Build popDF for stat calculation
    popDF <- data.frame(id = unique(fscout@data$strata),
                                grid.cell = sample_pops,
                                sample.size = sample_n,
                                col = pops$pophist$col[sample_pops],
                                row = pops$pophist$row[sample_pops])


    stats_out <- holoStats(out = fscout, 
                            popDF = popDF,
                            extent = c(pops$struct["xdim"], pops$struct["ydim"]), 
                            cores = 1) 

    #!# IS THIS OUTPUT IN THE FORMAT THAT WE WANT??
    all_out <- c(rep = repl, date = date(), parms_out, stats_out)
    all_out <- all_out[c(3,1,2,4:length(all_out))]
    all_out$refs <- paste(eval(parse(text=as.character(parms$refs))), collapse = ".")

    #Write a file, if none exists, or append to file, if present
    setwd(outdir)

    if(!file.exists(paste0("SpatStats_parmset-",i,".csv"))) {
        write.table(all_out, paste0("SpatStats_parmset-",i,".csv"), sep = ",", quote = FALSE, row.names = FALSE)
    } else {
        write.table(all_out, paste0("SpatStats_parmset-",i,".csv"), quote = FALSE, row.names = FALSE, sep=",", append = TRUE, col.names = FALSE)
    }

}
