#simrep.R
#Runs the holoSimCell model for Ash using an ENM to define habitat suitability

#Set the "node" you're simulating 1-500 possible
#Set the number of replicates to simulate
#Add the user's initials to keep this file separate from others...
#Usage: Rscript simrep.R 1 10 JDR
args <- commandArgs(TRUE)
i <- as.numeric(args[1])
nreps <- as.numeric(args[2]) 
who <- as.character(args[3])  

if(length(args) == 0){
    i <- 1
    nreps <- 10
    who <- "JDR"
}

setwd("../OUT")   

library(holoSimCell)

### imputed, popmap (individualID->pop mapping), pts (sample locations) and ashpred
### are now built into holoSimCell
### as built in dataframes (in data/ directory)

rownames(popmap) <- popmap[,1]
table(popmap[gsub("fp","",names(imputed)),2])
imputed.pruned=imputed[,-which(gsub("fp","",names(imputed))%in%popmap[popmap$abbrev=="Michigan","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="UNK","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="MO1","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="ON1","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="VA1","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="MB1","id"])]
removes <- c()
popids <- popmap[gsub("fp","",names(imputed.pruned)),2]
table(popids)
for (a in unique(popids))
{
    if (sum(popids==a)>14)
    {
        removes <- c(removes,sample(which(popids==a),1))
    }
}
imputed.pruned <- imputed.pruned[,-1*removes]

poptbl <- table(popmap[gsub("fp","",names(imputed.pruned)),2])

samppts <- pts[pts$abbrev %in% names(poptbl),]

for(repl in 1:nreps) {

    #Draw parms
    parms <- drawParms(control = system.file("extdata/csv","priors.csv",package="holoSimCell"))

    ##this should produce a landscape with 20 x 18 (x,y) _square_ cells that also have
    ##21empirical samples in separate grid cells (otherwise need to figure out something else)
    ## ashland is a stored R object as well
    if (!exists("ashland"))
    {

    } else landscape <- ashland

    ph = getpophist.cells(hab_suit=landscape,samptime=1,refs=parms$refs,refsz=parms$ref_Ne,
                      mix=parms$mix,
                      shortscale=parms$shortscale,shortshape=parms$shortshape,
                      longmean=parms$longmean,
                      sz=parms$sz,
                      K=parms$Ne)
    if(FALSE) {
        plothist(ph)
    }

    ###########################################################################
    ##Setting up stuff for fastsimcoal down here
    ##There may be better places to store this (in the new ph object?)
    loc_parms <- data.frame(marker = "snp",
                        nloci = parms$nloci,           
                        seq_length = parms$seq_length,
                        mu = parms$mu)

    loc_parms2 <- loc_parms     #Will use this object for simulations, increase # of loci until we get loc_parms$nloci SNPs
    loc_parms2$nloci <- round(loc_parms2$nloci*2) 	#Simulate more loci than needed to account for monomorphic sites

    preLGMparms <- data.frame(preLGM_t = parms$preLGM_t/parms$G,		#Time / GenTime
                    preLGM_Ne = parms$preLGM_Ne,
                    ref_Ne = parms$ref_Ne)

    parms_out <- as.data.frame(c(ph$struct[which(!names(ph$struct) %in% names(parms))], parms))

    #With smaller K, some populations have very very low N at the end of the simulation
    #In those cases, we need to inflate N a bit for the coalescent simulation
    ph$Nvecs[ph$Nvecs[,702] > 0 & ph$Nvecs[,702] < 1,702] <- 1


    #Run the coalescent simulation
    fscout <- NULL
    while(is.null(fscout)) {
        tmp <- runFSC_step(ph = ph, 
						l = landscape,
						sample_n = poptbl[1],
						preLGMparms = preLGMparms,
						label = paste0("Ash_", round(runif(1),5)),
						delete_files = TRUE,
						num_cores = 1,
						exec = "fsc26",
						loc_parms = loc_parms2,
						found_Ne = parms$found_Ne)

        SNPnum <- get.gSum(tmp)
        message(paste("Coalescent simulation complete:", SNPnum$nvar, "variable sites"))

        #Check to see if there are enough variable sites in the dataset...
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
    #Treated as a factor though, do not index with this padded value
    oldID <- sort(unique(c(ph$coalhist$src, which(ph$Nvecs[,1] > 0))))
    newID <- c(1:length(oldID))-1
    sample_pops <- plyr::mapvalues(landscape$sampled, oldID, newID, warn_missing = FALSE)

    #This is simpler, more reliable
    #also seems to be some error check in strataG that was causing problems (none of the specified ids were found in the specified strata)
    #nameStrat is deprecated
    #fscout <- nameStrat(fscout = fscout, pops = ph, sample_pops = sample_pops, sample_n = as.vector(poptbl))
    fscout@data$strata <- paste0("pop-", sapply(strsplit(fscout@data$ids,"_"), function(x){x[1]}))
    #The line above is probably more flexible than it looks, FSC gives the fsc popid as the first part of the individual name  

    #Build popDF for stat calculation
    strat_order <- order(as.character(sample_pops))
    popDF <- data.frame(id = unique(fscout@data$strata),
                            grid.cell = sample_pops[strat_order],
                            sample.size = as.vector(poptbl)[strat_order],
                            col = ph$pophist$col[sample_pops[strat_order]],
                            row = ph$pophist$row[sample_pops[strat_order]])
    popDF <- popDF[order(popDF$grid.cell),]


    stats_out <- holoStats(out = fscout, 
                       popDF = popDF,
                       extent = c(ph$struct["xdim"], ph$struct["ydim"]), 
                       cores = 1) 

    BVmax <- max(bioticVelocity(ph, metrics = "centroid")$centroidVelocity)
    
    #!# IS THIS OUTPUT IN THE FORMAT THAT WE WANT??
    all_out <- c(date = date(), i, repl, parms_out, BVmax, stats_out)
    all_out$refs <- paste(eval(parse(text=as.character(ph$struct["refs"]))), collapse = ".")

    #Write a file, if none exists, or append to file, if present
    if(!file.exists(paste0("Ash_", i, "_", who, ".csv"))) {
        write.table(all_out, paste0("Ash_", i, "_", who, ".csv"), sep = ",", quote = FALSE, row.names = FALSE)
    } else {
        write.table(all_out, paste0("Ash_", i, "_", who, ".csv"), quote = FALSE, row.names = FALSE, sep=",", append = TRUE, col.names = FALSE)
    }

}

