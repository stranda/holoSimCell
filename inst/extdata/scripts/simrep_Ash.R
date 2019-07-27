#simrep.R
#Runs the holoSimCell model for Ash using an ENM to define habitat suitability

#Set the "node" you're simulating 1-500 possible
#Set the number of replicates to simulate
#Add the refuge location (grid cell ID)
#Add the user's initials to keep this file separate from others...
#Usage: Rscript simrep.R 1 10 5 JDR

args <- commandArgs(TRUE)
i <- as.numeric(args[1])
nreps <- as.numeric(args[2]) 
who <- as.character(args[3])  
#refs <- as.character(args[3])  #!!# If we want to have the refuge location passed as part of the command line 
#who <- as.character(args[4])

if(length(args) == 0){
    i <- 1
    nreps <- 1
    who <- "JDR"
    #refs <- 5
}

simdir <- system("echo $TMPDIR", intern = TRUE)
outdir <- "/mnt/research/TIMBER/Ash/OUT" 

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


for(repl in 1:nreps) {

    landscape <- ashpollen
    
    #Draw parms
    parms <- drawParms(control = system.file("extdata/csv","priors.csv",package="holoSimCell"))
    #parms$refs <- refs   #!!# If we want to have the refuge location passed as part of the command line
    parms$mu <- 1e-6  #!!!# Bumping up mutation rate!!  Does this help with speeds?

    if (FALSE) #don't run any without suitability
    {
                                        #Logical parameter of teh simulation, use hab_suit or not...
        if(parms$use.hab_suit == 0) {
                                        #landscape = NULL   #Don't do it this way, entire matrix is habitable
            landscape$hab_suit[!is.na(landscape$hab_suit)] <- 1  #This way ignores the glacier
                                        #landscape$hab_suit[landscape$hab_suit > 0] <- 1   #This way maintains the 0 suitability for glaciated cells
        }
    } else { #but set the model number to 2 for use.hab_suit
            parms$use.hab_suit = 2
    }

    
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
    loc_parms2$nloci <- round(loc_parms2$nloci*1.25) 	#!!# #Simulate more loci than needed to account for monomorphic sites

    preLGMparms <- data.frame(preLGM_t = parms$preLGM_t/parms$G,		#Time / GenTime
                    preLGM_Ne = parms$preLGM_Ne,
                    ref_Ne = parms$ref_Ne)

    parms_out <- as.data.frame(c(ph$struct[which(!names(ph$struct) %in% names(parms))], parms))

    #With smaller K, some populations have very very low N at the end of the simulation
    #In those cases, we need to inflate N a bit for the coalescent simulation
    ph$Nvecs[ph$Nvecs[,702] > 0 & ph$Nvecs[,702] < 1,702] <- 1


    #Run the coalescent simulation
    setwd(simdir) 
    fscout <- NULL
    SNPnum <- data.frame(nvar = 0)
    MAF_filt <- 0.01
    while(SNPnum$nvar < loc_parms$nloci) {
        message(paste("Coalescent simulation with",loc_parms2$nloci,"loci"))
        tmp <- runFSC_step(ph = ph, 
                        l = landscape,
                        sample_n = poptbl[1],
                        preLGMparms = preLGMparms,
                        label = paste0("Ash_", round(runif(1, min=1, max=1e6))),
                        delete_files = TRUE,
                        num_cores = 2,
                        exec = "fsc26",
                        loc_parms = loc_parms2,
                        found_Ne = parms$found_Ne)

        if(is.null(fscout)) {
            fscout <- tmp
        } else {
            fscout@data <- cbind(fscout@data, tmp@data[,-c(1,2)])
            colnames(fscout@data) <- c(colnames(fscout@data)[1:2], paste0("Locus_",c(1:(dim(fscout@data)[2]-2))))        
        }
    
        rm(tmp)
    
        SNPnum <- get.gSum(fscout)
        message(paste("Coalescent simulation complete:", SNPnum$nvar, "variable sites"))
        variableloci <- which(SNPnum$nall == 2)
        fscout <- fscout[,variableloci]
        colnames(fscout@data) <- c(colnames(fscout@data)[1:2], paste0("Locus_",c(1:(dim(fscout@data)[2]-2))))

        pass_MAF_filt <- which((apply(fscout@data[,-c(1,2)],2,table)/dim(fscout@data)[1])[2,] > MAF_filt)
        fscout <- fscout[,pass_MAF_filt]
        colnames(fscout@data) <- c(colnames(fscout@data)[1:2], paste0("Locus_",c(1:(dim(fscout@data)[2]-2))))

        SNPnum <- get.gSum(fscout)
        message(paste(SNPnum$nvar, "variable sites passing MAF filter"))

        if(SNPnum$nvar < loc_parms$nloci) {
            loc_parms2$nloci <- max(3*(loc_parms$nloci - SNPnum$nvar),100)  #!!!# Increasing minimum number of simulations to 100 (from 50)
        }
      
    }

    if(SNPnum$nvar > loc_parms$nloci) {
        variableloci <- which(SNPnum$nall == 2)
        keeploci <- sample(variableloci,loc_parms$nloci,replace = FALSE)
        fscout <- fscout[,keeploci]
        colnames(fscout@data) <- c(colnames(fscout@data)[1:2], paste0("Locus_",c(1:(dim(fscout@data)[2]-2))))
    }

    #Naming the strata - uses grid cell # with left-padding to make all the same length
    #Treated as a factor though, do not index with this padded value
    oldID <- sort(unique(c(ph$coalhist$src, which(ph$Nvecs[,1] > 0))))
    newID <- c(1:length(oldID))-1
    sample_pops <- plyr::mapvalues(landscape$sampled, oldID, newID, warn_missing = FALSE)

    #This is simpler, more reliable - well, not exactly...
    #also seems to be some error check in strataG that was causing problems (none of the specified ids were found in the specified strata)
    #nameStrat is deprecated
    #fscout <- nameStrat(fscout = fscout, pops = ph, sample_pops = sample_pops, sample_n = as.vector(poptbl))
    #fscout@data$strata <- paste0("pop-", sapply(strsplit(fscout@data$ids,"_"), function(x){x[1]}))
    #The line above is probably more flexible than it looks, FSC gives the fsc popid as the first part of the individual name  

    #Build popDF for stat calculation
    #strat_order <- order(as.character(sample_pops))
    #popDF <- data.frame(id = unique(fscout@data$strata),
    #                        grid.cell = sample_pops[strat_order],
    #                        sample.size = as.vector(poptbl)[strat_order],
    #                        col = ph$pophist$col[sample_pops[strat_order]],
    #                        row = ph$pophist$row[sample_pops[strat_order]])
    #popDF <- popDF[order(popDF$grid.cell),]

    #Build popDF from landscape$sampdf
    popDF <- landscape$sampdf
    colnames(popDF)[colnames(popDF)=="abbrev"] <- "id"
    popDF$grid.cell <- plyr::mapvalues(popDF$cell, oldID, newID, warn_missing = FALSE)
    popDF$sample.size <- as.vector(poptbl[match(popDF$id,names(poptbl))])
    
    stats_out <- holoStats(out = fscout, 
                       popDF = popDF,
                       extent = c(ph$struct["xdim"], ph$struct["ydim"]), 
                       cores = 1) 

    BVmax <- max(bioticVelocity(ph, metrics = "centroid")$centroidVelocity)
    
    #!# IS THIS OUTPUT IN THE FORMAT THAT WE WANT??
    all_out <- c(date = date(), node=i, rep=repl, parms_out, BVmax=BVmax, stats_out)
    all_out$refs <- paste(eval(parse(text=as.character(ph$struct["refs"]))), collapse = ".")

    #Write a file, if none exists, or append to file, if present
    setwd(outdir)
    if(!file.exists(paste0("Ash_", i, "_", who, ".csv"))) {
        write.table(all_out, paste0("Ash_", i, "_", who, ".csv"), sep = ",", quote = FALSE, row.names = FALSE)
    } else {
        write.table(all_out, paste0("Ash_", i, "_", who, ".csv"), quote = FALSE, row.names = FALSE, sep=",", append = TRUE, col.names = FALSE)
    }

}

