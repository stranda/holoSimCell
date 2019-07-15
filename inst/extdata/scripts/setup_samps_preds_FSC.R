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



##this should produce a landscape with 20 x 18 (x,y) _square_ cells that also have
##21empirical samples in separate grid cells (otherwise need to figure out something else)
## ashland is a stored R object as well
if (!exists("ashland"))
{
    landscape <- def_grid_pred(pred=ashpred[,,701:1],samppts=samppts,init.ext=c(18,18),keep.thresh=0.05)
} else landscape <- ashland

ph = getpophist.cells(hab_suit=landscape,samptime=1,refs=(5),refsz=100,
                      mix=0.01,
                      shortscale=10,shortshape=1,
                      longmean=70,
                      sz=150)
plothist(ph)

###########################################################################
##Setting up stuff for fastsimcoal down here
##There may be better places to store this (in the new ph object?)
loc_parms <- data.frame(marker = "snp",
                        nloci = 100,           
                        seq_length = 80,
                        mu = 1e-7)

loc_parms2 <- loc_parms     #Will use this object for simulations, increase # of loci until we get loc_parms$nloci SNPs
loc_parms2$nloci <- round(loc_parms2$nloci*1.25) 	#Simulate more loci than needed to account for monomorphic sites

preLGMparms <- data.frame(preLGM_t = 100000/30,		#Time / GenTime
                    preLGM_Ne = 20000,
                    ref_Ne = 5000)

parms_out <- as.data.frame(c(ph$struct, loc_parms, preLGMparms))

#With smaller K, some populations have very very low N at the end of the simulation
#In those cases, we need to inflate N a bit for the coalescent simulation
if(FALSE) {
	ph$Nvecs[ph$Nvecs[,702] > 0 & ph$Nvecs[,702] < 2,702]
	ph$Nvecs[ph$Nvecs[,702] > 0 & ph$Nvecs[,702] < 2,702] <- 1
}

#Run the coalescent simulation
fscout <- NULL
while(is.null(fscout)) {
    tmp <- runFSC_step3(ph = ph, 
						l = landscape,
						sample_n = 14,
						preLGMparms = preLGMparms,
						label = "AshLand",
						delete_files = TRUE,
						num_cores = 1,
						exec = "fsc251",
						loc_parms = loc_parms2,
						found_Ne = 50)

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
fscout <- nameStrat(fscout = fscout, pops = ph, sample_pops = sample_pops, sample_n = as.vector(poptbl))

#Build popDF for stat calculation
popDF <- data.frame(id = unique(fscout@data$strata),
                            grid.cell = sample_pops,
                            sample.size = as.vector(poptbl),
                            col = ph$pophist$col[sample_pops],
                            row = ph$pophist$row[sample_pops])
popDF <- popDF[order(popDF$grid.cell),]


stats_out <- holoStats(out = fscout, 
                       popDF = popDF,
                       extent = c(ph$struct["xdim"], ph$struct["ydim"]), 
                       cores = 1) 

BVmax <- max(bioticVelocity(ph, metrics = "centroid")$centroidVelocity)
    
#!# IS THIS OUTPUT IN THE FORMAT THAT WE WANT??
all_out <- c(date = date(), parms_out, BVmax, stats_out)
all_out$refs <- paste(eval(parse(text=as.character(ph$struct["refs"]))), collapse = ".")

#Write a file, if none exists, or append to file, if present
if(!file.exists("AshLand.csv")) {
    write.table(all_out, "AshLand.csv", sep = ",", quote = FALSE, row.names = FALSE)
} else {
    write.table(all_out, "AshLand.csv", quote = FALSE, row.names = FALSE, sep=",", append = TRUE, col.names = FALSE)
}



