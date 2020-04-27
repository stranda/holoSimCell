#library(holoSimCell)
devtools::load_all()
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



##this should produce a landscape with (x,y) _square_ cells that also have
##21empirical samples in separate grid cells (otherwise need to figure out something else)
## ashland is a stored R object as well
if (!exists("ashland"))
{
###    ashland <- def_grid_pred(pred=ashpred[,,701:1],samppts=samppts,init.ext=c(40,36),keep.thresh=0.05)
    corners = matrix(c(-1702657,  1600000,
                       -1702657,  -1100000,
                       2908121,   -1100000,
                       2908212,  1600000),ncol=2,byrow=T)
    colnames(corners) <- c("x","y")
    rownames(corners) <- c("ul","ll","lr","ur")

    ashland <- def_grid_pred(pred=ashpred[,,701:1],samppts=samppts,init.ext=c(40,36),keep.thresh=0.05,
                            corners = corners
                             )
}

landscape <- ashland

###seed is based on time in seconds and the number of characters in the library path
###
###
sec=as.numeric(Sys.time())-1500000000
lp= as.numeric(as.character(nchar(paste(.libPaths(), collapse = " "))))
slp <- as.integer(floor(sec*lp))

set.seed(as.integer(sec))

ph = getpophist2.cells(hab_suit=landscape,
                       refs=(5),
                       refsz=100,
                       mix=0.001,  #note how small.
                       shortscale=6,  # scale parameter of weibull with shape below
                       shortshape=1, #weibull shape
                       longmean=75,  # mean of normal with sd = longmean
                       sz=150) #size of a cell (same units as longmean and shortscale)

gmap=make.gmap(ph$pophist,
               xnum=2, #number of cells to aggregate in x-direction
               ynum=2) #number of aggregate in the y-direction

ph2 <- pophist.aggregate(ph,gmap=gmap)

if(FALSE){
  pdf("aggregate_example.pdf")
  plothist(ph)
  plothist(ph2)
  dev.off()
}



outdir <- "~/Desktop"
simdir <- outdir
parms <- drawParms(control = system.file("extdata/csv","priors.csv",package="holoSimCell"))
parms$seq_length <- 80
parms$mu <- 1e-7 
loc_parms <- data.frame(marker = "snp",
                        nloci = parms$nloci,           
                        seq_length = parms$seq_length,
                        mu = parms$mu)


loc_parms2 <- loc_parms     #Will use this object for simulations, increase # of loci until we get loc_parms$nloci SNPs
loc_parms2$nloci <- round(loc_parms2$nloci*5) 	#!!# #Simulate more loci than needed to account for monomorphic sites

preLGMparms <- data.frame(preLGM_t = parms$preLGM_t/parms$G,		#Time / GenTime
                          preLGM_Ne = parms$preLGM_Ne,
                          ref_Ne = parms$ref_Ne)

parms_out <- as.data.frame(c(ph$struct[which(!names(ph$struct) %in% names(parms))], parms))

#With smaller K, some populations have very very low N at the end of the simulation
#In those cases, we need to inflate N a bit for the coalescent simulation
ph2$Nvecs[ph2$Nvecs[,702] > 0 & ph2$Nvecs[,702] < 1,702] <- 1


#Run the coalescent simulation
setwd(simdir) 

#For easy testing of runFSC_step_agg2() guts
if(FALSE) {
  phOLD <- ph
  ph <- ph2
  l <- landscape
  num_cores <- 1
  label <- "NewTest"
  exec <- "fsc26"
  sample_n <- 14
  found_Ne <- 50
}

out <- runFSC_step_agg2(ph = ph2,				#A new pophist object - (pophist, Nvecs, tmat, struct, hab_suit, coalhist)
                        l = landscape, 			#A new landscape object - (details, occupied, empty, sampled, hab_suit, sumrast, samplocsrast, samplocs)
                        sample_n = 14,		#Number of sampled individuals per population
                        preLGMparms = preLGMparms,		#This has parms for the refuge, preLGM size and timing
                        label = "test_agg",			#Label for FSC simulation files
                        delete_files = TRUE,	#Logical - clear out .par, .arp, and other FSC outputs?
                        num_cores = 1,			#Number of processors to use for FSC
                        exec = "fsc26",			#Executable for FSC (needs to be in a folder in the system $PATH)
                        loc_parms = loc_parms2,		#Vector of locus parameters
                        found_Ne = parms$found_Ne,			#Founding population size, required for STEP change model		
                        gmap = gmap,              #Mapping the original population onto aggregated grid
                        MAF = 0.01                #Minor allele frequency threshold, loci with minor allele frequencies below this value are excluded from sim
                        )

#Testing whether the MAF filter works...
if(FALSE) {
  snp.name <- colnames(out[, -(1:2)])
  chrom.names <- regmatches(snp.name, regexpr("^C[[:digit:]]+", snp.name))
  allele_per_loc <- c()
  maf_per_loc <- c()
  for(x in 1:length(chrom.names)) {
   tmpalltab <- table(c(out[,grep(chrom.names[x], colnames(out))[1]], out[,grep(chrom.names[x], colnames(out))[2]]))
   allele_per_loc[x] <- dim(tmpalltab)
   maf_per_loc[x] <- min(tmpalltab)/sum(tmpalltab)
   rm(tmpalltab)
  }
  sum(allele_per_loc == 2)/2  #Number of loci with 2 alleles
  sum(maf_per_loc >= 0.01)/2  #Number of loci with MAF > 0.01
}

popDF <- makePopdf(landscape,"cell")

#var_loc_cols <- which(allele_per_loc == 2) + 2
#FINALout <- out[,c(1,2,sample(var_loc_cols,parms$nloci,replace = FALSE))]
#dim(FINALout)
#sum(apply(FINALout[,-c(1,2)], 2, FUN=function(x) dim(table(x))) == 2)

####### REMAINING TO DO...
#1
#Need to test holoStats() with the new output format from strataG
#Some things will definitely need to change (popDF$row & popDF$col, for example)

#2
#Need a while loop to repeat simulation until we have enough polymorphic markers
#scale the number of loci being simulated based on polymorphism rates in first simulation
#ideally, never run more than 2 fastsimcoal simulations - although these are MUCH faster now...
