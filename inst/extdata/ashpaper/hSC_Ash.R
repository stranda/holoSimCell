#library(holoSimCell)
#devtools::load_all()
### imputed, popmap (individualID->pop mapping), pts (sample locations) and ashpred
### are now built into holoSimCell
### as built in dataframes (in data/ directory)
## check out the extra two command-line args: simdir and outdir
args <- commandArgs(TRUE)
i <- as.numeric(args[1])
nreps <- as.numeric(args[2]) 
who <- as.character(args[3])  
label <- as.character(args[4])
#simdir <- as.character(args[5])
outdir <- as.character(args[5])
if (FALSE) #TRUE means production
    {
        simdir <- system("echo $SCRATCH", intern = TRUE)
    } else {
        simdir <- "."
    }
if(length(args) == 0) {
  i <- 1
  nreps <- 2
  who <- "JDR"
  label <- "June_test"
  outdir <- "~/Desktop/hSC_testing/outdir"
  #simdir <- "~/Desktop/hSC_testing/simdir"
}

if ((is.na(simdir))|(is.na(outdir))) {stop("need to specify a correct simdir and/or outdir")}

cat(paste("Run Details:\ni=",i,"nreps =",nreps,"who=",who,"\nlabel=",label,"\nsimdir=",simdir,"\noutdir=",outdir,"\n"))

library(holoSimCell)

#Set the filename, simulation, and output directories for the run
fn <- paste0(label,"_",i,"_", who, ".csv")


###
### popmap is a built-in dataframe with all empirical individuals and their population membership
###
rownames(popmap) <- popmap[,1]
table(popmap[gsub("fp","",names(imputed)),2]) #printout popnames and samples


###
### There are some cells that contain two empirical populations.  Right now we are dropping one in
### each of these cells with the following code.  'imputed' is a built-in dataframe in the holosimcell
### package that has snp data for fraxinus pennsylvanica.
###
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

rs <- brick(paste0(system.file("extdata","rasters",package="holoSimCell"),"/","study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif"))

####
#### this function (newLandscapeDim) takes a rasterbrick and a proportion of cols to resample to.  So if there are 100 cols
#### and proportion is 0.5, the landscape is resampled to 50 columns (cells get twice as wide).  The rows are resampled to
#### make the cells as square as possible
####  it's easy to change this proportion and see the implications for execution times, etc.
#### on the current (Ash) problem, 0.45 gives a forward time simulation of about 2 minutes.

newrs <- newLandscapeDim(rs,0.45)


print("original rasterbrick dimensions:")
print(dim(rs))

print("new resampled rasterbrick dimensions:")
print(dim(newrs))

icelakesland <- def_grid_pred2(pred=1-newrs,
                                     samps=transSampLoc(samppts,
                                                          range.epsg=4326,
                                                          raster.proj=crs(rs)@projargs),
                                     raster.proj=crs(rs)@projargs
                                     )


landscape <- icelakesland
landscape$hab_suit[landscape$hab_suit > 0] <- 1 #Cells under the glacier have 0 suitability, not NA suitability


if (!uniqueSampled(landscape))
{
  stop("The landscape you are using is combining multiple sampled populations into a single raster cell")
}

###seed is based on time in seconds and the number of characters in the library path
###
###
repl <- 1
while(repl <= nreps) {
  sec=as.numeric(Sys.time())-1500000000
  lp= as.numeric(as.character(nchar(paste(.libPaths(), collapse = " "))))
  slp <- as.integer(floor(sec*lp))
  
  set.seed(as.integer(sec))
  
  if(file.exists("Ash_priors.csv")) {
      parms <- drawParms(control = "Ash_priors.csv")
      print("Choosing local prior param file: Ash_priors.csv")
  } else if (file.exists(system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell")))
  {
      parms <- drawParms(control = system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell"))
      print("Choosing distributed prior param file: extdata/ashpaper/Ash_priors.csv")
  } else 
  {
      parms <- drawParms(control = system.file("extdata/csv","priors.csv",package="holoSimCell"))
      print("Choosing distributed prior param file: extdata/csv/priors.csv")
  }
  

  parms$refs <- sample(c("PA","TX","GA","ALL"), 1, replace = FALSE)
  
  if(parms$refs == "PA") {
    refpops <- c(1000, 1001, 999, 1051, 949)
    #refpops <- 1000
  } else if(parms$refs == "TX") {
    refpops <- c(526, 527, 525, 577, 475)
    #refpops <- 526
  } else if(parms$refs == "GA") {
    refpops <- c(536, 537, 535, 587, 485)
    #refpops <- 536
  } else if(parms$refs == "ALL") {
    refpops <- c(536, 537, 535, 587, 485,
                 526, 527, 525, 577, 475,
                 1000, 1001, 999, 1051, 949)
    #refpops <- c(526,536,1000)
  }

### NEW COMMENT IN SEPT 2020
### Because the size of the cells is now carried through the simulation better, it seems like the
### dispersal parameters should scale with it as well.  The next line calcs the average of the width
### and height of the cells.  This is then multiplied times the params (shortscale and longmean) pulled from
### priors during the call to getpophist2.cells.  So in this parameterization we are still estimating in
### proportion of cell width, but under the hood, all the dimensions are correct.

### Would be easy to change the prior file instead....

avgCellsz <- mean(c(res(landscape$sumrast)))  ### if the cells are not square, then use the average of the cell width and length
                                          #Forward simulation
ph = getpophist2.cells(h = landscape$details$ncells, xdim = landscape$details$x.dim, ydim = landscape$details$y.dim,
                       hab_suit=landscape,
                       refs=refpops,  #set at cell 540 right now 
                       refsz=parms$ref_Ne,
                       mix=parms$mix,  #note how small.
                       shortscale=parms$shortscale*avgCellsz,  # scale parameter of weibull with shape below
                       shortshape=parms$shortshape, #weibull shape
                       longmean=parms$longmean*avgCellsz,  # mean of normal with sd = longmean
                       ysz=res(landscape$sumrast)[2], #height of cell in raster (same units as longmean and shortscale)
                       xsz=res(landscape$sumrast)[1], #width of cell in raster
                       K = parms$Ne) #maximum population size in a grid cell, scaled with hab_suit from landscape object



if (!testPophist(ph,landscape))
{
    print("here is where we could do something about non-colonized sample pops")
}


print ("Creating aggregation map")
###Cell aggregation for coalescent
gmap=make.gmap(ph$pophist,
               xnum=2, #number of cells to aggregate in x-direction
               ynum=2) #number of aggregate in the y-direction

if (doesGmapCombine(gmap,landscape))
{
    stop("Need to look at the resolutions because this gmap combines sampled populations")
}

print("Running the aggregation")

ph2 <- pophist.aggregate(ph,gmap=gmap)

  
  #Parameters specific to coalescent model
  loc_parms <- data.frame(marker = "snp",
                          nloci = parms$nloci,           
                          seq_length = parms$seq_length,
                          mu = parms$mu)
  
  preLGMparms <- data.frame(preLGM_t = parms$preLGM_t/parms$G,		#Time / GenTime
                            preLGM_Ne = parms$preLGM_Ne,
                            ref_Ne = parms$ref_Ne)
  
  parms_out <- as.data.frame(c(ph$struct[which(!names(ph$struct) %in% c(names(parms), names(ph$struct)[grep("refs", names(ph$struct))]))], parms))
  
  #With smaller K, some populations have very very low N at the end of the simulation
  #In those cases, we need to inflate N a bit for the coalescent simulation
  ph2$Nvecs[ph2$Nvecs[,702] > 0 & ph2$Nvecs[,702] < 1,702] <- 1
  
  #Run the coalescent simulation
  setwd(simdir) 
  out <- NULL
  out <- tryCatch({runFSC_step_agg3(ph = ph2,				#A new pophist object - (pophist, Nvecs, tmat, struct, hab_suit, coalhist)
                          l = landscape, 			#A new landscape object - (details, occupied, empty, sampled, hab_suit, sumrast, samplocsrast, samplocs)
                          sample_n = 14,		#Number of sampled individuals per population
                          preLGMparms = preLGMparms,		#This has parms for the refuge, preLGM size and timing
                          label = paste0(label,"_",repl,".",round(runif(1, min=1, max = 1e6),0)),			#Label for FSC simulation files
                          delete_files = TRUE,	#Logical - clear out .par, .arp, and other FSC outputs?
                          num_cores = 1,			#Number of processors to use for FSC
                          exec = "fsc26",			#Executable for FSC (needs to be in a folder in the system $PATH)
                          loc_parms = loc_parms,		#Vector of locus parameters
                          found_Ne = parms$found_Ne,			#Founding population size, required for STEP change model		
                          gmap = gmap,              #Mapping the original population onto aggregated grid
                          MAF = 0.01,                #Minor allele frequency threshold, loci with minor allele frequencies below this value are excluded from sim
                          maxloc = 50000           #Maximum number of marker loci to attempt in a fastsimcoal simulation
  )}, error = function(err) {
    print(err)
    return(NULL)
  })
  
  #Calculate summary statistics
  if(!is.null(out)) {
    popDF <- makePopdf(landscape,"cell")
    stats <- holoStats(out, popDF, cores = 1)
    
    ####Calculate several measures of biotic velocity####
    times_1k <- seq(-21000,0,by=990)
    times_1G <- seq(-21000,0,by=30)
    
    pharray <- pophistToArray(NvecNAs(ph, landscape), times = times_1G)
    
    # metrics to use… I just added “sum”, which will give you total population size
    metrics <- c('centroid', 'nsQuants', 'mean', 'prevalence', 'sum')
    
    # note that this will calculate BV across 990-yr intervals for all four metrics listed at once
    # should be run twice, once with onlyInSharedCells TRUE and once FALSE
    BV_pergen_shared <- bioticVelocity(
      x=pharray$pophistAsArray,
      times = times_1G,
      atTimes = times_1G,
      longitude=pharray$longitude,
      latitude=pharray$latitude,
      metrics=metrics,
      quants=c(0.05, 0.1, 0.9, 0.95),
      onlyInSharedCells = TRUE)
    
    #now with onlyInSharedCells FALSE
    BV_pergen_all <- bioticVelocity(
      x=pharray$pophistAsArray,
      times = times_1G,
      atTimes = times_1G,
      longitude=pharray$longitude,
      latitude=pharray$latitude,
      metrics=metrics,
      quants=c(0.05, 0.1, 0.9, 0.95),
      onlyInSharedCells = FALSE)
    
    #now per mill with onlyInSharedCells TRUE
    BV_permill_shared <- bioticVelocity(
      x=pharray$pophistAsArray,
      times = times_1G,
      atTimes = times_1k,
      longitude=pharray$longitude,
      latitude=pharray$latitude,
      metrics=metrics,
      quants=c(0.05, 0.1, 0.9, 0.95),
      onlyInSharedCells = TRUE)
    
    #now per mill with onlyInSharedCells FALSE
    BV_permill_all <- bioticVelocity(
      x=pharray$pophistAsArray,
      times = times_1G,
      atTimes = times_1k,
      longitude=pharray$longitude,
      latitude=pharray$latitude,
      metrics=metrics,
      quants=c(0.05, 0.1, 0.9, 0.95),
      onlyInSharedCells = FALSE)
    
    #Then calculate the metrics you want to record from these 4 objects and add names
    #Maximum per-generation biotic velocity attained during the simulation
    BVmaxSHARED <- max(BV_pergen_shared$centroidVelocity)
    BVmaxALL <- max(BV_pergen_all$centroidVelocity)
    #prevalence per millenium
    BVprevMILL <- BV_permill_all$prevalence
    names(BVprevMILL) <- paste0("BVprev_",times_1k[-1])
    #mean abundance per millenium
    BVmeanMILL <- BV_permill_all$mean
    names(BVmeanMILL) <- paste0("BVmean_",times_1k[-1])
    #total abundance per millenium
    BVtotMILL <- BV_permill_all$sum
    names(BVtotMILL) <- paste0("BVtotal_",times_1k[-1])
    #biotic velocity per millenium
    BVsharedMILL <- BV_permill_shared$centroidVelocity
    names(BVsharedMILL) <- paste0("BVshared1k_",times_1k[-length(times_1k)],"to",times_1k[-1])
    BVallMILL <- BV_permill_all$centroidVelocity
    names(BVallMILL) <- paste0("BVall1k_",times_1k[-length(times_1k)],"to",times_1k[-1])
    #quantiles of per-generation biotic velocity
    BVsharedQUANT <- quantile(BV_pergen_shared$centroidVelocity)
    names(BVsharedQUANT) <- paste0("BVshared1G_",names(BVsharedQUANT))
    BVallQUANT <- quantile(BV_pergen_all$centroidVelocity)
    names(BVallQUANT) <- paste0("BVall1G_",names(BVallQUANT))
    #North quantile movement per millenium
    BVNQsharedMILL <- BV_permill_shared$nsQuantVelocity_quant0p95
    names(BVNQsharedMILL) <- paste0("BVNQshared1k_",times_1k[-length(times_1k)],"to",times_1k[-1])
    BVNQallMILL <- BV_permill_all$nsQuantVelocity_quant0p95
    names(BVNQallMILL) <- paste0("BVNQall1k_",times_1k[-length(times_1k)],"to",times_1k[-1])
    #South quantile movement per millenium
    BVSQsharedMILL <- BV_permill_shared$nsQuantVelocity_quant0p05
    names(BVSQsharedMILL) <- paste0("BVSQshared1k_",times_1k[-length(times_1k)],"to",times_1k[-1])
    BVSQallMILL <- BV_permill_all$nsQuantVelocity_quant0p05
    names(BVSQallMILL) <- paste0("BVSQall1k_",times_1k[-length(times_1k)],"to",times_1k[-1])
    #Quantiles of north quantile movement per generation
    BVNQsharedQUANT <- quantile(BV_pergen_shared$nsQuantVelocity_quant0p95)
    names(BVNQsharedQUANT) <- paste0("BVNQshared1G_", names(BVNQsharedQUANT))
    BVNQallQUANT <- quantile(BV_pergen_all$nsQuantVelocity_quant0p95)
    names(BVNQallQUANT) <- paste0("BVNQall1G_", names(BVNQallQUANT))
    #Quantiles of south quantile movement per generation
    BVSQsharedQUANT <- quantile(BV_pergen_shared$nsQuantVelocity_quant0p05)
    names(BVSQsharedQUANT) <- paste0("BVSQshared1G_", names(BVSQsharedQUANT))
    BVSQallQUANT <- quantile(BV_pergen_all$nsQuantVelocity_quant0p05)
    names(BVSQallQUANT) <- paste0("BVSQall1G_", names(BVSQallQUANT))
    
    #Combine parameters and sumstats into one vector
    all_out <- c(date = date(), node=i, rep=repl, parms_out, 
                 BVmaxALL=BVmaxALL, BVmaxSHARED = BVmaxSHARED, BVprevMILL, BVmeanMILL, BVtotMILL,
                 BVsharedMILL, BVsharedQUANT, BVNQsharedMILL, BVSQsharedMILL, BVNQsharedQUANT, BVSQsharedQUANT,
                 BVallMILL, BVallQUANT, BVNQallMILL, BVSQallMILL, BVNQallQUANT, BVSQallQUANT, stats)
    
    #Write output
    setwd(outdir)
    if(!file.exists(fn)) {
      write.table(all_out, fn, sep = ",", quote = FALSE, row.names = FALSE)
    } else {
      write.table(all_out, fn, quote = FALSE, row.names = FALSE, sep=",", append = TRUE, col.names = FALSE)
    }
    
    rm(all_out)
    
    repl <- repl+1
  }
  
}
