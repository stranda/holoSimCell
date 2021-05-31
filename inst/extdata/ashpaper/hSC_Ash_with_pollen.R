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
if (TRUE) #TRUE means production
    {
        simdir <- system("echo $SCRATCH", intern = TRUE)
    } else {
        simdir <- "."
    }
if(length(args) == 0) {
  i <- 1
  nreps <- 2
  who <- "JDR"
  label <- "Ashpaper_test"
  outdir <- "~/Desktop/hSC_testing/outdir"
  #simdir <- "~/Desktop/hSC_testing/simdir"
}

if ((is.na(simdir))|(is.na(outdir))) {stop("need to specify a correct simdir and/or outdir")}

cat(paste("Run Details:\nPollen surfaces\ni=",i,"nreps =",nreps,"who=",who,"\nlabel=",label,"\nsimdir=",simdir,"\noutdir=",outdir,"\n"))

tmp = .libPaths()
.libPaths(tmp[!grepl("home",tmp)])  ##remove any paths that have "home" in them, like: /home/f0007250  

library(holoSimCell)

#Set the filename, simulation, and output directories for the run
fn <- paste0(label,"_",i,"_", who, ".csv")

###
###All of the genetic subsampling and the landscape processing is now abstracted into this function (and one it calls)
###This file replicates the old setup that we have used, but there needs to be a new setup for the multiple hab_suits
### I think we should run a separate routine at package creation or installation that converts all of the rasters to
### landscapes and pushes them into a list that we can read off the disk once per simulation rather than repeatedly reading them
### I'm assuming the memory cost is not too high

lnum=nrow(pollenSurfaces$pulls) #number of landscapes to make (number of enm rasters)


if (FALSE)  #logic used to create landscapes from pollen objects--don't run, built-in to package for speed
{
    ##get the suitabilities
    tmp <- readRDS("~/GoogleDrive/doc/proposals/nsf/2017/NSF_ABI_2018_2021/data_and_analyses/pg_pollen/preds_for_ABC_n50.RDS")
    for (m in 1:length(tmp) )
    {
        print(m)
        landscape <- ashSetupLandscape(brickname=raster::brick(tmp[[m]]),cellreduce=0.45,partialsuit=T)
        save(file=paste0("../landscapes/pollenPull_",m,".rda"),landscape)
    }
    ##get the refuge numbers
    tmp <- readRDS("~/GoogleDrive/doc/proposals/nsf/2017/NSF_ABI_2018_2021/data_and_analyses/pg_pollen/refuge_rasters_n50.RDS")
    

    
}


 
###seed is based on time in seconds and the number of characters in the library path
###
###
repl <- 1
while(repl <= nreps) {
  sec=as.numeric(Sys.time())-1500000000
  slp <- sec+(as.numeric(i))
  
  set.seed(as.integer(slp))
  
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
  

###choose enm model and refugia
###
# modchoice <- as.integer(round(runif(1,min=1,max=nrow(enmScenarios$enms))))
  modchoice <- sample(c(1:nrow(enmScenarios$enms)),1,replace=FALSE)
  
  parms$refs <- paste0("ENM_",modchoice)

###read the pre-calculated landscape off the disk.  stored in inst/extdata/landscapes/*.rda
  load(file=paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",enmScenarios$enms$rasterStackName[modchoice],".rda"))
  
  refpops <- enmScenarios$refugeCellNum[[modchoice]] 

#  rast_grid <-
#      matrix(data = c(1:landscape$details$ncells), nrow = landscape$details$y.dim, ncol = landscape$details$x.dim, byrow = TRUE)
#  hSC_grid <- rast_grid[nrow(rast_grid):1,]
#  refpops <- hSC_grid[which(rast_grid %in% enmScenarios$refugeCellIds[[modchoice]])]

  
  ### NEW COMMENT IN SEPT 2020
  ### Because the size of the cells is now carried through the simulation better, it seems like the
  ### dispersal parameters should scale with it as well.  The next line calcs the average of the width
  ### and height of the cells.  This is then multiplied times the params (shortscale and longmean) pulled from
  ### priors during the call to getpophist2.cells.  So in this parameterization we are still estimating in
  ### proportion of cell width, but under the hood, all the dimensions are correct.

  ### Would be easy to change the prior file instead....

  avgCellsz <- mean(c(res(landscape$sumrast)))  ### if the cells are not square, then use the average of the cell width and length
                                        #Forward simulation
  save(file="parms.rda",parms,landscape,refpops,avgCellsz)
  ph = getpophist2.cells(h = landscape$details$ncells, xdim = landscape$details$x.dim, ydim = landscape$details$y.dim,
                       hab_suit=landscape,
                       refs=refpops,  #set at cell 540 right now 
                       refsz=parms$ref_Ne,
                       lambda=parms$lambda,
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
    incomplete <- TRUE
  } else {
    incomplete <- FALSE
  }

  if(incomplete == FALSE) {
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
      metrics <- c('centroid', 'nsQuants', 'summary')
    
      # note that this will calculate BV across 990-yr intervals for all four metrics listed at once
      # should be run twice, once with onlyInSharedCells TRUE and once FALSE
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
    
      #velocity from 21,000 ybp to today
      BV_21kyr_shared <- bioticVelocity(
        x=pharray$pophistAsArray,
        times = times_1G,
        atTimes = c(-21000,0),
        longitude=pharray$longitude,
        latitude=pharray$latitude,
        metrics=metrics,
        quants=c(0.05, 0.95),
        onlyInSharedCells = TRUE)

      #Then calculate the metrics you want to record from these 4 objects and add names
      #prevalence per millenium
      BVprevMILL <- BV_permill_all$prevalence
      names(BVprevMILL) <- paste0("BVprev_",-1*times_1k[-1],"ybp")
      #mean abundance per millenium
      BVmeanMILL <- BV_permill_all$mean
      names(BVmeanMILL) <- paste0("BVmean_",-1*times_1k[-1],"ybp")
      #total abundance per millenium
      BVtotMILL <- BV_permill_all$sum
      names(BVtotMILL) <- paste0("BVtotal_",-1*times_1k[-1],"ybp")
      #biotic velocity per millenium
      BVsharedMILL <- BV_permill_shared$centroidVelocity
      names(BVsharedMILL) <- paste0("BVshared1k_",-1*times_1k[-length(times_1k)],"to",-1*times_1k[-1],"ybp")
      BVallMILL <- BV_permill_all$centroidVelocity
      names(BVallMILL) <- paste0("BVall1k_",-1*times_1k[-length(times_1k)],"to",-1*times_1k[-1],"ybp")
      #North quantile movement per millenium
      BVNQsharedMILL <- BV_permill_shared$nsQuantVelocity_quant0p95
      names(BVNQsharedMILL) <- paste0("BVNQshared1k_",-1*times_1k[-length(times_1k)],"to",-1*times_1k[-1],"ybp")
      BVNQallMILL <- BV_permill_all$nsQuantVelocity_quant0p95
      names(BVNQallMILL) <- paste0("BVNQall1k_",-1*times_1k[-length(times_1k)],"to",-1*times_1k[-1],"ybp")
      #South quantile movement per millenium
      BVSQsharedMILL <- BV_permill_shared$nsQuantVelocity_quant0p05
      names(BVSQsharedMILL) <- paste0("BVSQshared1k_",-1*times_1k[-length(times_1k)],"to",-1*times_1k[-1],"ybp")
      BVSQallMILL <- BV_permill_all$nsQuantVelocity_quant0p05
      names(BVSQallMILL) <- paste0("BVSQall1k_",-1*times_1k[-length(times_1k)],"to",-1*times_1k[-1],"ybp")
      #Velocity over 21kyr - centroid, northern and southern margins
      BV_21kyr <- BV_21kyr_shared$centroidVelocity
      BVNQ_21kyr <- BV_21kyr_shared$nsQuantVelocity_quant0p95
      BVSQ_21kyr <- BV_21kyr_shared$nsQuantVelocity_quant0p05
    
      #Combine parameters and sumstats into one vector
      all_out <- c(date = date(), node=i, rep=repl, parms_out, 
                 BVprevMILL, BVmeanMILL, BVtotMILL,
                 BVsharedMILL, BVNQsharedMILL, BVSQsharedMILL, 
                 BVallMILL, BVNQallMILL, BVSQallMILL, 
                 BV_21kyr = BV_21kyr, BVNQ_21kyr = BVNQ_21kyr, BVSQ_21kyr= BVSQ_21kyr,
                 stats)
    
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
}
