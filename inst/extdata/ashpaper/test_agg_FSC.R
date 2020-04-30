library(holoSimCell)
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

if ((!exists("icenolakesland")))
{
    ##this should produce a landscape with (x,y) _square_ cells that also have
    ##21empirical samples in separate grid cells (otherwise need to figure out something else)
    ## ashland is a stored R object as well
    
    ##read in the correct raster stack
    rs <- stack("study_region_daltonIceMask_noLakes_linearIceSheetInterpolation.tif")
    e <- extent(rs)
    corners <- (matrix(c( e[1], e[4],
                         e[1], e[3],
                         e[2], e[3],
                         e[2],e[4]),ncol=2,byrow=T))
    colnames(corners) <- c("x","y")
    rownames(corners) <- c("ul","ll","lr","ur")
    
    icenolakes <- as.array(1-rs) #really slow! thats one reason these layers get stored
    
    icenolakesland <- def_grid_pred(pred=icenolakes,samppts=samppts,
                                    init.ext=c(45,65),
                                    keep.thresh=0.05,corners=corners)
}

if (!exists("icelakesland"))
{
 rs <- stack("study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif")
    e <- extent(rs)
    corners <- (matrix(c( e[1], e[4],
                         e[1], e[3],
                         e[2], e[3],
                         e[2],e[4]),ncol=2,byrow=T))
    colnames(corners) <- c("x","y")
    rownames(corners) <- c("ul","ll","lr","ur")
    
    icelakes <- as.array(1-rs) #really slow! thats one reason these layers get stored
    
    icelakesland <- def_grid_pred(pred=icelakes,samppts=samppts,
                                  init.ext=c(45,65),
                                  keep.thresh=0.01,corners=corners)
}

landscape <- icelakesland

if (!uniqueSampled(landscape))
{
    stop("The landscape you are using is combining multiple sampled populations into a single raster cell")
}



###seed is based on time in seconds and the number of characters in the library path
###
###
sec=as.numeric(Sys.time())-1500000000
lp= as.numeric(as.character(nchar(paste(.libPaths(), collapse = " "))))
slp <- as.integer(floor(sec*lp))

set.seed(as.integer(sec))

ph = getpophist2.cells(hab_suit=landscape,
                       refs=(540),
                       refsz=100,
                       mix=0.005,  #note how small.
                       shortscale=0.04,  # scale parameter of weibull with shape below
                       shortshape=1, #weibull shape
                       longmean=0.15,  # mean of normal with sd = longmean
                       sz=1) #size of a cell (same units as longmean and shortscale)

if (!testPophist(ph,landscape))
{
    print("here is where we could do something about non-colonized sample pops")
}


gmap=make.gmap(ph$pophist,
               xnum=3, #number of cells to aggregate in x-direction
               ynum=3) #number of aggregate in the y-direction


if (doesGmapCombine(gmap,landscape))
{
    stop("Need to look at the resolutions because this gmap combines sampled populations")
}


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
                        loc_parms = loc_parms,		#Vector of locus parameters
                        found_Ne = parms$found_Ne,			#Founding population size, required for STEP change model		
                        gmap = gmap,              #Mapping the original population onto aggregated grid
                        MAF = 0.01                #Minor allele frequency threshold, loci with minor allele frequencies below this value are excluded from sim
)

popDF <- makePopdf(landscape,"cell")

stats <- holoStats(out, popDF, cores = 1)
