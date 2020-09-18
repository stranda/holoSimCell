###
### testing out the landscape creation and forward-time simulation
###
library(holoSimCell)



###make sure the sampled populations are correct
###imputed is a built-in dataset based on fraxinus pennsylvanica
rownames(popmap) <- popmap[,1]
table(popmap[gsub("fp","",names(imputed)),2])

###deleting populations to ensure a single pop per cell.  need to fix this and use the full empirical data
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

newrs <- enmSdm:::squareRasterCells(newLandscapeDim(rs,0.45))

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
landscape$hab_suit[!is.na(landscape$hab_suit)] <- 1   #Naive habitat suitability, habitable cells all have suitability of 1

if (!uniqueSampled(landscape))
{
  stop("The landscape you are using is combining multiple sampled populations into a single raster cell")
}

repl <- 1

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
    refpops <- c(1301, 1302, 1300, 1345, 1257)
                                        #refpops <- 1301
} else if(parms$refs == "TX") {
    refpops <- c(586, 587, 585, 630, 542)
                                        #refpops <- 586
} else if(parms$refs == "GA") {
    refpops <- c(638, 639, 637, 682, 594)
                                        #refpops <- 638
} else if(parms$refs == "ALL") {
    refpops <- c(1301, 1302, 1300, 1345, 1257,
                 586, 587, 585, 630, 542,
                 638, 639, 637, 682, 594)
                                        #refpops <- c(1301, 587, 638)
}


### Because the size of the cells is now carried through the simulation better, it seems like the
### dispersal parameters should scale with it as well.  The next line calcs the average of the width
### and height of the cells.  This is then multiplied times the params pulled from
### priors during the call to getpophist2.cells

avgCellsz <- mean(c(res(landscape$sumrast)))

print(parms)
                                        #Forward simulation
ph = getpophist2.cells(hab_suit=landscape,
                       refs=refpops,  #set at cell 540 right now 
                       refsz=parms$found_Ne,
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

