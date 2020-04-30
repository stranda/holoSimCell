##
## function to take a landscape and make sure that all sampled sites are in unique cells
##
#' @export
uniqueSampled <- function(landscape)
{
    ok <-  TRUE
    sdf <- landscape$sampdf
    sdf <- sdf[order(sdf$cell),]
    r <- rle(sdf$cell)
    if (max(r$lengths)>1)
    {
        ok <- FALSE
        message(paste("This landscape combines mutiple sampled populations into one raster cell",
                      paste(sdf$abbrev[sdf$cell %in% r$values[r$lengths>1]],collapse=", ")))
    }

    ok
}



##
## function to take output from getpophist_cells and def_grid and test
##  for cell occupancy
##

#' @export
testPophist <- function(ph,landscape)
{
    ok <- T
    popdf <- landscape$sampdf
    nonfilled <- popdf[which(ph$Nvec[popdf$cell,ncol(ph$Nvec)]<1),]
    if (nrow(nonfilled)>0)
    {
        ok <- F
        message(paste("These pops have size < 1 at end of simulation",paste(nonfilled$abbrev,collapse=",")))
        message(paste("Sizes, respectively:", paste(ph$Nvec[nonfilled$cell,ncol(ph$Nvec)],collapse=",")))
    }

    sampsuit <- landscape$hab_suit[,popdf$cell]
    if (min(sampsuit[nrow(landscape$hab_suit),])==0)
    {
        ok <- F
        message(paste("habitat suitability for these sites zero at current time",
                      popdf$cell[which(sampsuit[nrow(landscape$hab_suit),]==0)])
                )
    }
    
    if (!all.equal (extent(landscape$sumrast),extent(landscape$samplocsrast)))
    {
        message("the extent of the landscape and the samples is different")
        ok <- F
    }
    
    ok
}


##
## test the gmap selection against a landscape to see if sampled cells are
##   being combined
##
#' @export
#' 
doesGmapCombine <- function(gmap,landscape)
{
    ok <- FALSE
    sdf <- landscape$sampdf
    tgm <- gmap[sdf$cell,]
    tgm <- tgm[order(tgm$gpop),]
    r <- rle(tgm$gpop)
    if (max(r$lengths)>1)
    {
        ok <- TRUE
        pops <- r$values[which(r$lengths>1)]
        message(paste("this gmap combines populations",paste(sdf$abbrev[sdf$cell%in%gmap$pop[gmap$gpop==pops]],collapse=", ")))
    }
    

    ok
}
