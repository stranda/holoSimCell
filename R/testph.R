##
## function to take a landscape and make sure that all sampled sites are in unique cells
##

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



#' Test pophist object before coalescent simulation
#'
#' Function to take output from getpophist_cells and def_grid and test for cell occupancy
#' 
#' @param ph a pophist object, output by getpophist2.cells()
#' @param landscape the landscape object used in the forward simulation
#'
#' @details
#' All cells with population genetic samples must be colonized during the forward simulation.  Simulations that do not fully colonize sampled populations are discarded.
#'
#' @return
#' Returns a logical (T or F) indicating whether all sampled cells are occupied at the end of the forward simulation.
#'
#' @examples
#' library(holoSimCell)
#' parms <- drawParms(control = system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell"))
#' modchoice <- 1
#' load(file=paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",pollenPulls[[modchoice]]$file))
#' refpops <- pollenPulls[[modchoice]]$refs
#' avgCellsz <- mean(c(res(landscape$sumrast)))
#'
#' ph = getpophist2.cells(h = landscape$details$ncells, xdim = landscape$details$x.dim, ydim = landscape$details$y.dim,
#'                        hab_suit=landscape,
#'                        refs=refpops,   
#'                        refsz=parms$ref_Ne,
#'                        lambda=parms$lambda,
#'                        mix=parms$mix,  
#'                        shortscale=parms$shortscale*avgCellsz,  
#'                        shortshape=parms$shortshape, 
#'                        longmean=parms$longmean*avgCellsz,  
#'                        ysz=res(landscape$sumrast)[2], 
#'                        xsz=res(landscape$sumrast)[1], 
#'                        K = parms$Ne) 
#' testPophist(ph, landscape)
#'
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


#' Tests aggregation scheme before coalescent simulation
#'
#' Test the gmap aggregation scheme against a landscape to see if cells with genetic samples are being combined
#'
#' @param gmap data frame that maps forward time populations to genetic populations (from make.gmap())
#' @param landscape the landscape object used in the forward simulation
#'
#' @details
#' Avoid combining cells with two separate genetic samples during cell aggregation
#'
#' @return
#' Returns a logical (TRUE or FALSE) indicating whether the aggregation scheme specified by gmap would combine two cells with genetic samples.
#'
#' @examples
#' library(holoSimCell)
#' parms <- drawParms(control = system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell"))
#' modchoice <- 1
#' load(file=paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",pollenPulls[[modchoice]]$file))
#' refpops <- pollenPulls[[modchoice]]$refs
#' avgCellsz <- mean(c(res(landscape$sumrast)))
#'
#' ph = getpophist2.cells(h = landscape$details$ncells, xdim = landscape$details$x.dim, ydim = landscape$details$y.dim,
#'                        hab_suit=landscape,
#'                        refs=refpops,   
#'                        refsz=parms$ref_Ne,
#'                        lambda=parms$lambda,
#'                        mix=parms$mix,  
#'                        shortscale=parms$shortscale*avgCellsz,  
#'                        shortshape=parms$shortshape, 
#'                        longmean=parms$longmean*avgCellsz,  
#'                        ysz=res(landscape$sumrast)[2], 
#'                        xsz=res(landscape$sumrast)[1], 
#'                        K = parms$Ne) 
#' 
#' gmap=make.gmap(ph$pophist,
#'                xnum=2, #number of cells to aggregate in x-direction
#'                ynum=2) #number of aggregate in the y-direction
#'
#' doesGmapCombine(gmap, landscape)
#'
#' @seealso \code{\link{make.gmap}}, \code{\link{pophist.aggregate}}
#' @export
 
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
