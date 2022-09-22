##' Subset built-in ash data based on populations
##' @param pops a vector of population names

ashRemoveGeneticPops <- function(popmap,pops)
{
    imputed.pruned=imputed

    for (p in pops)
        imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev==p,"id"])]


    removes <- c()
    rownames(popmap) <- popmap$id
    popids <- popmap[gsub("fp","",names(imputed.pruned)),2]
    
    ##table(popids)
    
    for (a in unique(popids))
    {
        if (sum(popids==a)>14)
        {
            removes <- c(removes,sample(which(popids==a),1))
        }
    }
    imputed.pruned[,-1*removes]
}

#' Establish simulation landscape
#'
#' Setup a landscape for Ash for our simulations. everything but the surface is baked into holoSimCell
#'
#' @param brickname name of a file of a geotiff object.  layers correspond to time clicks in simulations
#' @param equalsuit logical (TRUE / FALSE) indicating whether suitability is assumed to be equal in all cells across the landscape
#' @param partialsuit logical (TRUE / FALSE) indicating whether fractional habitat suitabilities should be used 
#' @param cellreduce scalar specifying the proportional reduction in landscape dimensions from the input geotiff file (new dimensions = cellreduce x original dimensions)
#' @param xlim optionally used to define the longitudinal dimension of the output landscape (if no geotiff object is supplied)
#' @param ylim optionally used to define the latitudinal dimension of the output landscape (if no geotiff object is supplied)
#' @param timesteps specifies the temporal dimension (number of time steps) of the output landscape
#'
#' @details
#' This function reads in a geotiff raster brick and converts the individual rasters into matrices of habitat suitability. Returns a landscape object to be used in forward demographic simulations.
#'
#' @return
#' Returns a landscape object with the following components:
#' \itemize{
#' \item{\code{details}} {Data frame with the spatial extent of the landscape grid}
#' \item{\code{occupied}} {Vector of population IDs for occupied cells in the landscape grid}
#' \item{\code{empty}} {Vector of population IDs for empty cells in the landscape grid}
#' \item{\code{sampled}} {Vector of population IDs for sampled populations in the landscape grid}
#' \item{\code{hab_suit}} {Matrix with species habitat suitability (ranging from zero to one) through time. Rows in the matrix correspond to discrete time units and columns correspond to cells in the landscape}
#' \item{\code{sumrast}} {Raster that stores the extent and resolution of the simulation landscape}
#' \item{\code{samplocsrast}} {Raster showing the locations of cells that correspond to sampled populations}
#' \item{\code{samplocs}} {Simple feature encoding spatial vector data related to the sampling populations and a data frame with metadata of the sampling populations (population id, number of individuals, type of spatial data and coordinates)}
#' \item{\code{sampdf}} {Data frame with the spatial location of the sampling populations in the landscape grid}
#' \item{\code{NAmask}} {A RasterBrick object used to mask cells that are unsuitable (e.g., covered by glaciers, out of the study region, in large lakes or oceans)}
#' }
#'
#' @examples
#' library(holoSimCell)
#' m <- 1
#' landscape <- ashSetupLandscape(brickname=system.file("extdata", "rasters", "ccsm_160kmExtent_maxent.tif", package = "holoSimCell"),cellreduce=0.45,partialsuit=T)
#'
#' @seealso \code{\link{getpophist2.cells}}, \code{\link{def_grid_pred2}}
#'
#' @export
ashSetupLandscape <- function(brickname=system.file("extdata", "rasters", "ccsm_160kmExtent_maxent.tif", package = "holoSimCell"),equalsuit=F,partialsuit=F,cellreduce=0.45,xlim=NULL,ylim=NULL,timesteps=NULL)
{

    rownames(popmap) <- popmap[,1]
###
### There are some cells that contain two empirical populations.  Right now we are dropping one in
### each of these cells with the following code.  'imputed' is a built-in dataframe in the holosimcell
### package that has snp data for fraxinus pennsylvanica.
###
    imputed.pruned <- ashRemoveGeneticPops(popmap=popmap,pops=c("Michigan","UNK","MO1","ON1","VA1","MB1")) #these are the pops that are removed
    poptbl <- table(popmap[gsub("fp","",names(imputed.pruned)),2])

    samppts <- pts[pts$abbrev %in% names(poptbl),]
    


####
#### this function (newLandscapeDim) takes a rasterbrick and a proportion of cols to resample to.  So if there are 100 cols
#### and proportion is 0.5, the landscape is resampled to 50 columns (cells get twice as wide).  The rows are resampled to
#### make the cells as square as possible
####  it's easy to change this proportion and see the implications for execution times, etc.
#### on the current (Ash) problem, 0.45 gives a forward time simulation of about 2 minutes.

    if (!equalsuit)
    {
        if ("RasterBrick" %in% class(brickname))
        {
            rs <- brickname #passed a brick instead of a path of a brick in a file

        } else {
            rs <- raster::brick(brickname) # passed a brick filename
        }
        
        newrs <- newLandscapeDim(rs,cellreduce)
        
    } else {
        if (is.null(xlim)|is.null(ylim))
            newrs=raster::raster(nrows=50,ncols=50,xmn=0,xmx=5000,ymn=0,ymx=5000,vals=1)
        else
            newrs=raster::raster(nrows=ylim,ncols=xlim,xmn=0,xmx=5000,ymn=0,ymx=5000,vals=1)
        newrs = raster::brick(newrs,nl=701)
    }


    
land <- def_grid_pred2(pred=newrs,
                                     samps=transSampLoc(samppts,
                                                          range.epsg=4326,
                                                          raster.proj=crs(rs)@projargs),
                                     raster.proj=crs(rs)@projargs
                                     )


landscape <- land
if (partialsuit==F) landscape$hab_suit[landscape$hab_suit > 0] <- 1 #Cells under the glacier have 0 suitability, not NA suitability
if (!uniqueSampled(landscape))
{
  stop("The landscape you are using is combining multiple sampled populations into a single raster cell")
}

landscape

}
