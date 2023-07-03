#' Convert raster suitability to habitat suitability surface for simulations
#'
#' Takes multiple spatial inputs (range shapefiles, enm rasters) and makes a grid for the genetic simulations
#'
#' @param pred is a raster stack with layers of hab suitability (assumes [0-1] with larger numbers better suitability)
#' @param samps locations of actual genetic samples.  SpatialPoints object assumed in same proj as raster input
#' @param raster.proj the proj.4 string to impose on the raster enm (if not already specified) (defaults albers)
#'  
#' @details Takes a raster stack and converts into the internal habitat suitability matrix (rows are time and columns represent cell IDs). This object is returned along with rasters delineating the locations of genetic samples and sum of suitabilities for each cell through time. These are all returned from the function as a new 'landscape' object, which formst he main input into /code{getpophist2.cells}.
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
#' rs <- raster::brick(system.file("extdata/rasters/ccsm_160kmExtent_maxent.tif", package = "holoSimCell"))
#' newrs <- newLandscapeDim(rs,0.45)
#' land <- def_grid_pred2(pred=newrs, samps=transSampLoc(pts, range.epsg=4326, raster.proj=crs(rs)@projargs), raster.proj=crs(rs)@projargs)
#'
#' @seealso \code{\link{ashSetupLandscape}}, \code{\link{getpophist2.cells}}
#'
#' @export
def_grid_pred2 = function(pred=NULL,
                          samps = NULL,
                          raster.proj='+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
                          )
{
    if (is.null(pred))
    {
        stop("need a 'pred' rasterstack (or 3d matrix with corners) in the final resolution")
    }
    
    if (!class(pred) %in% c("RasterBrick"))
    {
            stop("need a 'pred' RasterBrick in the final resolution")           
    } 

    if (is.null(samps))
    {
        stop("need to specify sample points in samps, in a spatial points object in the correct projection")
    }

    nc = dim(pred)[1]
    nr = dim(pred)[2]
    
    ###this step converts rasters to suitabilities slow because raster is slow at reporting values
    habSuit <- do.call(rbind,lapply(1:dim(pred)[3],function(m){
        c((t(as.matrix(raster::subset(pred,m)))[,nc:1]))
        }))


    r <- raster::subset(pred,1)
    values(r) <- T
    samplocs <- mask(r,samps,updateValue=F)
    
    occ_pops = which(colSums(habSuit,na.rm=T)>0)
    empty_pops = which(colSums(habSuit,na.rm=T)==0)

    samprowcol <-  as.data.frame(rowColFromCell(samplocs,cellFromXY(samplocs,as(samps,"Spatial"))))

    samprowcol$row <- dim(samplocs)[1]-samprowcol$row + 1
    samprowcol$abbrev <- as.character(samps$abbrev)
    samprowcol$cell <- (samprowcol$row -1)* dim(samplocs)[2] +  samprowcol$col

    samp_pops = samprowcol$cell

    
    sampstruct = list()
    
    sampstruct[["details"]] = data.frame(ncells = dim(habSuit)[2],
                                         npops = length(occ_pops),
                                         nsamp = length(samp_pops),
                                         nempty= length(empty_pops),
                                         x.dim=ncol(r), y.dim=nrow(r))
    sampstruct[["occupied"]] = occ_pops
    sampstruct[["empty"]] = empty_pops
    sampstruct[["sampled"]] = samp_pops
    sampstruct[["hab_suit"]] = habSuit
    sampstruct[["sumrast"]]=r
    sampstruct[["samplocsrast"]]=samplocs
    sampstruct[["samplocs"]]=samps
    sampstruct[["sampdf"]]=samprowcol
    sampstruct[["NAmask"]]=is.na(pred) #slow (converting raster to values)
    
    sampstruct
}

#' Re-project sample locations into simulation projection
#'
#' transform samppts to raster projection matching coordinate reference system of the simulation landscape.
#'
#' @param samppts a dataframe with 5 cols: pop, abbrev, lat, long, and N
#' @param range.epsg epsg number for the projection in which the points are encoded (default 4326 -> wgs84)
#' @param raster.proj proj4 string for the underlying raster
#'
#' @details Takes a dataframe of sample locations in a projection given by \code{range.epsg} (for example, \code{range.epsg = 4326} corresponds to WGS84; see \url{https://epsg.io}) and reprojects into the projection given in \code{raster.proj}. Depends on the \code{sf} package to make the conversion.
#'
#' @return Feature collection with location of sampled populations converted to the coordinate reference system used in the simulation landscape
#'
#' @examples
#' library(holoSimCell)
#' rs <- raster::brick(system.file("extdata/rasters/ccsm_160kmExtent_maxent.tif", package = "holoSimCell"))
#' newrs <- newLandscapeDim(rs,0.45)
#' samps <- transSampLoc(pts, range.epsg=4326, raster.proj=crs(rs)@projargs)
#' land <- def_grid_pred2(pred=newrs, samps=samps, raster.proj=crs(rs)@projargs)
#'
#' @seealso \code{\link{ashSetupLandscape}}, \code{\link{def_grid_pred2}}, \code{\link{getpophist2.cells}}, \url{https://epsg.io}
#'
#' @export

transSampLoc <- function(samppts=NULL,
                         range.epsg=4326,
                         raster.proj='+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
{
    if (!is.null(samppts))
        {
            samps.sf <- st_as_sf(samppts,coords=c("long","lat"),crs=range.epsg)
            samps.st <- st_transform(samps.sf,crs=raster.proj,type="proj")
        } else { stop("must specify samppts") }
    samps.st
}


#' Resample a raster brick into a new resolution
#'
#' resample the landscape to new dimensions given by continuous scaling number \code{fac}
#'
#' @param rs a raster brick object to resample
#' @param fac is a number that is multiplied times the number of cols to calculate new dims
#'
#' @details Takes a raster brick and resamples into a new resolution more coarse than the current raster brick resolution. \code{fac} is a factor on the interval 0-1 that determines how many cells in the new raster brick, relative to the original raster brick. For example, \code{fac = 0.5} would take the number of cells in each dimension in the original raster brick and reduce the same extent into 50% of the cells. This function depends on enmSdm's /code{rastWithSquareCells} function to do the resampling to a square-celled landscape of a particular resolution. IOW the output raster has square cells. Or at least close to square cells.
#'
#' @return Raster brick with resolution specified by the combination of the input landscape (\code{rs}) and the rescaling parameter (\code{fac})
#'
#' @examples
#' library(holoSimCell)
#' rs <- raster::brick(system.file("extdata/rasters/ccsm_160kmExtent_maxent.tif", package = "holoSimCell"))
#' newrs <- newLandscapeDim(rs,0.45)
#' land <- def_grid_pred2(pred=newrs, samps=transSampLoc(pts, range.epsg=4326, raster.proj=crs(rs)@projargs), raster.proj=crs(rs)@projargs)
#'
#' @seealso \code{\link{ashSetupLandscape}}, \code{\link{def_grid_pred2}}, \code{\link{getpophist2.cells}}
#'
#' @export


newLandscapeDim <- function(rs, fac=1.0)
{
    asp = dim(rs)[1]/dim(rs)[2]

    ncol.new=ceiling(fac * dim(rs)[2])

    nrow.new=ceiling(asp * ncol.new)

    ncells = nrow.new * ncol.new
    
    r = enmSdm::rastWithSquareCells(raster::subset(rs,1),ncells)
    
    raster::resample(rs,raster::brick(r,dim(rs)[3]))
    
}
