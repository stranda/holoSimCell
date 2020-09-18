#Function to impose a grid on a shapefile (as in species range)

#Shapefiles used to be downloaded from: https://gec.cr.usgs.gov/data/little/
#But that link is now dead

#It appears that they can be downloaded from: https://www.fs.fed.us/nrs/atlas/littlefia/#
#But not from this site (all dead links): http://www.atozmapsdata2.com/downloads/Country/Modern/C-USA-USGSTreeSpecies-index.html
#There is also a GitHub repository with archived data from the top site above: https://github.com/wpetry/USTreeAtlas

#GitHub is likely the preferred solution -- https://github.com/wpetry/USTreeAtlas

#There may be some CRS issues in here.  Where can I find the CRS in a shapefile???  
                                        #Doesn't appear to be in the fraxpenn folder
#' Takes multiple spatial inputs (range shapefiles, enm rasters) and makes a grid for the genetic simulations
#' @param pred is a raster stack with layers of hab suitability (assumes [0-1] with larger numbers better suitability)
#' @param samps locations of actual genetic samples.  SpatialPoints object assumed in same proj as raster input
#' @param raster.proj the proj.4 string to impose on the raster enm (if not already specified) (defaults albers)
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

#' transform samppts to raster projection.
#' @param samppts a dataframe with 5 cols: pop, abbrev, lat, long, and N
#' @param range.epsg epsg number for the projection in which the points are encoded (default 4326 -> wgs84)
#' @param raster.proj proj4 string for the underlying raster
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


#### resample to lower resolution
#### Adam may have a function for this
#' resample the landscape to new dimensions given by continuous scaling number 'fac'
#' @param rs a raster brick object to resample
#' @param fac is a number that is multiplied times the number of cols to calculate new dims
#' @description This function depends on enmSdm's 'rastWithSquareCells()' function to do the resampling to
#' a square-celled landscape of a particular resolution.  IOW the output raster has square cells.  Or at least close to square cells.
#' @export
newLandscapeDim <- function(rs, fac=1.0)
{
    asp = dim(rs)[1]/dim(rs)[2]

    ncol.new=ceiling(fac * dim(rs)[2])

    nrow.new=ceiling(asp * ncol.new)

    ncells = nrow.new * ncol.new
    
    r = enmSdm::rastWithSquareCells(raster::subset(rs,1),ncells)
    
    resample(rs,raster::brick(r,dim(rs)[3]))
    
}
