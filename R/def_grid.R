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
#' @param pred is a 3d matrix with layers of hab suitability indexed by the third dimension
#' @param samppts locations of actual genetic samples.  Dataframe and must have 'long' and 'lat' columns in wgs84
#' @param init.ext vector with the number of columns and rows hoping for in output (usually is altered in run)
#' @param range.epsg the epsg number to impose on the shapefile containing range (if not already specified) defaults wgs84)
#' @param raster.proj the proj.4 string to impose on the raster enm (if not already specified) (defaults albers)
#' @param sim.epsg the epsg number that the sim runs in and upon which the resulting grid is calculated
#' @param keep.thresh the proportion of a cell that has to be occupied before it's included in simulation model
#' 
#' @export

def_grid_pred= function(pred=NULL, samppts = NULL,
                        init.ext = NULL,
                        raster.proj='+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs',
                        keep.thresh = 0.1,
                        sim.epsg=5070,
                        range.epsg=4326,
                        corners = matrix(c( -1702657, 3970378,
                                           -1702657, -100000,
                                           2908121, -100000,
                                           2908212,3970378),ncol=2,byrow=T))
{

    if (is.null(pred))
    {
        stop("need a 'pred', (3d object with x,y equal space and z equals time)")
    }
    
    if (!is.null(samppts))
    {
        samps <- st_as_sf(samppts,coords=c("long","lat"),crs=range.epsg)
        samps <- st_transform(samps,sim.epsg)
    }
    
    sumextent <- apply(pred,c(1,2),sum,na.rm=T)
    pred <- pred[rowSums(sumextent)!=0,colSums(sumextent)!=0,]
    r <- setRasterExtent(pred[dim(pred)[1]:1,,1],corners=corners)
    
    e <- projectExtent(r,raster.proj)    
    r <- projectRaster(r,e)
    
    rast.extent <- extent(r)
    rste <-sapply(rast.extent,c) 
    
    if (!is.null(samppts)) {
        samps.extent <- extent(samps)
        se <- sapply(samps.extent,c)
    } else se=NA
    
###using the raster as the base, need to spend more time on logic above to use all sources. or downsample
    
    if (!is.null(init.ext)) ### there is a size to resample to.  Will do for every 3rd dim layer
    {

        ##try to get new cells square
        r <- setRasterExtent(pred[dim(pred)[1]:1,,1],corners=corners)
        asp=aspRaster(r)
        er <- c("x"=abs(diff(extent(r)[1:2])),"y"=abs(diff(extent(r)[3:4])))
        cx=init.ext[1]
        cy=init.ext[2]
        dx=er["x"]/cx
        dy=er["y"]/cy
        if (asp>1)
        {
            cx=round(er["x"]/dy)
        } else
            if (asp<1)
            {
                cy=round(er["y"]/dx)
            }
        out.ext <- c(cx,cy)

        newr <- raster(ncol=out.ext[1],nrow=out.ext[2],crs=raster.proj,
                       xmn=rste[1],xmx=rste[2],ymn=rste[3],ymx=rste[4])
        newpred <- array(NA,dim=c(out.ext[2],out.ext[1],dim(pred)[3]))
        for (i in 1:dim(pred)[3])
        {
            p <- pred[,,i]
            r <- setRasterExtent(p[dim(p)[1]:1,],corners=corners)
            rr <- resample(r,newr)
            newpred[,,i] <- as.matrix(rr)
        }
        pred <- newpred[dim(newpred)[1]:1,,]
    } else {
        pred <- pred[dim(pred)[1]:1,,]
    }

    ###this next bit of code subsets cells based on the keep.threshold
    sumextent <- apply(pred,c(1,2),sum,na.rm=T)
    sumextent <- sumextent/max(sumextent)
    inds <- which(sumextent<=keep.thresh,arr.ind=T) #these are the inds to leave out (if keep.thresh>=0, then zeros will be removed)
    for (rw in 1:nrow(inds))
    {
        row=inds[rw,1]
        col=inds[rw,2]
        pred[row,col,] <- NA 
    }

    r <- setRasterExtent(sumextent,corners=corners)
    samplocs <- mask(r,samps)
    
    habSuit <- do.call(rbind,lapply(1:dim(pred)[3],function(m){
        c(t(pred[,,m]))
        }))
    
    
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
    
    sampstruct
}


setRasterExtent <- function(mat,corners=NULL)
{
    if (is.null(corners))
    {
        corners <- (matrix(c( -1702657, 3970378,
                             -1702657, -100000,
                             2908121, -100000,
                             2908212,3970378),ncol=2,byrow=T))
        colnames(corners) <- c("x","y")
        rownames(corners) <- c("ul","ll","lr","ur")
    }
    if (is.null(colnames(corners))) colnames(corners) <- c("x","y")
    if (is.null(rownames(corners))) rownames(corners) <- c("ul","ll","lr","ur")
    xr <- abs(corners["ul","x"]-corners["ur","x"])/dim(mat)[2]
    yr <- abs(corners["ul","y"]-corners["ll","y"])/dim(mat)[1]
    r=raster(mat[dim(mat)[1]:1,],xmn=corners["ul","x"],xmx=corners["ur","x"],
             ymx=corners["ul","y"],ymn=corners["ll","y"],
             crs='+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    r
}

trimExtent <- function(e,cols=c(1:4,55:66),rows=c(1:3,38:39))
{
    e[-rows,-cols,]
}

###both of these funtions take the enm matrices and make sure they
###convert to our cell/population numbering system
cell2coords <- function(x,nrow,ncol)
{
   c(( x %/% ncol )+1,x %% ncol)
}

coords2cell <- function(r,c,ncol)
{
    (r-1) * ncol +  c
}

#' gets the aspect ratio of a raster
#' @param r a raster object
#' @export
aspRaster <- function(r)
{
 res(r)[2]/res(r)[1]
}
