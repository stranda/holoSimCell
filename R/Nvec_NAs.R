#' Convert population size in unsuitable cells to NA
#'
#' converts cells in Nvec form to NA if their underlying hab suit is NA for that time point
#'
#' @param ph pophist object
#' @param l landscape containing the rasters
#'
#' @details returns a pophist object where the Nvecs matrix (population size through time) has been altered such that population size in glaciated cells is recorded as NA. All other components of the returned pophist object are unaltered from \code{ph}. Useful for calculation of biotic velocity using \code{enmSdm::bioticVelocity()}.
#'
#' @return a pophist object where population size in unsuitable cells at each time point is recorded as NA
#'
#' @examples
#' library(holoSimCell)
#' parms <- drawParms(control = system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell"))
#' load(file=paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",pollenPulls[[1]]$file))
#' refpops <- pollenPulls[[1]]$refs
#' avgCellsz <- mean(c(res(landscape$sumrast)))
#'
#' ph = getpophist2.cells(h = landscape$details$ncells, xdim = landscape$details$x.dim, ydim = landscape$details$y.dim,
#'                        landscape=landscape,
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
#' times_1k <- seq(-21000,0,by=990)
#'    
#' pharray <- pophistToArray(NvecNAs(ph, landscape), times = times_1G)
#'    
#' metrics <- c('centroid', 'nsQuants', 'summary')
#'
#' BV_permill_all <- bioticVelocity(
#'   x=pharray$pophistAsArray,
#'   times = times_1G,
#'   atTimes = times_1k,
#'   longitude=pharray$longitude,
#'   latitude=pharray$latitude,
#'   metrics=metrics,
#'   quants=c(0.05, 0.1, 0.9, 0.95),
#'   onlyInSharedCells = FALSE)
#'
#' @seealso \code{\link{getpophist2.cells}}, \code{\link{pophistToArray}}, \code{\link[enmSdm]{bioticVelocity}}
#'
#' @export
NvecNAs <- function(ph,l)
{
    ph$Nvecs <- ph$Nvecs[,-1]
    if (is.null(ph)) {stop("must supply a population history")}
    if (is.null(l)) {stop("must supply a landscape")}
    if (dim(ph$Nvecs)[1] != (prod(dim(l$NAmask)[1:2]))) {stop("the number of cells in ph and l are different")}
    if (dim(ph$Nvecs)[2] != (dim(l$NAmask)[3])) {stop("the number of time clicks in ph and l are different")}

    nvo <- ph$Nvecs
    m=raster::as.array(l$NAmask)
    m2=base::array(dim=c(dim(m)[c(2,1,3)]))
    for (tm in 1:dim(m)[3])
    {
        m2[,,tm] <- t(m[,,tm])
        m2[,,tm] <- m2[,dim(m2)[2]:1,tm]
    }
    hs = t(l$hab_suit)
    hs = hs[,c(1:ncol(hs))] #changed from hs = hs[,c(1,1:ncol(hs))]
    nvo[is.na(hs)] <- NA
    ph$Nvecs <- nvo
    ph
}
