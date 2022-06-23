#' Create a popDF object
#'
#' Generates a population data frame with information on site IDs and cell or geographic coordinates of sampled populations
#'
#' @param l the landscape object used in the forward simulation
#' @param system specifies the coordinate system to use, either cells (x and y positions on the landscape) or geographic coordinates (latitude and longitude)
#'
#' @details
#' The data frame produced by this function is used as an input for \code{holoStats}, particularly for spatial summaries of genetic variation on the landscape (e.g., isolation by distance analysis, relationships between genetic diversity and latitude / longitude)
#'
#' @return
#' Returns a three-column data frame with population IDs (pop) and column (x) / row (y) locations on the simulated landscape for each sampled population
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
#' ph2 <- pophist.aggregate(ph,gmap=gmap)
#'
#' loc_parms <- data.frame(marker = "snp",
#'                         nloci = parms$nloci,           
#'                         seq_length = parms$seq_length,
#'                         mu = parms$mu)
#'   
#' preLGMparms <- data.frame(preLGM_t = parms$preLGM_t/parms$G,   
#'                           preLGM_Ne = parms$preLGM_Ne,
#'                          ref_Ne = parms$ref_Ne)
#' 
#' out <- runFSC_step_agg3(ph = ph2,
#'                         l = landscape,
#'                         sample_n = 14,
#'                         preLGMparms = preLGMparms,
#'                         label = "test",
#'                         delete_files = TRUE,
#'                         num_cores = 1,
#'                         exec = "fsc26",
#'                         loc_parms = loc_parms,
#'                         found_Ne = parms$found_Ne,
#'                         gmap = gmap,
#'                         MAF = 0.01,
#'                         maxloc = 50000)
#' popDF <- makePopdf(landscape,"cell")
#'
#' @seealso \code{\link{holoStats}}
#' @export


makePopdf <- function(l,system=c("geo","cell")[1])
{
    if (!is.null(l))
    {
        if (system=="geo")
        {
            sl <- l$samplocs
            crd <- t(sapply(sl$geometry,function(x){c(x=x[[1]],y=x[[2]])}))
            ret <- data.frame(pop=sl$abbrev,x=crd[,1],y=crd[,2])
        } else if (system=="cell")
        {
            sdf <- l$sampdf
            ret <- data.frame(pop=sdf$abbrev,x=sdf$col,y=sdf$row)
        } else stop("'system' has to be either 'geo' or 'cell'")
    } else stop("input landscape does not exist")
    ret[order(ret$pop),]
}
