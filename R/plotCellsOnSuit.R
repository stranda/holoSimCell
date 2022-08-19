#' plot the location of cells on landscape
#'
#' Illustration of location of specific cells on the forward simulation landscape.
#'
#' @param cells a vector of cell ids to plot
#' @param landscape the landscape object that drives the simulation (see \code{ashSetupLandscape} for a full description of the landscape object)
#' @param timeslice (default = 1), the layer of the landscape that corresponds to a time step (generation) in the simulation
#'
#' @details
#' This function allows the user to visualize the location of specific cells on the simulation landscape.  We have found this useful for confirming the location of refugial cells in our simulations.
#'
#' @return
#' Returns a plot illustrating the simulated landscape at a particular point in time.  The simulation landscape is shown in light yellow (corresponding to open, or unglaciated, habitat) and cells queried are shown in maroon.
#' 
#' @examples
#' library(holoSimCell)
#' parms <- drawParms(control = system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell"))
#' modchoice <- 1
#' load(file=paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",pollenPulls[[modchoice]]$file))
#' refpops <- pollenPulls[[modchoice]]$refs
#' plotCellsOnSuit(refpops, landscape, timeslice = 1)
#'
#' @seealso \code{\link{ashSetupLandscape}}
#'
#' @export
plotCellsOnSuit <- function(cells,landscape,timeslice=1)
{
    hs=landscape$hab_suit[timeslice,]
    hs=as.numeric(hs>0)
    hs[cells] <- 2
    m <- matrix(hs,nrow=landscape$details$x.dim,ncol=landscape$details$y.dim,byrow=F)
    image(x=c(1:landscape$detail$x.dim),y=c(1:landscape$details$y.dim),z=m,xlab="Columns",ylab="Rows")
}
