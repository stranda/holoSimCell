#' Shiny interface for population model
#'
#' Interactive interface to visualize simulation histories under different parameterizations
#' 
#' @details
#' \itemize{
#' Opens a shiny interface that allows the user to vary model parameters and visualize a simulated colonization history.  Parameters include:
#' \item{\code{Scale of short-dist}} {The scale parameter for the short-distance component of the dispersal kernel, expressed as a percentage of cell width.}
#' \item{\code{Shape of short-dist}} {The shape parameter for the short-distance component of the dispersal kernel.}
#' \item{\code{Scale of long-dist}} {The mean distance for long-distance dispersal events, expressed as a fraction of cell width.}
#' \item{\code{proportion of long-dist}} {The mixture parameter, the proportion of dispersal events that involve long-distance movements.}
#' \item{\code{lambda}} {Growth rate of populations within cells.  Growth is exponential until carrying capacity of the cell is reached.}
#' \item{\code{carry capacity}} {The carrying capacity of grid cells on the landscape.  Assumes all cells are equally suitable (no ENM or pollen surface incorporated).}
#' \item{\code{Refuge?}} {One of several options for the location of glacial refugia.  Populations in 5 cells in each specified area are present at the start of demographic simulations.  Options include: Pennsylvania (PA), Georgia (GA), Texas (TX), or all of the above (ALL).}
#' }
#'
#' @export

forwardGUI <- function()
{
    shiny::runApp(paste0(system.file("extdata","shiny",package="holoSimCell")))
}

