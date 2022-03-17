#' Choose refuge cells from a plot of the landscape used in the simulations.
#'
#' Interactive function allowing a user to identify specific cells on the landscape, which are subsequently used as starting locations for a forward-time demographic simulation.
#'
#' @param r this is the surface (a single raster layer)
#'
#' @details Interactive function using locator() to identify cell IDs from user input. Assumes refugia include a central cell and adjacent cells in each cardinal direction (5 total).
#'
#' @return A vector of five cell IDs associated with the user-specified refuge location
#'
#' @examples
#' library(holoSimCell)
#' rs <- raster::brick(paste0(system.file("extdata","rasters",package="holoSimCell"),"/","study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif"))
#' chooseRefugeCells(rs[[1]])
#'
#' @seealso \code{\link[graphics]{locator}}, \code{\link{getpophist2.cells}}
#'
#' @export

chooseRefugeCells <- function(r)
{
    if ("RasterLayer" %in% class(r))
    {
        pr = r
        raster::plot(pr,main="Press mouse #1 for center of refuge point")
        p <- locator(type="p",n=1)
        rw = seq(nrow(pr),1)[rowFromY(pr, p$y)]
        cl = colFromX(pr, p$x)

        rcells = c(
            cellFromRowCol(pr, rw,cl),
            cellFromRowCol(pr, rw,cl+1),
            cellFromRowCol(pr, rw,cl-1),
            cellFromRowCol(pr, rw+1,cl),
            cellFromRowCol(pr, rw-1,cl)
        )
        pr[c((seq(nrow(pr),1)[rw]-1):(seq(nrow(pr),1)[rw]+1)),cl] <- 1
        pr[seq(nrow(pr),1)[rw],c((cl-1):(cl+1))] <- 1
        raster::plot(pr,main="Refugium cells")
        return(rcells)
    } else {
        stop("r must be a raster layer")
    }

}


#