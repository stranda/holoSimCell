##' hack to show the location of cells on landscape
##' @param cells a vector of cell ids to plot
##' @param landscape object containing habitat suitabiliities
##' @export
plotCellsOnSuit <- function(cells,landscape,timeslice=1)
{
    hs=landscape$hab_suit[timeslice,]
    hs=as.numeric(hs>0)
    hs[cells] <- 2
    m <- matrix(hs,nrow=landscape$details$x.dim,ncol=landscape$details$y.dim,byrow=F)
    image(x=c(1:landscape$detail$x.dim),y=c(1:landscape$details$y.dim),z=m,xlab="Columns",ylab="Rows")
}
