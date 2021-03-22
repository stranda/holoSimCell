##
## hack to show the location of cells on landscape
##

plotCellsOnSuit <- function(cells,landscape,timeslice=1)
{
    hs=landscape$hab_suit[timeslice,]
    hs=as.numeric(hs>0)
    hs[cells] <- 2
    m <- matrix(hs,nrow=landscape$details$x.dim,ncol=landscape$details$y.dim,byrow=F)
    image(m)
}
