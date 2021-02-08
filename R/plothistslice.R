###
### a plothist like function but one where you can specify a time point for the landscape history
### and plot all of the colonizations up to that point (as well as the habitat suitability at that point)
###
###

##' plot os slide of the landscape history
##' @param timeslice this is the layer of the landscape that corresponds to a time-click
##' @param ph this is the population history object
##' @param landscape this is the landscape object that drives the simulation
##' @export
plotHistSlice <- function(timeslice,ph,landscape,window=c(0,timeslice))
{
    if ((nrow(landscape$hab_suit)<timeslice) |( timeslice<1))
        stop ("need to specify a timeslice within the range of time clicks for this landscape")
    r <- landscape$sumrast
    values(r) <- landscape$hab_suit[timeslice,]
    pophist <- ph$pophist
    nv <- ph$Nvec
    xy <- data.frame(xyFromCell(r,unique(sort(pophist$pop))))
    names(xy) <- c("plotx","ploty")
    xy$pop <- unique(sort(pophist$pop))
    xy$ploty <- c(t(t(matrix(xy$ploty,nrow=dim(r)[2],ncol=dim(r)[1]))[dim(r)[1]:1,]))
    pophist <- merge(pophist,xy,att.x=T)
    pophist <- pophist[((!is.na(pophist$arrive)) & (pophist$arrive<=timeslice)),]
    pophist <- pophist[order(pophist$arrive),]

    opar=par()
    par(mar=c(2,2,2,0.5)+0.1)
    plot(flip(r,2), legend=F)

    szfunc <- function(x) {log(x+1)/6}
    
    for (i in 1:nrow(pophist)) #plot population locations
    {
        
        points(x=pophist$plotx[i],y=pophist$ploty[i],pch=16,col=1,cex=szfunc(nv[pophist$pop[i],timeslice]))
    }

    for (i in 1:nrow(pophist)) #plot arrows
        if ((pophist$arrive[i]>=min(window)) & (pophist$arrive[i]<=max(window)))
        {
        x1=pophist$plotx[i]
        y1=pophist$ploty[i]

        x0=pophist$plotx[pophist$pop==pophist$source[i]]
        y0=pophist$ploty[pophist$pop==pophist$source[i]]

        arrows(x0,y0,x1,y1, length=0.07, angle=15, code=2)
        }
    suppressWarnings(par(opar))
}


