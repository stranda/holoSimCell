#'
#' plot the population history
#' @param pops population history
#' @export
#' 
plothist <- function(pops)
{
    layout(matrix(c(1,1,1,1,2,3),nrow=3,ncol=2,byrow=T))
    rows <- min(pops$row):max(pops$row)
    cols <- min(pops$col):max(pops$col)
    pops$coldist=NA
    plot(1,type="n",ylab="",xlab="",ylim=range(rows),xlim=range(cols),axes=F,asp=1)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray") #gray background
    points(row~col,pops,pch=16,cex=log(dim(pops)[1])/6)
    for (i in 1:dim(pops)[1])
    {
        if ((!is.na(pops$arrive[i]))&(pops$arrive[i]!=0))
        {
            
            x1=pops[i,"col"]
            y1=pops[i,"row"]
            x0= pops[pops$pop==pops[i,"source"],"col"]
            y0= pops[pops$pop==pops[i,"source"],"row"]
#            print(c(x0,y0,x1,y1))
            arrows(x0,y0,x1,y1,
                   col=heat.colors(max(pops$arrive,na.rm=T))[pops$arrive[i]],
                   lwd=2,
                   length=0.1
               )
            pops$coldist[i] <- sqrt((y0-y1)^2 + (x0-x1)^2)
        }
    }

    hist(pops$coldist,xlab="Colonization distance",main="")
    hist(pops$arrive,xlab="Colonization time",main="")
    layout(matrix(1))
}
