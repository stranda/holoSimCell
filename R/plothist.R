#'
#' plot the population history
#' @param ph population history object (direct output from gepophist.cells)
#' @export
#' 
plothist <- function(ph, maxtime=NULL)
{
    ch=ph$coalhist
    pops=ph$pophist
    pops=data.frame(do.call(rbind,lapply(unique(pops$pop),function(p)  #only keep the most recent col event
    {
        if (length(which(pops==p))>1)
        {
            tmp <- pops[pops$pop==p,]
            tmp <- tmp[order(-tmp$arrive),][1,]
        } else pops[pops$pop==p,]
    })))
    
    ch=merge(ch,pops,by.x="src",by.y="pop")[,c(-6,-7)]
    ch=ch[with(ch,order(-time,src,snk)),]
    ch=rbind(data.frame(src=pops$pop[(!is.na(pops$arrive))&(pops$arrive==0)],
                        time=pops$arrive[(!is.na(pops$arrive))&(pops$arrive==0)],
                        snk=pops$source[(!is.na(pops$arrive))&(pops$arrive==0)],
                        col=pops$col[(!is.na(pops$arrive))&(pops$arrive==0)],
                        row=pops$row[(!is.na(pops$arrive))&(pops$arrive==0)]),
             ch)
                        
   
    layout(matrix(c(2,1,1,1,
                    2,1,1,1,
                    2,1,1,1,
                    0,3,4,5),nrow=4,ncol=4,byrow=T),
           widths=c(0.22,0.26,0.26,0.26)) #plot 3 as null

    
    rows <- min(pops$row):max(pops$row)
    rw <- max(rows)
    
    cols <- min(pops$col):max(pops$col)
    cl <- max(cols)
    
    ch$coldist=NA

    
    mai = par("mai")
    par(mai=c(0.1,0.1,0.1,0.1))
    
    plot(1,type="n",ylab="",xlab="",ylim=range(rows),xlim=range(cols),axes=F,asp=1)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray") #gray background
###    points(row~col,pops,pch=16,cex=log(dim(pops)[1])/6)

    ((pops$pop-1) %/% cl) +1->row
    ((pops$pop-1) %% cl) +1 ->col
    
    if ((!is.null(ph$Nvecs))&(min(dim(ph$Nvecs))>0))
    {
        harmmean <- function(x){1/mean(1/x)}
        sz <- apply(ph$Nvecs+1,1,harmmean)-1
        pops$sz <- sz
    } else {sz <- log(dim(pops)[1])/6}

    points(row~col,pops,pch=16,cex=log(sz+1))

    for (i in 1:dim(ch)[1])
    {
        if (!is.na(ch$snk[i]))
            {
                x1=ch[i,"col"]
                y1=ch[i,"row"]
                x0= ch[ch$src==ch[i,"snk"],"col"][1]
                y0= ch[ch$src==ch[i,"snk"],"row"][1]
                                        #            print(c(x0,y0,x1,y1))
                if (is.null(maxtime))
                    clrs = heat.colors(max(ch$time,na.rm=T))[ch$time[i]]  else
                    clrs = heat.colors(maxtime)[ch$time[i]]
                
                arrows(x0,y0,x1,y1,
                       col=clrs,
                       lwd=2,
                       length=0.1
                       )
                ch$coldist[i] <- sqrt((y0-y1)^2 + (x0-x1)^2)
                                        #           print(ch$coldist[i])
            }
    }
    mxt <- ifelse(is.null(maxtime),max(ch$time,na.rm=T),maxtime)
steps=5
    lgd <- round(seq(1,mxt,length.out=steps))

    par(mai=c(0.5,0.6,0.5,0.6))
    
    barplot(matrix(rep(diff(lgd)[1],steps),ncol=1),col=heat.colors(mxt)[lgd],axes=F,
            main="Event Age \n(from start)")
    for (lvl in 1:steps)
    {
        text(x=0.5,y=lgd[lvl]+(diff(lgd)[1]/2),labels=lgd[lvl],cex=2,adj=0.5)
        }
#    plot(1~1,type="p",axes=F,xlab="",ylab="")
#    legend(x=-2, y=max(ch[,"row"])/2, legend=lgd, col=heat.colors(mxt)[lgd],lwd=2)
par(mai=mai)

    hist(ch$coldist,xlab="Colonization distance",main="Realized colonization distances")
    hist(ch$time,xlab="Colonization time",main="Realized colonization times")


    cent <- bioticVelocity(ph,metrics="centroid")$centroidVelocity
    ncent <- bioticVelocity(ph,metrics="nCentroid")$nCentroidVelocity
    both <- c(cent,ncent)
    both <- both[!is.na(both)|is.finite(both)]

    plot(cent,type="l",
         xlab="Time (bigger more recent)",ylab="Centroid velocity",
         main="Range movement",ylim=range(c(both)))
    points(ncent,type="l",col="red")
    legend(x=0.4*length(cent),y=0.6*max(both),legend=c("Center","North margin"),col=c(1,2),lty=c(1,1))

    layout(matrix(1))
    
}


#'
#' plot the population history
#' @param pops population history
#' @export
#' 
plothist.old <- function(ph)
{
    pops=ph$pophist
    
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
