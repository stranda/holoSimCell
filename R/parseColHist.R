#' parse the colhist object that comes out of getpophist.cells into coal history
#' 
#' @param ch colhist object
#'
#' @export

parseColhist <- function(ch)
{
    if (length(ch) == sum(sapply(ch,is.null))) #there are no colonizations
    {
        NULL
    } else {
    
    coaldemo <- data.frame(do.call(rbind,lapply(length(ch):1,function(x)
    {
        tch <- ch[[x]]
        if (!is.null(tch))
            if(dim(tch)[1]>0)
            {
                cbind(time=rep(x,dim(tch)[1]),
                      src=tch$pop,
                      snk=tch$source)
            } else NULL
    })))

    rl <- rle(coaldemo[with(coaldemo,order(src)),"src"])
    repeaters <- rl$values[rl$lengths>1]
    cd.norep <- coaldemo[!coaldemo$src %in% repeaters,]
    cd.rep <- coaldemo[coaldemo$src %in% repeaters,]
    cd.rep <- cd.rep[with(cd.rep,order(src,-time,snk)),]
    needed.rep <- data.frame(do.call(rbind,lapply(repeaters,function(r)
    {
        tdf <- cd.rep[cd.rep$src==r,]
        tdf$keep <- T
        if (dim(tdf)[1]>0)
        {
            tndf <- coaldemo[coaldemo$snk==r,]
            if (dim(tndf)[1]>0)
            {
                tndf$used <- F
                for (i in 1:dim(tdf)[1])
                {
                    before <- tndf[(tndf$time>=tdf$time[i])&(!tndf$used),]
                    if (dim(before)[1]>0)  ##there are some descendant pops from this one
                    {
                        tdf$keep[i] <- T
                        tndf$used[tndf$time%in%before$time] <- T
                    } else {
                        tdf$keep[i] <- F
                    }                    
                }
                ret <- tdf[tdf$keep,]
            } else {  #not ever an ancestor population
                ret <- tdf[tdf$time==min(tdf$time),]
            }
        }
        ret
    })))
    cd <- rbind(cd.norep,needed.rep[,-which(names(needed.rep)=="keep")])
        cd [with(cd,order(-time,src,snk)),]
        }
}


