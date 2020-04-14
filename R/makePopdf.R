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
