#' converts cells in Nvec form to NA if their underlying hab suit is NA for that time point
#'
#' @param ph pophist object
#' @param l landscape containing the rasters
#'
#' @export
NvecNAs <- function(ph,l)
{
    if (is.null(ph)) {stop("must supply a population history")}
    if (is.null(l)) {stop("must supply a landscape")}
    if (dim(ph$Nvecs)[1] != (prod(dim(l$NAmask)[1:2]))) {stop("the number of cells in ph and l are different")}
    if (dim(ph$Nvecs)[2] != (dim(l$NAmask)[3])) {stop("the number of time clicks in ph and l are different")}

    m=raster::as.array(l$NAmask)
    m2=base::array(dim=c(dim(m)[c(2,1,3)]))
    for (tm in 1:dim(m)[3])
    {
        m2[,,tm] <- t(m[,,tm])
        m2[,,tm] <- m2[,dim(m2)[2]:1,tm]
    }
    hs = t(l$hab_suit)
    hs = hs[,c(1,1:ncol(hs))]
    nvo[is.na(hs)] <- NA
    nvo
}
