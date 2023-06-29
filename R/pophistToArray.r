#' Convert "pophist" object to an array
#'
#' This function converts a "pophist" object to an array. This function is helpful for converting a "pophist" object to input for the \code{\link[enmSdm]{bioticVelocity}} function.
#' @param x A pophist object (a list).
#' @param longitude Either \code{NULL} (default) or a vector of longitudes, one per column in the output object. Note that if \code{longitude} is \code{NULL} and the "pophist" object in \code{x} contains a field named \code{x$pophist$longitude}, then this will be used (and will override the values in \code{longitude} if it is specified).
#' @param latitude Either \code{NULL} (default) or a vector of latitudes, one per row in the output object. Note that if \code{latitude} is \code{NULL} and the "pophist" object in \code{x} contains a field named \code{x$pophist$latitude}, then this will be used (and will override the values in \code{latitude} if it is specified).
#' @param times Either \code{NULL} (default) or a numeric vector. This specifies the time period represented by each column in \code{x}. Times \emph{must} appear in sequential order. For example, if time periods are 24 kybp, 23 kybp, 22 kybp, use \code{c(-24, -23, -22)}, \emph{not} \code{c(24, 23, 22)}. Note that if \code{times} is \code{NULL} and the "pophist" object in \code{x} contains a field named \code{x$pophist$times}, then this will be used (and will override the values in \code{times} if it is specified).
#' @param warn Logical, if \code{TRUE} then display function-specific warnings.
#'
#' @details Utility function for pre-processing an Nvecs matrix from a pophist object for biotic velocity calculations.
#'
#' @return A list object with:
#' \itemize{
#' 	\item \code{pophistAsArray}: Array representing values in \code{x$pophist} with one "layer" (3rd dimension) per time period.
#' 	\item \code{longitude}: Array representing longitudes.
#' 	\item \code{latitude}: Array representing latitudes.
#' 	\item \code{times}: Vector representing times.
#' }
#'
#' @examples
#' library(holoSimCell)
#' parms <- drawParms(control = system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell"))
#' load(file=paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",pollenPulls[[1]]$file))
#' refpops <- pollenPulls[[1]]$refs
#' avgCellsz <- mean(c(res(landscape$sumrast)))
#'
#' ph = getpophist2.cells(h = landscape$details$ncells, xdim = landscape$details$x.dim, ydim = landscape$details$y.dim,
#'                        landscape=landscape,
#'                        refs=refpops,   
#'                        refsz=parms$ref_Ne,
#'                        lambda=parms$lambda,
#'                        mix=parms$mix,  
#'                        shortscale=parms$shortscale*avgCellsz,  
#'                        shortshape=parms$shortshape, 
#'                        longmean=parms$longmean*avgCellsz,  
#'                        ysz=res(landscape$sumrast)[2], 
#'                        xsz=res(landscape$sumrast)[1], 
#'                        K = parms$Ne) 
#'
#' times_1k <- seq(-21000,0,by=990)
#'    
#' pharray <- pophistToArray(NvecNAs(ph, landscape), times = times_1G)
#'    
#' metrics <- c('centroid', 'nsQuants', 'summary')
#'
#' BV_permill_all <- bioticVelocity(
#'   x=pharray$pophistAsArray,
#'   times = times_1G,
#'   atTimes = times_1k,
#'   longitude=pharray$longitude,
#'   latitude=pharray$latitude,
#'   metrics=metrics,
#'   quants=c(0.05, 0.1, 0.9, 0.95),
#'   onlyInSharedCells = FALSE)
#'
#' @seealso \code{\link{getpophist2.cells}}, \code{\link{NvecNAs}}, \code{\link[enmSdm]{bioticVelocity}}
#'
#' @export
pophistToArray <- function(
	x,
	longitude = NULL,
	latitude = NULL,
	times = NULL,
	warn = TRUE
) {

	# time periods
	totalTimes <- ncol(x$Nvecs)

	if (!is.null(x$pophist$times)) {
		times <- x$pophist$times
	} else if (is.null(times)) {
		if (warn) warning('No values for time were supplied and the population history object\ndoes not contain a field named ".$pophist$times". Arbitrary values of\ntime will be assigned.')
		times <- 1:totalTimes
	}

	# populate array
	nRows <- max(x$pophist$row)
	nCols <- max(x$pophist$col)
	xNew <- array(NA, dim=c(nRows, nCols, totalTimes))
	for (i in 1L:totalTimes) xNew[, , i] <- matrix(x[['Nvecs']][ , i], nrow=nRows, ncol=nCols, byrow=TRUE)
	
	# longitude
	if (is.null(longitude[1L]) & is.null(x$pophist$longitude)) {

		longs <- x$pophist$col[!duplicated(x$pophist$pop)]
		longitude <- matrix(longs, nrow=nRows, ncol=nCols, byrow=TRUE)
		if (warn) warning('No values for latitude were supplied and the population history object\ndoes not contain a field named ".$pophist$longitude". Arbitrary values of\nlongitude will be assigned.')
	
	} else if (is.null(longitude[1L]) & !is.null(x$pophist$longitude)) {
		longs <- x$pophist$longitude[!duplicated(x$pophist$pop)]
		longitude <- matrix(longs, nrow=nRows, ncol=nCols, byrow=TRUE)
	} else if (!is.null(longitude[1L]) & is.null(x$pophist$longitude)) {
		longitude <- matrix(longitude, nrow=nRows, ncol=nCols, byrow=TRUE)
	}
	
	# latitude
	if (is.null(latitude[1L]) & is.null(x$pophist$latitude)) {

		lats <- x$pophist$row[!duplicated(x$pophist$pop)]
		latitude <- matrix(lats, nrow=nRows, ncol=nCols, byrow=TRUE)
		if (warn) warning('No values for latitude were supplied and the population history object\ndoes not contain a field named ".$pophist$latitude". Arbitrary values of\nlatitude will be assigned.')
	
	} else if (is.null(latitude[1L]) & !is.null(x$pophist$latitude)) {
		lats <- x$pophist$latitude[!duplicated(x$pophist$pop)]
		latitude <- matrix(lats, nrow=nRows, ncol=nCols, byrow=TRUE)
	} else if (!is.null(latitude[1L]) & is.null(x$pophist$latitude)) {
		latitude <- matrix(latitude, nrow=nRows, ncol=nCols, byrow=TRUE)
	}
	
	list(pophistAsArray=xNew[c(nrow(xNew[,,1]):1),,], longitude=longitude, latitude=latitude, times=times)
	
}

