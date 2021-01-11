#' Rescale carrying capacity rasters and assign refugia
#'
#' This function takes an "abundance" raster (i.e., from an ENM or from a pollen surface) and identifies refugia and starting (relative) abundances for each refugium. It rescales this to the extent and resolution of a "simulation" raster which typically has coarser spatial resolution than the abundance raster.
#'
#' @param abund Abundance raster.
#' @param sim Simulation raster.
#' @param quant Numeric or \code{NULL}. Quantile at which to threshold the abundance raster. Values that fall above this threshold will be assumed to represent a refuge.
#' @param threshold Numeric or \code{NULL}. Value at which to threshold the abundance raster. Values that fall above this threshold will be assumed to represent a refuge. Note that you can either use \code{quant} or \code{threshold}, but if you specify both then \code{threshold} will be used instead (with a warning).
#' @param intermediate Logical. If \code{FALSE} (default), return a raster stack with refuge ID and abundance with the same extent and resolution as the simulation raster. If \code{TRUE}, return a list. One item in the list is the raster stack previously described, and the other is a stack with refuge ID and abundance at the same extent and resolution as the abundance raster.
#'
#' @details
#' This function attempts to rescale the identify of refugia and abundances obtained from the abundance raster to a new spatial resolution and extent. Since the simulation raster often has cells that are larger than the cells in the abundance raster, in some cases it cannot faithfully retain abundances or even all unique refugia identified in the abundance raster. The procedure first thresholds the abundance raster using either a quantile or a threshold.  All cells that are equal to or above this threshold are assumed to be a potential refuge.  It then creates a unique "abundId" raster by assigning a unique integer number for each block of contiguous cells (using Moore neighborhood adjacency). \cr
#' Then, a "simAbund" raster is created with the same extent and resolution as the simulation raster. For each each cell in this raster, the function determines if it contains at least one cell in the "abundId" raster that is assigned to a refuge. The challenge here is that a single "simAbund" cell can contain cells that are assigned to multiple refugia in the "abundId" raster, and that "simAbund" cell can also include cells that have abundances that are assigned to no refugia in the abundance raster. Thus, if we simply assigned abundances to the "simAbund" cell by resampling the abundance raster, we would in some cases be too generous because a single "simAbund" cell can include cells that do not belong to this refuge. \cr
#' The procedure assigns abundances by first calculating a proportionality scalar where the numerator is the sum of abundances of abundance raster cells in this refugium and in the "simAbund" cell, and the denominator the sum of all abundances of all cells in this "simAbund" cell. The abundance assigned to this "simAbund" cell for this particular refugium is this scalar times the abundance from the resampling of the abundance raster to the extent/resolution of the simulation raster. Thus, abundances assigned to any particular cell in a refuge will be equal to or less than the abundance of the resampled values. \cr
#' The procedure then assigns each cell an integer number identifying which refugium to belongs to and an abundance corresponding to the given refuge. When cells contain more than one "abundId" refuge cell, the refuge with the greater abundance is assigned to the cell.  As a consequence, a refuge that appears in the "abundId" raster could be trimmed in extent or even eliminated if it is only represented by a few cells that have small abundances relative to a more "massive" refugium in the same cell. Also, as a result, it is possible to have distinct refugia in cells that are adjacent to one another when rescaled to the extent/resolution of the simulation raster but are spatially distinct at the scale of the abundance raster.
#' @return Raster stack or a list of raster stacks.
#' @examples
#' 
#' abund <- stack('C:/Ecology/Drive/Research/ABC vs #' Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/predictions/ccsm_320kmExtent_brt.tif')
#' abund <- abund[[1]]
#' 
#' sim <- raster('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_resampled_to_genetic_demographic_simulation_resolution.tif')
#' 
#' rasts <- assign_refugia_from_abundance_raster(abund, sim, 0.95, intermediate=TRUE)
#' 
#' cols <- c('red', 'orange', 'yellow', 'green', 'blue', 'purple', 'gray',
#' 	'chartreuse', 'darkgreen', 'cornflowerblue', 'goldenrod3', 'black',
#' 	'steelblue3', 'forestgreen', 'pink', 'cyan', 'darkred')
#' 
#' par(mfrow=c(2, 2))
#' 
#' # plot(abund, main='abundance at native resolution')
#' col <- cols[1:cellStats(rasts$abund[['id']], 'max')]
#' plot(rasts$abund[['id']], col=col, main='refugia ID at native resolution')
#' plot(rasts$abund[['abundRefuge']], main='relative abundance in refugia at native resolution')
#' col <- cols[1:cellStats(rasts$sim[['id']], 'max')]
#' plot(rasts$sim[['id']], col=col, main='refuge ID at simulation resolution')
#' plot(rasts$sim[['abundance']], main='relative abundance at simulation resolution')
#' 
#' @export
assign_refugia_from_abundance_raster <- function(
	abund,
	sim,
	quant = NULL,
	threshold = NULL,
	intermediate = FALSE
) {

	if (!is.null(quant) & !is.null(threshold)) warning('Both "quant" and "threshold" are not NULL.\nThresholding abundance raster using "threshold".')

	### abundance raster
	######################

	# create data frame for abundance raster with:
	# cell number in abundance frame
	# long, lat
	# a mask of refugia
	# ID number of each refuge
	# abundance (abundance)
	# cell number in simulation raster
	
	# "abundance"
	names(abund) <- 'abundAbund'
	
	# mask refugia
	if (is.null(threshold)) threshold <- quantile(abund, quant)
	maskAbund <- abund >= threshold
	names(maskAbund) <- 'mask'
	
	# identify refugia
	idsAbund <- raster::clump(maskAbund, directions=8, gaps=FALSE)
	names(idsAbund) <- 'id'
	numRefugia <- raster::cellStats(idsAbund, 'max')
	
	# abundance in refugial cells
	abundRefugeAbund <- abund * maskAbund
	names(abundRefugeAbund) <- 'abundRefuge'
	
	# long/lat
	ll <- enmSdm::longLatRasters(abund)

	# cell numbers
	cellsAbund <- raster::setValues(abund, 1:raster::ncell(abund))
	names(cellsAbund) <- 'cellNumAbund'
	
	suitStack <- raster::stack(cellsAbund, ll, maskAbund, idsAbund, abund, abundRefugeAbund)
	suitFrame <- as.data.frame(suitStack)

	suitFrame$cellNumSim <- raster::cellFromXY(sim, suitFrame[ , c('longitude', 'latitude')])

	### simulation raster
	#####################
	
	# create data frame with:
	# simulation raster cell number
	# abundance resampled from abundance raster
	# sum of abundances of cells from abundance raster, by refuge ID
	# refuge ID
	# rescaled abundances (to match resampled abundances) for the refuge to which this cell is assigned accounting for:
	#	* empty cells (outside a refuge in the abundance raster but with non-NA values)
	#   * abundance cells that fall into other refugia (should not be counted toward this cell's abundance)

	# "abundance"
	abundSim <- raster::resample(abund, sim)
	names(abundSim) <- 'abundSim'
	
	# cell numbers
	cellsSim <- raster::setValues(sim, 1:raster::ncell(sim))
	names(cellsSim) <- 'cellNumSim'
	
	simStack <- raster::stack(cellsSim, abundSim)
	
	simFrame <- as.data.frame(simStack)
	
	### calculate abundances in simulation raster for each refuge
	for (countRefuge in 1:numRefugia) {
	
		simFrame$DUMMY <- NA
		names(simFrame)[ncol(simFrame)] <- paste0('refuge', countRefuge)
		
		suitFrameOutsideRefuge <- suitFrame[!is.na(suitFrame$abundAbund) & (is.na(suitFrame$id) | omnibus::naCompare('!=', suitFrame$id, countRefuge)), ]
		suitFrameInsideRefuge <- suitFrame[omnibus::naCompare('==', suitFrame$id, countRefuge), ]
		
		# assign scaled abundances to each simulation raster cell in this refuge
		# For each simulation cell that overlaps with at least one cell in the abundance raster assigned to a particular refuge, find:
		#		* a proportionality factor as a proportion of the sum of suitabilities in the abundance raster cells in this refuge divided by the sum of all suitabilities in all cells that fall into this simulation raster cell (regardless of whether they're in a refuge or not)
		#		* the resampled abundance (resampling the abundance raster to the simulation raster)
		# Final abundance in this cell for a particular refuge refuge is the product of these two values.
		simCellsInRefuge <- sort(unique(suitFrameInsideRefuge$cellNumSim))
		for (simCell in simCellsInRefuge) {
		
			simCellInterpAbund <- simFrame$abundSim[simFrame$cellNumSim == simCell]
			abundInSimCellOutsideRefuge <- suitFrameOutsideRefuge$abundAbund[suitFrameOutsideRefuge$cellNumSim == simCell]
			abundInSimCellInsideRefuge <- suitFrameInsideRefuge$abundAbund[suitFrameInsideRefuge$cellNumSim == simCell]
			
			if (length(abundInSimCellOutsideRefuge) == 0) {
				simAbund <- simCellInterpAbund
			} else {
				insideAbund <- sum(abundInSimCellInsideRefuge)
				outsideAbund <- sum(abundInSimCellOutsideRefuge)
				totalAbund <- insideAbund + outsideAbund
				simAbund <- simCellInterpAbund * (insideAbund / totalAbund)
			}
			
			simFrame[simFrame$cellNumSim == simCell, paste0('refuge', countRefuge)] <- simAbund
			
		
		} # next simulation raster cell
	
	} # next refuge
	
	### assign refuge ID and abundance
	simIdsList <- apply(simFrame[ , paste0('refuge', 1:numRefugia), drop=FALSE], 1, which.max)
	
	simFrame$id <- NA
	for (i in seq_along(simIdsList)) simFrame$id[i] <- ifelse(length(simIdsList[[i]]) == 0, NA, simIdsList[[i]])
	
	simFrame$refugeAbund <- NA
	for (countRefuge in 1:numRefugia) {
	
		refugeIndex <- which(omnibus::naCompare('==', simFrame$id, countRefuge))
		simFrame$refugeAbund[refugeIndex] <- simFrame[refugeIndex, paste0('refuge', countRefuge)]
	
	}

	idsSim <- abundSim <- NA * sim
	idsSim <- raster::setValues(idsSim, simFrame$id)
	abundSim <- raster::setValues(abundSim, simFrame$refugeAbund)
	names(idsSim) <- 'id'
	names(abundSim) <- 'abundance'
	
	if (intermediate) {
		list(sim=raster::stack(idsSim, abundSim), abund=suitStack[[c('id', 'abundRefuge')]])
	} else {
		raster::stack(idsSim, abundSim)
	}
	
}
