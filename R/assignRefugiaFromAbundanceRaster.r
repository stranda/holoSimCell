#' Rescale carrying capacity rasters and assign refugia
#'
#' This function takes an "abundance" raster (i.e., from an ENM or from a pollen surface) and identifies refugia and starting (relative) abundances for each refugium. It rescales this to the extent and resolution of a "simulation" raster which typically has coarser spatial resolution than the abundance raster. It then generates a vector of cell IDs that correspond to refugia cells and calculates the average relative abundance across refugial cells.
#'
#' @param abund Abundance raster.
#' @param sim Simulation raster.
#' @param threshold Numeric. Value at which to threshold the abundance raster. Values that fall above this threshold will be assumed to represent a refuge and values below will be assumed to be outside (i.e., zero abundance).
#'
#' @details
#' OLD, no longer valid: This function rescales abundances obtained from the abundance raster to a new spatial resolution and extent used for demographic simulations. Since the demographic simulation raster often has cells that are larger than the cells in the abundance raster, it cannot faithfully retain abundances or even all unique refugia identified in the abundance raster. The procedure first thresholds the abundance raster using a user-defined value above which it is assumed a cell was inside a refuge and below which it was assumed to be outside any refuge (zero abundance).  A unique "abundId" raster is then created by assigning a unique integer number for each block of contiguous cells (using Moore neighborhood adjacency). This raster is then resampled to the resolution used by the simulation raster. The renumbering is redone so that cells that refugial were not adjacent at the original resolution but are in the new resolution are assigned to the same refuge. \cr
#' OLD, no longer valid: Then, a "simAbund" raster is created with the same extent and resolution as the simulation raster. For each each cell in this raster, the function determines if it contains at least one cell in the "abundId" raster that is assigned to a refuge. The challenge here is that a single "simAbund" cell can contain cells that are assigned to multiple refugia in the "abundId" raster, and that "simAbund" cell can also include cells that have abundances that are assigned to no refugia in the abundance raster. Thus, if we simply assigned abundances to the "simAbund" cell by resampling the abundance raster, we would in some cases be too generous because a single "simAbund" cell can include cells that do not belong to this refuge. \cr
#' OLD, no longer valid: The procedure assigns abundances by first calculating a proportionality scalar where the numerator is the sum of abundances of abundance raster cells in this refugium and in the "simAbund" cell, and the denominator the sum of all abundances of all cells in this "simAbund" cell. The abundance assigned to this "simAbund" cell for this particular refugium is this scalar times the abundance from the resampling of the abundance raster to the extent/resolution of the simulation raster. Thus, abundances assigned to any particular cell in a refuge will be equal to or less than the abundance of the resampled values. \cr
#' OLD, no longer valid: The procedure then assigns each cell an integer number identifying which refugium to belongs to and an abundance corresponding to the given refuge. When cells contain more than one "abundId" refuge cell, the refuge with the greater abundance is assigned to the cell.  As a consequence, a refuge that appears in the "abundId" raster could be trimmed in extent or even eliminated if it is only represented by a few cells that have small abundances relative to a more "massive" refugium in the same cell. Also, as a result, it is possible to have distinct refugia in cells that are adjacent to one another when rescaled to the extent/resolution of the simulation raster but are spatially distinct at the scale of the abundance raster.
#' @return A list with: 1) A raster stack representing refuge ID numbers and abundances at the \emph{simulation} resolution and extent; 2) a vector of cell numbers for refugial cells at the \emph{simulation} scale; and 3) a single numeric value representing mean refugial abundance across cells at the \emph{simulation} extent.
#'
#' @examples
#' 
#' library(raster)
#' abund <- stack('E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/predictions/ccsm_320kmExtent_brt.tif')
#' abund <- abund[[1]]
#' 
#' sim <- raster('E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_resampled_to_genetic_demographic_simulation_resolution.tif')
#' 
#' threshold <- 0.6
#' refs <- assignRefugiaFromAbundanceRaster(abund, sim, threshold)
#' 
#' cols <- c('red', 'orange', 'yellow', 'green', 'blue', 'purple', 'gray',
#' 	'chartreuse', 'darkgreen', 'cornflowerblue', 'goldenrod3', 'black',
#' 	'steelblue3', 'forestgreen', 'pink', 'cyan', 'darkred')
#' 
#' par(mfrow=c(1, 2))
#' 
#' col <- cols[1:maxValue(refs$sim[['refugiaId']])]
#' plot(refs$simulationScale[['refugiaId']], col=col, main='refuge ID')
#' plot(refs$simulationScale[['refugiaAbund']], main='refuge abundance')
#' 
#' @export
assignRefugiaFromAbundanceRaster <- function(
	abund,
	sim,
	threshold
) {

	# ### abundance raster
	# ######################

	# # create data frame for abundance raster with:
	# # cell number in abundance frame
	# # long, lat
	# # a mask of refugia
	# # ID number of each refuge
	# # abundance (abundance)
	# # cell number in simulation raster

	# # "abundance"
	# names(abund) <- 'origAbund'
	
	# # mask refugia
	# origMask <- abund >= threshold
	# names(origMask) <- 'origMask'
	
	# # identify refugia
	# idsOrig <- raster::clump(origMask, directions=8, gaps=FALSE)
	# names(idsOrig) <- 'idOrig'

	# # abundance in refugial cells
	# abundRefugeAbund <- abund * origMask
	# names(abundRefugeAbund) <- 'abundOrigRefuge'
	
	# # long/lat
	# ll <- enmSdm::longLatRasters(abund)

	# # re-assess ID number based on resolution of simulation raster
	# idsSim <- raster::resample(idsOrig, sim)
	# idsSim <- raster::clump(idsSim, directions=8, gaps=FALSE)
	# idsSim <- raster::extract(idsSim, raster::as.data.frame(ll))
	# numRefugia <- max(idsSim, na.rm=TRUE)

	# # cell numbers
	# cellsAbund <- raster::setValues(abund, 1:raster::ncell(abund))
	# names(cellsAbund) <- 'cellNumOrig'
	
	# suitStack <- raster::stack(cellsAbund, ll, origMask, idsOrig, abund, abundRefugeAbund)
	# suitFrame <- as.data.frame(suitStack)

	# suitFrame$idSim <- idsSim
	
	# # cell numbers of simulation raster in the layout used by the demo/genetic simulations:
	# # cell in lower left is 1, to its right is 2, etc, then wraps to next row so cell [nrow - 1, 1] is next in line
	# nrows <- nrow(sim)
	# ncols <- ncol(sim)
	# ncells <- raster::ncell(sim)

	# v <- rep(seq(nrows * (ncols - 1) - 1, 1, by=-ncols), each=ncols) + 0:(ncols - 1)
	# cellNumSim <- matrix(v, nrow=nrows, ncol=ncols, byrow=TRUE)
	# cellNumSim <- raster::raster(cellNumSim, template=sim)
	# suitFrame$cellNumSim <- raster::extract(cellNumSim, suitFrame[ , c('longitude', 'latitude')])
		
	# ### simulation raster
	# #####################
	
	# # create data frame with:
	# # simulation raster cell number
	# # abundance resampled from abundance raster
	# # sum of abundances of cells from abundance raster, by refuge ID
	# # refuge ID
	# # rescaled abundances (to match resampled abundances) for the refuge to which this cell is assigned accounting for:
	# #	* empty cells (outside a refuge in the abundance raster but with non-NA values)
	# #   * abundance cells that fall into other refugia (should not be counted toward this cell's abundance)

	# # "abundance"
	# abundSim <- raster::resample(abund, sim)
	# names(abundSim) <- 'abundSim'
	
	# # cell numbers
	# cellsSim <- raster::setValues(sim, 1:raster::ncell(sim))
	# names(cellsSim) <- 'cellNumSim'
	
	# simStack <- raster::stack(cellsSim, abundSim)
	# simFrame <- as.data.frame(simStack)
	
	# ### calculate abundances in simulation raster for each refuge
	# for (countRefuge in 1:numRefugia) {
	
		# simFrame$DUMMY <- NA
		# names(simFrame)[ncol(simFrame)] <- paste0('refuge', countRefuge)
		
		# suitFrameOutsideRefuge <- suitFrame[!is.na(suitFrame$origAbund) & (is.na(suitFrame$idSim) | omnibus::naCompare('!=', suitFrame$idSim, countRefuge)), ]
		# suitFrameInsideRefuge <- suitFrame[omnibus::naCompare('==', suitFrame$idSim, countRefuge), ]
		
		# # assign scaled abundances to each simulation raster cell in this refuge
		# # For each simulation cell that overlaps with at least one cell in the abundance raster assigned to a particular refuge, find:
		# #		* a proportionality factor as a proportion of the sum of suitabilities in the abundance raster cells in this refuge divided by the sum of all suitabilities in all cells that fall into this simulation raster cell (regardless of whether they're in a refuge or not)
		# #		* the resampled abundance (resampling the abundance raster to the simulation raster)
		# # Final abundance in this cell for a particular refuge refuge is the product of these two values.
		# simCellsInRefuge <- sort(unique(suitFrameInsideRefuge$cellNumSim))
		# for (simCell in simCellsInRefuge) {
		
			# simCellInterpAbund <- simFrame$abundSim[simFrame$cellNumSim == simCell]
			# abundInSimCellOutsideRefuge <- suitFrameOutsideRefuge$origAbund[suitFrameOutsideRefuge$cellNumSim == simCell]
			# abundInSimCellInsideRefuge <- suitFrameInsideRefuge$origAbund[suitFrameInsideRefuge$cellNumSim == simCell]
			
			# if (length(abundInSimCellOutsideRefuge) == 0) {
				# simAbund <- simCellInterpAbund
			# } else {
				# insideAbund <- sum(abundInSimCellInsideRefuge)
				# outsideAbund <- sum(abundInSimCellOutsideRefuge)
				# totalAbund <- insideAbund + outsideAbund
				# simAbund <- simCellInterpAbund * (insideAbund / totalAbund)
			# }
			
			# simFrame[simFrame$cellNumSim == simCell, paste0('refuge', countRefuge)] <- simAbund
			
		
		# } # next simulation raster cell
	
	# } # next refuge
	
	# ### assign refuge ID and abundance
	# simIdsList <- apply(simFrame[ , paste0('refuge', 1:numRefugia), drop=FALSE], 1, which.max)
	
	# simFrame$id <- NA
	# for (i in seq_along(simIdsList)) simFrame$id[i] <- ifelse(length(simIdsList[[i]]) == 0, NA, simIdsList[[i]])
	
	# simFrame$refugeAbund <- NA
	# for (countRefuge in 1:numRefugia) {
	
		# refugeIndex <- which(omnibus::naCompare('==', simFrame$id, countRefuge))
		# simFrame$refugeAbund[refugeIndex] <- simFrame[refugeIndex, paste0('refuge', countRefuge)]
	
	# }

	# # set values of ID and abundance rasters
	# idsSimRast <- abundSim <- NA * sim
	# idsSimRast <- raster::setValues(idsSimRast, simFrame$id)
	# abundSim <- raster::setValues(abundSim, simFrame$refugeAbund)
	
	# # flip up/down because we've renumbered simulation raster cell numbers as per demographic/genetic simulations
	# idsSimRast <- raster::as.matrix(idsSimRast)
	# abundSim <- raster::as.matrix(abundSim)
	
	# idsSimRast <- idsSimRast[nrows:1, ]
	# abundSim <- abundSim[nrows:1, ]
	
	# idsSimRast <- raster::raster(idsSimRast, template=sim)
	# abundSim <- raster::raster(abundSim, template=sim)
	
	# names(idsSimRast) <- 'id'
	# names(abundSim) <- 'abundance'
	
	# # cell numbers: renumber so bottom left is (1, 1), increments to the right, then wraps around to next row up
	# refugeCellIds <- simFrame$cellNumSim[which(!is.na(simFrame$refugeAbund))] # gets "raster" cell numbers
	
	# # mean abundance in all refugia
	# meanRefugeAbund <- raster::cellStats(abundSim, 'mean')

	# rescale
	origRefugia <- abund >= threshold
	origRefugiaAbund <- abund * origRefugia
	simAbund <- raster::resample(origRefugiaAbund, sim)
	
	simAbund <- raster::calc(simAbund, fun=function(x) ifelse(x > 1, 1, x))
	simAbund <- raster::calc(simAbund, fun=function(x) ifelse(x < 0, 0, x))
	
	# identify refugia
	simRefugia <- simAbund >= threshold
	simRefugiaId <- raster::clump(simRefugia, directions=8, gaps=FALSE)
	names(simRefugiaId) <- 'refugiaId'
	
	# calculate abundance in refugia
	simAbund <- simAbund * simRefugia
	names(simAbund) <- 'refugiaAbund'
	
	# refuge cell IDs
	# cell in lower left is 1, to its right is 2, etc, then wraps to next row so cell [nrow - 1, 1] is next in line
	nrows <- nrow(sim)
	ncols <- ncol(sim)
	ncells <- raster::ncell(sim)

	v <- rep(seq(nrows * (ncols - 1) - 1, 1, by=-ncols), each=ncols) + 0:(ncols - 1)
	cellNumSim <- matrix(v, nrow=nrows, ncol=ncols, byrow=TRUE)
	cellNumSim <- raster::raster(cellNumSim, template=sim)
	cellNumSim <- as.vector(cellNumSim)
	simRefugiaBinary <- as.vector(simRefugia)
	refugeCellNum <- cellNumSim[simRefugiaBinary]
	if (any(is.na(refugeCellNum))) refugeCellNum <- refugeCellNum[!is.na(refugeCellNum)]
	
	# mean refuge abundance
	meanRefugeAbund <- raster::cellStats(simAbund, 'sum') / length(refugeCellNum)
	
	out <- list(
		simulationScale = raster::stack(simRefugiaId, simAbund),
		refugeCellNum = refugeCellNum,
		meanRefugeAbund = meanRefugeAbund
	)
	
	out
	
}
