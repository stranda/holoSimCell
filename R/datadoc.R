#' Description of alternative ENMs
#' 
#' A two-element list that describes 24 different ENMs used in simulations from Castilla et al. (2023). 
#' 
#' \itemize{
#'	\code{enms} is a data frame with the following columns:
#' 	\item numRefugia. The number of spatially separated refugia at the LGM
#' 	\item gcm. The Global Circulation Model used in the ENM
#' 	\item ext. The buffer around occurrence points used to define background regions (in km)
#' 	\item algo. The niche modeling algorithm used in the model
#' 	\item threshold. The threshold used to define refuge cells, cells with suitability above this value at the LGM are occupied at the start of a simulation
#'	\item cbi. Continuous Boyce Index - a measure of ENM performance
#'	\item tss. True Skill Statistic - a measure of ENM performance
#'	\item se. ?????
#'	\item sp. ?????
#'	\item auc. The area under the receiver-operator curve - a measure of ENM performance
#'	\item meanRefugeAbund. Average habitat suitability of cells classified as refugia (given the defined threshold)
#'	\item rasterStackName. A text string name for the ENM, describes GCM, extent, and algorithm used. Refers to landscapes located in \code{system.file(package="holoSimCell","extdata/landscapes")}
#' }
#' 
#' \itemize{
#'	\code{refugeCellNum} is a list with 24 elements, 1 per ENM, with each element storing a vector of IDs for refuge cells.  These ids correspond to the rasters named in 'rasterStackName' in the \code{enms} object
#' }
#' @docType data
#' @keywords datasets
#' @name enmScenarios
#' @usage data(enmScenarios)
#' @format two element list (data frame and 24 element list)
NULL



#' Genotypes for green ash (\emph{Fraxinus pennsylvanica}) 
#' 
#' A set of SNP genotypes collected for green ash (\emph{Fraxinus pennsylvanica} Marshall) throughout North america.  This data frame has columns that correspond to individuals and rows that correspond to SNP loci.  The individual IDs are used in column names and mapped back to populations in \code{data(popmap)}. Genotype calls were made following alignment to the \emph{F. pennsylvanica} genome from Huff et al. 2021 - \url{https://doi.org/10.5281/zenodo.5176117}.
#' 
#' @references \url{https://doi.org/10.5281/zenodo.5176117}
#' @docType data
#' @keywords datasets
#' @name imputed
#' @usage data(imputed)
#' @format data frame with 1306 rows (loci) and 362 columns (individuals)
NULL



#' Prediction surfaces for green ash pollen
#' 
#' A list of 200 pointers to files containing predictions of green ash abundance through time based on a Bayesian pollen vegetation model.  Each represents a pull from the posterior distribution of green ash abundances estimated from the model.
#' 
#' \itemize{
#'	Each element is a list of two elements:	
#' 	\item file. A text \code{.rda} file that stores the landscape object associated with this draw from the posterior.
#' 	\item refs. A vector of refuge cell IDs.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name pollenPulls
#' @usage data(pollenPulls)
#' @format A list of 200 elements, each a list of 2 elements (character string and vector)
NULL


#' Association between green ash IDs and populations
#' 
#' This data frame provides the information to connect individual IDs (\code{id}; also found in the \code{imputed} data frame) to population IDs (\code{abbrev}) and their geographic coordinates found in \code{pts}.
#' 
#' @docType data
#' @keywords datasets
#' @name popmap
#' @usage data(popmap)
#' @format data frame with 2 columns and 761 rows
NULL


#' Geographic coordinates of green ash populations
#' 
#' A data frame that includes columns that link population abbreviations found in \code{popmap} to the full population names and their coordinates in WGS84 Latitude and Longitude.  Not all populations included here are found in the imputed dataset due to variability in data quality.
#' 
#' \itemize{
#' 	\item pop. Full population name
#' 	\item abbrev. short population name
#' 	\item lat. Latitude in decimal degrees
#' 	\item long. Longitude in decimal degrees
#' 	\item N. sample size collected (this may overstate the number of individuals compared to \code{imputed} due to variable data quality.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name pts
#' @usage data(pts)
#' @format data frame with 5 columns and 49 rows
NULL


