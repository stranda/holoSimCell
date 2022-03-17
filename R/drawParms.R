#' Draw simulation parameter values from a control file with specified prior distributions
#'
#' Takes an input file with some fixed variables and some prior distributions, outputs a list of parameter values for the simulation
#'
#' @param control The full path to a .csv formatted file specifying prior distributions for model parameters.  
#' 
#' @details
#' \itemize{ 
#' The control file includes five required columns:
#' \item{\code{param}} {the case-sensitive parameter name.}
#' \item{\code{min}} {the lower bound of uniform and loguniform prior distributions (NA if type = "beta", parameter name if type = "conditional").}
#' \item{\code{max}} {the upper bound of uniform and loguniform prior distributions (scalar for beta distribution if type = "beta", parameter name if type = "conditional").}
#' \item{\code{type}} {a character string specifying the distribution to draw from (options are "fixed", "uniform", "loguniform", "conditional", "logical", and "beta").}
#' \item{\code{is.int}} {a TRUE/FALSE indicating whether the parameter is an integer (if TRUE, randomly-drawn values are rounded to the nearest whole number).}
#' }
#'
#' \itemize{
#' The control file includes the following required parameters as rows: 
#' \item{\code{shortscale}} {the scale parameter for the distribution of short-distance dispersal events.  \emph{HELP HERE!!}}
#' \item{\code{shortshape}} {the shape parameter for the distribution of short-distance dispersal events.  \emph{HELP HERE!!}}
#' \item{\code{sz}} {the cell size in arbitrary units.  \emph{HELP HERE!!}}
#' \item{\code{nloci}} {the number of genetic marker loci included in the observed dataset.}
#' \item{\code{seq_length}} {the sequence length of the marker, used along with mu below in the fastsimcoal DNA model.}
#' \item{\code{mu}} {the mutation rate (per bp, per generation) of the simulated locus, used with seq_length above in the fastsimcoal DNA model.}
#' \item{\code{G}} {the generation time of the species in years.}
#' \item{\code{longmean}} {the average distance of long-distance dispersal events.}
#' \item{\code{lambda}} {the rate of population growth within cells.}
#' \item{\code{mix}} {the mixture parameter for the distribution of dispersal distances, fraction of dispersal events that are long-distance.}
#' \item{\code{Ne}} {the maximum effective population size of a cell in the landscape.}
#' \item{\code{preLGM_t}} {the time at which refugial populations diverged from one another (e.g., the last interglacial).}
#' \item{\code{preLGM_Ne}} {the effective population size of the species prior to refuge divergence.}
#' \item{\code{found_Ne}} {the effective population size of newly colonized populations (used in coalescent simulations only).}
#' \item{\code{ref_Ne}} {the effective population size of cells included in refugia (scaled by cell-specific habitat suitability in the first time step).}
#' }
#'
#' @return One-row, named data frame of parameter values for a simulation
#'
#' @examples
#' library(holoSimCell)
#' drawParms(control = system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell"))
#' drawParms
#'
#' @export

drawParms <- function(control = NULL) {
	priors <- read.csv(control, header = TRUE, row.names = 1)
	parms <- c()
	for(p in 1:length(priors[,1])) {
		if(priors$type[p] == "uniform") {
			min <- as.numeric(as.character(priors$min[p]))
			max <- as.numeric(as.character(priors$max[p]))
			parms[p] <- runif(1, min, max)
			rm(min)
			rm(max)
		} else if(priors$type[p] == "loguniform") {
		  min <- as.numeric(as.character(priors$min[p]))
		  max <- as.numeric(as.character(priors$max[p]))
		  parms[p] <- exp(runif(1, log(min), log(max)))
		  rm(min)
		  rm(max)
		} else if(priors$type[p] == "beta") {
		  max <- as.numeric(as.character(priors$max)[p])
		  parms[p]<- max*rbeta(1,1,4)
		  rm(max)
		} else if(priors$type[p] == "conditional") {
			if(length(which(rownames(priors)==as.character(priors$min[p]))) == 0) {
				min <- as.numeric(as.character(priors$min[p]))
				max <- parms[which(rownames(priors)==priors$max[p])]
				parms[p] <- runif(1, min, max)
				rm(min)
				rm(max)
			} else if(length(which(rownames(priors)==as.character(priors$max[p]))) == 0) {
				min <- parms[which(rownames(priors)==priors$max[p])]
				max <- as.numeric(as.character(priors$min[p]))
				parms[p] <- runif(1, min, max)
				rm(min)
				rm(max)
			}	
		} else if(priors$type[p] == "fixed") {
			parms[p] <- as.numeric(as.character(priors$min[p]))
		} else if(priors$type[p] == "character") {
			parms[p] <- as.character(priors$min[p])
		} else if(priors$type[p] == "logical") {
			parms[p] <- sample(c(0,1), 1, replace = FALSE)
		}	
		
		if(priors$is.int[p] == TRUE) {
			parms[p] = round(parms[p])
		}

	} 
	
	parms <- as.data.frame(t(parms))

	names(parms) <- rownames(priors)

	parms

}
