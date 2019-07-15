#' draw parameter values from a control file
#'
#' takes an input file with some fixed variables and some prior distributions, outputs a list of parameter values for the simulation
#'
#' @export

#Draw parameters for the simulation... output is parms object, should have a LOT of the stuff in it
#Want to use something similar to what you have for round gobies
drawParms <- function(control = NULL) {
	priors <- read.csv(control, header = TRUE, row.names = 1)
	parms <- c()
	for(p in 1:length(priors[,1])) {
		#print(rownames(priors)[p])
		if(priors$type[p] == "uniform") {
			min <- as.numeric(as.character(priors$min[p]))
			max <- as.numeric(as.character(priors$max[p]))
			parms[p] <- runif(1, min, max)
			rm(min)
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
