#' run fastsimcoal
#'
#' takes a population history and creates sequences based on the coalscent
#'
#' @export
runFSC_step = function(
  ph = ph,				#A new pophist object - (pophist, Nvecs, tmat, struct, hab_suit, coalhist)
  l = landscape, 			#A new landscape object - (details, occupied, empty, sampled, hab_suit, sumrast, samplocsrast, samplocs)
  sample_n = NULL,		#Number of sampled individuals per population
  preLGMparms = NULL,		#This has parms for the refuge, preLGM size and timing
  label = "test",			#Label for FSC simulation files
  delete_files = TRUE,	#Logical - clear out .par, .arp, and other FSC outputs?
  num_cores = 1,			#Number of processors to use for FSC
  exec = "fsc25",			#Executable for FSC (needs to be in $PATH variable)
  loc_parms = NULL,		#Vector of locus parameters
  found_Ne = NULL			#Founding population size, required for STEP change model		
) {
  
  ############################################################	
  #Error checks
  if(is.null(found_Ne)) {
    stop("Must specify founding Ne if using the STEP change growth model")
  }
  ############################################################
  
  
  ############################################################
  #Create objects used downstream as FSC inputs
  #Initialize objects for FSC
  pops <- ph$pophist
  simhist <- ph$Nvecs
  tmat <- ph$tmat
  struct <- ph$struct
  coalhist <- ph$coalhist
  
  migmat <- vector("list", 1)
  migmat[[1]] <- round(tmat,5)
  diag(migmat[[1]]) <- 0
  
  allevents <- data.frame(time=NA, source = NA, sink = NA, prop_m = NA, new_size = NA, new_growth=NA, new_mig = NA)
  
  ############################################################
  
  
  ############################################################
  #Deal with empty populations (never colonized,
  # outside of suitable area at all times in the sim, 
  # or below threshold suitability)
  
  #Remove them from the migration matrix
  #empty_pops <- sort(unique(c(l$empty, which(rowSums(simhist) == 0))))
  empty_pops <- l$empty
  empty_pops2 <- sort(unique(c(l$empty, which(rowSums(simhist) == 0))))    #Don't allow migration with populations that went extinct!
  if(length(empty_pops2) > 0) {
    simhist <- simhist[-empty_pops,]
    migmat[[1]][empty_pops2,] <- 0		#Turns migration off for all empty pops, and any pops that went extinct (e.g., GoM populations)
    migmat[[1]][,empty_pops2] <- 0
    migmat[[1]] <- migmat[[1]][-empty_pops,]
    migmat[[1]] <- migmat[[1]][,-empty_pops]
  } 
  ############################################################
  
  
  ############################################################
  #Build objects with IDs fixed to account for empty populations and the different starting index (0 for FSC, 1 for forward sim)
  #Correct arrival times so you're thinking backwards in time
  #Using the new coalhist and landscape objects
  ngens <- struct["maxtime"]
  oldID <- sort(unique(c(coalhist$src, which(ph$Nvecs[,1] > 0))))
  newID <- c(1:length(oldID))-1
  coalhist$snk <- plyr::mapvalues(coalhist$snk, oldID, newID, warn_missing=FALSE)
  coalhist$src <- plyr::mapvalues(coalhist$src, oldID, newID, warn_missing=FALSE)
  coalhist$time <- ngens - coalhist$time
  ############################################################
  
  
  ##################################################
  # POPULATION INFORMATION
  ##################################################
  #Definitions for pop.info
  h <- struct["xdim"]*struct["ydim"]   
  
  pop_size <- simhist[,length(simhist[1,])]
  pop_size[pop_size == 0] <- 500			#!# What should we do here?  set to 1?  Set to K?  Set to max attained in this pop?  Can't be 0!!
  
  #Sample size
  sample_size <- rep(0,h)
  #Growth rate = 0
  growth_rate <- rep(0,h)
  
  if(length(empty_pops > 0)) {
    sample_size <- sample_size[-empty_pops]
    growth_rate <- growth_rate[-empty_pops]
  }
  
  sample_pops <- plyr::mapvalues(l$sampled, oldID, newID, warn_missing=FALSE)
  sample_size[sample_pops] <- sample_n
  
  pop_info <- fscPopInfo(pop.size = pop_size, sample.size = sample_size, growth.rate = growth_rate)	
  
  
  ##################################################
  # PARAMETERS FOR SIMULATED LOCI
  ##################################################
  if(loc_parms$marker == "snp") {
    locus_params <- fscLocusParams(locus.type = as.character(loc_parms$marker), num.loci = loc_parms$nloci, mut.rate = loc_parms$seq_length*loc_parms$mu)
    attr(locus_params, "opts") <- "-s 0"
  } else if(loc_parms$marker == "dna") {
    locus_params <- fscLocusParams(locus.type = as.character(loc_parms$marker), sequence.length = loc_parms$seq_length, num.chrom = loc_parms$nloci, mut.rate = loc_parms$mu)
  } else if(loc_parms$marker == "msat") {
    locus_params <- fscLocusParams(locus.type = as.character(loc_parms$marker), num.loci = 1, num.chrom = loc_parms$nloci, mut.rate = loc_parms$mu)
  }
  
  
  
  ##################################################	
  # HISTORY OF DEMOGRAPHIC EVENTS
  ##################################################
  usemigmat <- 0
  ev <- 1
  
  for(G in 0:ngens) {
    tmp_pops <- coalhist[which(coalhist$time == G),]
    #tmp_pops <- FSCpops[which(FSCpops$arrive == G),]
    tmp_leaving <- c()
    
    if(length(tmp_pops[,1]) > 0) {
      for(fp in 1:length(tmp_pops[,1])) {
        if(!is.na(tmp_pops$snk[fp])) {
          allevents[ev,] <- c(G-1, tmp_pops$src[fp], tmp_pops$src[fp], 0, as.numeric(found_Ne)/as.numeric(pop_size[tmp_pops$src[fp]+1]), 0, usemigmat)
          ev <- ev+1
          allevents[ev,] <- c(G, tmp_pops$src[fp], tmp_pops$snk[fp], 1, 1, 0, usemigmat+1)
          ev <- ev+1
          tmp_leaving <- c(tmp_leaving, tmp_pops$src[fp]+1)
        }
      }
      
      migmat[[usemigmat+2]] <- migmat[[usemigmat+1]]
      migmat[[usemigmat+2]][tmp_leaving,] <- 0
      migmat[[usemigmat+2]][,tmp_leaving] <- 0
      usemigmat <- usemigmat + 1
      
      rm(tmp_pops)
      rm(tmp_leaving)
      
    }
  }
  
  ##################################################
  # COMBINE EVERYTHING, SORT, AND RUN FSC 
  ##################################################
  allevents <- allevents[order(allevents$time, allevents$new_mig, allevents$prop_m, allevents$source),]	
  
  hist_ev <- fscHistEv(num.gen = allevents$time, source.deme = allevents$source, sink.deme = allevents$sink, prop.migrants = allevents$prop_m, new.sink.size = allevents$new_size, new.sink.growth = allevents$new_growth, new.mig.mat = allevents$new_mig)
  hist_ev1 <- fscHistEv(num.gen = ngens, source.deme = newID[which(simhist[,1] > 0)], sink.deme = newID[which(simhist[,1] > 0)], prop.migrants = 0, new.sink.size = 1, new.sink.growth = 0, new.mig.mat = "nomig")
  
  if(length(newID[which(simhist[,1] > 0)]) == 1) {
    hist_ev2 <- fscHistEv(num.gen = preLGMparms$preLGM_t[1], newID[which(simhist[,1] > 0)], sink.deme = newID[which(simhist[,1] > 0)], prop.migrants = 0, new.sink.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.sink.growth = 0, new.mig.mat = "nomig")
    hist_ev <- rbind(hist_ev, hist_ev1, hist_ev2)
  } else {
    arb.preLGM.pop <- newID[which(simhist[,1] > 0)][1]
    hist_ev2 <- fscHistEv(num.gen = preLGMparms$preLGM_t[1]-1, source.deme = newID[which(simhist[,1] > 0)][-1], sink.deme = arb.preLGM.pop, prop.migrants = 1, new.sink.size = 1, new.sink.growth = 0, new.mig.mat = "nomig")
    hist_ev3 <- fscHistEv(num.gen = preLGMparms$preLGM_t[1], source.deme = arb.preLGM.pop, sink.deme = arb.preLGM.pop, prop.migrants = 0, new.sink.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.sink.growth = 0, new.mig.mat = "nomig")
    hist_ev <- rbind(hist_ev, hist_ev1, hist_ev2, hist_ev3)
  }
  
  fscout <- fastsimcoal(label = label, pop.info = pop_info, locus.params = locus_params, mig.rates = migmat, hist.ev = hist_ev, num.cores = num_cores, delete.files = delete_files, exec = exec)
  
  fscout
  
}

#' run fastsimcoal
#'
#' takes a population history and creates sequences based on the coalscent
#'
#' @export
runFSC_step_OLD = function(
	popsobj = NULL,  		#Pass the 'pops' object from getpophist.cells()
	preLGMparms = NULL,		#This has parms for the refuge, preLGM size and timing
	sample_pops = NULL,		#The IDs of sampled populations
	sample_n = NULL,		#The sample size per population
	label = "test",			#Label for FSC simulation files
	delete_files = TRUE,	#Logical - clear out .par, .arp, and other FSC outputs?
	num_cores = 1,			#Number of processors to use for FSC
	exec = "fsc25",			#Executable for FSC (needs to be in $PATH variable)
	loc_parms = NULL,		#Vector of locus parameters
	found_Ne = NULL			#Founding population size, required for STEP change model		
	) {

	#require(strataG)

############################################################	
#Error checks
	if(is.null(found_Ne)) {
		stop("Must specify founding Ne if using the STEP change growth model")
	}
############################################################


############################################################
#Create objects used downstream as FSC inputs
	#Initialize objects for FSC
	pops <- popsobj$pophist
	simhist <- popsobj$Nvecs
	tmat <- popsobj$tmat
	struct <- popsobj$struct

	migmat <- vector("list", 1)
	migmat[[1]] <- round(tmat,5)
	diag(migmat[[1]]) <- 0

	allevents <- data.frame(time=NA, source = NA, sink = NA, prop_m = NA, new_size = NA, new_growth=NA, new_mig = NA)

############################################################


############################################################
#Deal with empty populations (never colonized or went extinct before the end of the sim)
	#Look for empty popualtions
	empty_pops <- pops$pop[is.na(pops$arrive) | simhist[,length(simhist[1,])] == 0] 

	#Remove empty populations from all FSC inputs
	if(length(empty_pops) > 0) {
		simhist <- simhist[-empty_pops,]
		FSCpops <- pops[-empty_pops,]
		migmat[[1]] <- migmat[[1]][-empty_pops,]
		migmat[[1]] <- migmat[[1]][,-empty_pops]
	} else {
		FSCpops <- pops
	}
############################################################


############################################################
#Build objects with IDs fixed to account for empty populations and the different starting index (0 for FSC, 1 for forward sim)
#Correct arrival times so you're thinking backwards in time
	ngens <- struct["maxtime"]
	FSCpops$pop <- c(1:length(FSCpops$pop))-1
	FSCpops$arrive <- ngens - FSCpops$arrive
	for(r in 1:length(FSCpops$arrive)) {
		popsrc <- FSCpops$source[r]
	
	#This fixes the IDs of the source pops (accounting for empty_pops and the different index starting point 0v1)
		if(!is.na(popsrc)) {
			srcrow <- pops$row[pops$pop == popsrc]
			srccol <- pops$col[pops$pop == popsrc]
			FSCpops$source[r] <- FSCpops$pop[FSCpops$row == srcrow & FSCpops$col == srccol]
			rm(popsrc, srcrow, srccol)
		}
	} 
############################################################


##################################################
# POPULATION INFORMATION
##################################################
#Definitions for pop.info
	h <- struct["xdim"]*struct["ydim"]   
	
	pop_size <- simhist[,length(simhist[1,])]
	
	#Sample size
	sample_size <- rep(0,h)
	sample_size[sample_pops] <- sample_n
	
	#Growth rate = 0
	growth_rate <- rep(0,h)

	if(length(empty_pops > 0)) {
	  sample_size <- sample_size[-empty_pops]
	  growth_rate <- growth_rate[-empty_pops]
	}

	pop_info <- fscPopInfo(pop.size = pop_size, sample.size = sample_size, growth.rate = growth_rate)	


##################################################
# PARAMETERS FOR SIMULATED LOCI
##################################################
	if(loc_parms$marker == "snp") {
		locus_params <- fscLocusParams(locus.type = as.character(loc_parms$marker), num.loci = loc_parms$nloci, mut.rate = loc_parms$seq_length*loc_parms$mu)
	} else if(loc_parms$marker == "dna") {
		locus_params <- fscLocusParams(locus.type = as.character(loc_parms$marker), sequence.length = loc_parms$seq_length, num.chrom = loc_parms$nloci, mut.rate = loc_parms$mu)
	} else if(loc_parms$marker == "msat") {
		locus_params <- fscLocusParams(locus.type = as.character(loc_parms$marker), num.loci = 1, num.chrom = loc_parms$nloci, mut.rate = loc_parms$mu)
	}
	
	attr(locus_params, "opts") <- "-s 0"

##################################################	
# HISTORY OF DEMOGRAPHIC EVENTS
##################################################
	usemigmat <- 0
	ev <- 1

	for(G in 0:ngens) {
		tmp_pops <- FSCpops[which(FSCpops$arrive == G),]
		tmp_leaving <- c()

		if(length(tmp_pops[,1]) > 0) {
			for(fp in 1:length(tmp_pops[,1])) {
				if(!is.na(tmp_pops$source[fp])) {
					allevents[ev,] <- c(G-1, tmp_pops$pop[fp], tmp_pops$pop[fp], 0, as.numeric(found_Ne)/as.numeric(pop_size[tmp_pops$pop[fp]+1]), 0, usemigmat)
					ev <- ev+1
					allevents[ev,] <- c(G, tmp_pops$pop[fp], tmp_pops$source[fp], 1, 1, 0, usemigmat+1)
					ev <- ev+1
					tmp_leaving <- c(tmp_leaving, tmp_pops$pop[fp]+1)
				}
			}
			
			migmat[[usemigmat+2]] <- migmat[[usemigmat+1]]
			migmat[[usemigmat+2]][tmp_leaving,] <- 0
			migmat[[usemigmat+2]][,tmp_leaving] <- 0
			usemigmat <- usemigmat + 1

			rm(tmp_pops)
			rm(tmp_leaving)

		}
	}

##################################################
# COMBINE EVERYTHING, SORT, AND RUN FSC 
##################################################
	allevents <- allevents[order(allevents$time, allevents$new_mig, allevents$prop_m, allevents$source),]	

	hist_ev <- fscHistEv(num.gen = allevents$time, source.deme = allevents$source, sink.deme = allevents$sink, prop.migrants = allevents$prop_m, new.sink.size = allevents$new_size, new.sink.growth = allevents$new_growth, new.mig.mat = allevents$new_mig)
	hist_ev1 <- fscHistEv(num.gen = ngens, source.deme = FSCpops$pop[FSCpops$arrive == ngens], sink.deme = FSCpops$pop[FSCpops$arrive == ngens], prop.migrants = 0, new.sink.size = 1, new.sink.growth = 0, new.mig.mat = "nomig")
	
	if(length(FSCpops$pop[FSCpops$arrive == ngens]) == 1) {
		hist_ev2 <- fscHistEv(num.gen = preLGMparms$preLGM_t[1], source.deme = FSCpops$pop[FSCpops$arrive == ngens], sink.deme = FSCpops$pop[FSCpops$arrive == ngens], prop.migrants = 0, new.sink.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.sink.growth = 0, new.mig.mat = "nomig")
		hist_ev <- rbind(hist_ev, hist_ev1, hist_ev2)
	} else {
		arb.preLGM.pop <- FSCpops$pop[FSCpops$arrive == ngens][1]
		hist_ev2 <- fscHistEv(num.gen = preLGMparms$preLGM_t[1]-1, source.deme = FSCpops$pop[FSCpops$arrive == ngens][-1], sink.deme = arb.preLGM.pop, prop.migrants = 1, new.sink.size = 1, new.sink.growth = 0, new.mig.mat = "nomig")
		hist_ev3 <- fscHistEv(num.gen = preLGMparms$preLGM_t[1], source.deme = arb.preLGM.pop, sink.deme = arb.preLGM.pop, prop.migrants = 0, new.sink.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.sink.growth = 0, new.mig.mat = "nomig")
		hist_ev <- rbind(hist_ev, hist_ev1, hist_ev2, hist_ev3)
	}
	
	fscout <- fastsimcoal(label = label, pop.info = pop_info, locus.params = locus_params, mig.rates = migmat, hist.ev = hist_ev, num.cores = num_cores, delete.files = delete_files, exec = exec)

	fscout

}
