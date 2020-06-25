#To Do List...
#1. set population size to 1 at the present for all populations that have gone extinct but need to be kept in the simulation
#2. set migration rates to 0 for all extinct populations
#3. identify the time of extinction for all zero_pops
#4. keep a "base_tmat" object that gets updated with populations that go extinct, but still has migration rates associated with zero pops
#5. add historical events at the time of extinction that expand the population size
#6. substitute row for base_tmat into the appropriate tmat at the time of extinction

#' run fastsimcoal
#'
#' takes a population history and creates sequences based on the coalscent
#'
#' @export
runFSC_step_agg3 = function(
  ph = ph,				#A new pophist object - (pophist, Nvecs, tmat, struct, hab_suit, coalhist)
  l = landscape, 			#A new landscape object - (details, occupied, empty, sampled, hab_suit, sumrast, samplocsrast, samplocs)
  sample_n = NULL,		#Number of sampled individuals per population
  preLGMparms = NULL,		#This has parms for the refuge, preLGM size and timing
  label = "test",			#Label for FSC simulation files
  delete_files = TRUE,	#Logical - clear out .par, .arp, and other FSC outputs?
  num_cores = 1,			#Number of processors to use for FSC
  exec = "fsc25",			#Executable for FSC (needs to be in $PATH variable)
  loc_parms = NULL,		#Vector of locus parameters
  found_Ne = NULL,			#Founding population size, required for STEP change model		
  gmap = NULL,         #gmap relating spatially aggregated cells for coalescent sim to fine-grained cells from forward sim
  MAF = NULL,           #MAF - minor allele frequency filter... loci with minor allele frequencies below this value are excluded from output
  maxloc = 200000       #max number of marker loci to attempt in fastsimcoal
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
  #migmat[[1]] <- round(tmat,5)
  migmat[[1]] <- tmat
  diag(migmat[[1]]) <- 0
  
  allevents <- data.frame(time=NA, source = NA, sink = NA, prop_m = NA, new_size = NA, new_growth=NA, new_mig = NA)
  
  ############################################################
  
  
  ############################################################
  #Deal with empty populations (never colonized,
  # outside of suitable area at all times in the sim, 
  # or below threshold suitability)
  
  #Remove them from the migration matrix
  #empty_pops <- sort(unique(c(l$empty, which(rowSums(simhist) == 0))))
  #empty_pops <- l$empty
  #empty_pops <- sort(unique(c(l$empty, pops$pop[!pops$pop %in% coalhist$src & !pops$pop %in% coalhist$snk & !pops$pop %in% l$empty])))
  #empty_pops2 <- sort(unique(c(l$empty, which(rowSums(simhist) == 0))))    #Don't allow migration with populations that went extinct!
  #New Definition of empty populations... JDR 4/1/2020
  empty <- rep(FALSE, nrow(gmap))
  empty[l$empty] <- TRUE
  empty_gpops <- sort(unique(c(                                    #This was causing a problem with multi-cell refugia - JDR 6/12/2020
    gmap$gpop[!gmap$gpop %in% gmap$gpop[empty == FALSE]],
    which(rowSums(ph$Nvecs) == 0),
    gmap$gpop[!gmap$gpop %in% coalhist$src & !gmap$gpop %in% coalhist$snk])))

  #empty_gpoops <- sort(unique(c(
  #  gmap$gpop[!gmap$gpop %in% gmap$gpop[empty == FALSE]],
  #  which(rowSums(ph$Nvecs) == 0))))
  
  if(length(empty_gpops) > 0) {
    simhist <- simhist[-empty_gpops,]
    #migmat[[1]][empty_gpops,] <- 0		#Turns migration off for all empty pops, and any pops that went extinct (e.g., GoM populations)
    #migmat[[1]][,empty_gpops] <- 0
    migmat[[1]] <- migmat[[1]][-empty_gpops,]
    migmat[[1]] <- migmat[[1]][,-empty_gpops]
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
  #h <- struct["xdim"]*struct["ydim"]   
  h <- max(gmap$gpop)
  
  
  pop_size <- simhist[,length(simhist[1,])]
  #pop_size[pop_size == 0] <- 500			#!# What should we do here?  set to 1?  Set to K?  Set to max attained in this pop?  Can't be 0!!
  extinct_pops <- which(pop_size == 0)
  base_tmat <- migmat[[1]]
  if(length(extinct_pops) > 0) {
    pop_size[extinct_pops] <- 1
    extinctIDs <- extinct_pops - 1
    migmat[[1]][extinct_pops,] <- 0
    migmat[[1]][,extinct_pops] <- 0
    extinctions <- data.frame(pop = extinct_pops, id = extinctIDs, time = NA)
    for(pop in extinctions$pop) {
  	  extinctions$time[extinctions$pop == pop] <- ngens - max(which(simhist[pop,] != 0))
    }
  } else {
    extinctions <- data.frame(pop = NULL, id = NULL, time = NULL)
  }

  #Sample size
  sample_size <- rep(0,h)
  #Growth rate = 0
  growth_rate <- rep(0,h)
  
  if(length(empty_gpops > 0)) {
    sample_size <- sample_size[-empty_gpops]
    growth_rate <- growth_rate[-empty_gpops]
  }
  
  ##!## Slight change here... map first to the gpop column of the gmap object, then to the FSC id (old, new) - JDR 4/1/2020
  sample_ids <- plyr::mapvalues(l$sampled, gmap$pop, gmap$gpop, warn_missing=FALSE)
  if(sum(sample_ids %in% empty_gpops) > 0) {
    stop("Some locations where genetic samples have been collected were not colonized during the forward simulation.  Try other parameter values and rerun getpophist!")
  }
  sample_pops <- plyr::mapvalues(sample_ids, oldID, newID, warn_missing=FALSE)
  sample_size[sample_pops+1] <- sample_n
 
 ###!### THIS STUFF IS NEW!! 
  demelist <- vector("list",length(sample_size))
  for(p in 1:length(sample_size)) {
    demelist[[p]] <- fscDeme(deme.size = pop_size[p], sample.size = sample_size[p], sample.time = 0, inbreeding = 0, growth = growth_rate[p])
  }

  demename <- c()
  for(popu in 1:length(demelist)) {
    if(!newID[popu] %in% sample_pops) {
      demename[popu] <- paste0("unsamp-",popu)
    } else {
      demename[popu] <- l$sampdf$abbrev[l$sampdf$cell %in% gmap$pop[gmap$gpop == oldID[popu]]]
    }
  }
  names(demelist) <- demename

  demes <- do.call(fscSettingsDemes, c(demelist, ploidy = 2))
#  pop_info <- fscPopInfo(pop.size = pop_size, sample.size = sample_size, growth.rate = growth_rate)	
  
  
  ##################################################
  # PARAMETERS FOR SIMULATED LOCI
  ##################################################
  #if(loc_parms$marker == "snp") {
  #  locus_params <- fscLocusParams(locus.type = as.character(loc_parms$marker), num.loci = loc_parms$nloci, mut.rate = loc_parms$seq_length*loc_parms$mu)
  #  attr(locus_params, "opts") <- "-s 0"
  #} else if(loc_parms$marker == "dna") {
  #  locus_params <- fscLocusParams(locus.type = as.character(loc_parms$marker), sequence.length = loc_parms$seq_length, num.chrom = loc_parms$nloci, mut.rate = loc_parms$mu)
  #} else if(loc_parms$marker == "msat") {
  #  locus_params <- fscLocusParams(locus.type = as.character(loc_parms$marker), num.loci = 1, num.chrom = loc_parms$nloci, mut.rate = loc_parms$mu)
  #}
  
  ###!### This stuff is new!!
  num_markers_sim <- loc_parms$nloci*1.25
  if(loc_parms$marker == "snp") {
    genetics <- fscSettingsGenetics(fscBlock_snp(sequence.length = loc_parms$seq_length, mut.rate = loc_parms$mu), num.chrom = num_markers_sim)  
  } else {
    stop("runFSC_step is only built for SNP data (at the moment)")
  }

  ##################################################	
  # HISTORY OF DEMOGRAPHIC EVENTS
  ##################################################
  usemigmat <- 0
  ev <- 1
  gonepops <- c()

  for(G in 0:ngens) {
    tmp_pops <- coalhist[which(coalhist$time == G),]
    #tmp_pops <- FSCpops[which(FSCpops$arrive == G),]
    tmp_leaving <- c()
    
    if(length(tmp_pops[,1]) > 0) {
      for(fp in 1:length(tmp_pops[,1])) {
        if(!is.na(tmp_pops$snk[fp])) {
          ##!## Adding this if statement to check whether this is the ORIGINAL founding event - JDR 4/1/2020
          if(tmp_pops$prop[fp] == 1) {
            allevents[ev,] <- c(G-1, tmp_pops$src[fp], tmp_pops$src[fp], 0, as.numeric(found_Ne)/as.numeric(pop_size[tmp_pops$src[fp]+1]), 0, usemigmat)
            ev <- ev+1
            tmp_leaving <- c(tmp_leaving, tmp_pops$src[fp]+1)
          }
          ##!## Adding in tmp_pops$prop for migrant proportion - JDR 4/1/2020
          allevents[ev,] <- c(G, tmp_pops$src[fp], tmp_pops$snk[fp], tmp_pops$prop[fp], 1, 0, usemigmat+1)
          ev <- ev+1
        }
      }
  
	    gonepops <- c(gonepops, tmp_leaving)    
      migmat[[usemigmat+2]] <- migmat[[usemigmat+1]]
      if(G %in% extinctions$time) {
  	  	for(pop in extinctions$pop[extinctions$time == G]) {
  	  		subrow <- base_tmat[pop,]
  	  		subcol <- base_tmat[,pop]
  	  		subrow[gonepops] <- 0
  	  		subcol[gonepops] <- 0
  	  		migmat[[usemigmat+2]][pop,] <- subrow
  	  		migmat[[usemigmat+2]][,pop] <- subcol
  	  		allevents[ev,] <- c(G, extinctions$id[extinctions$pop == pop], extinctions$id[extinctions$pop == pop], 0, max(simhist[pop,])/pop_size[pop], 0, usemigmat+1)
  	  		ev <- ev+1
  	  		rm(subrow, subcol)
  	  	}
  	  }
      migmat[[usemigmat+2]][tmp_leaving,] <- 0
      migmat[[usemigmat+2]][,tmp_leaving] <- 0
      usemigmat <- usemigmat + 1
      
      rm(tmp_pops)
      rm(tmp_leaving)
      
    }
  }

  migmat[[length(migmat)+1]] <- migmat[[length(migmat)]]
  migmat[[length(migmat)]][migmat[[length(migmat)]] > 0] <- 0
  nomig_migmat <- length(migmat)-1

  migration <- do.call(fscSettingsMigration, migmat)

  ##################################################
  # ADD HISTORICAL EVENTS TO RESURRECT EXTINCT POPS 
  ##################################################



  ##################################################
  # COMBINE EVERYTHING, SORT, AND RUN FSC 
  ##################################################
  allevents <- allevents[order(allevents$time, allevents$new_mig, allevents$prop_m, allevents$source),]	
  
  ###!### New stuff here!!
  eventlist <- vector("list", nrow(allevents))
  for(ev in 1:nrow(allevents)) {
    eventlist[[ev]] <- fscEvent(event.time = allevents$time[ev], source = allevents$source[ev], sink = allevents$sink[ev], prop.migrants = allevents$prop_m[ev], new.size = allevents$new_size[ev], new.growth = allevents$new_growth[ev], migr.mat = allevents$new_mig[ev])
  }
  #eventlist[[length(eventlist)+1]] <- fscEvent(event.time = ngens, source = newID[which(simhist[,1] > 0)], sink = newID[which(simhist[,1] > 0)], prop.migrants = 0, new.size = 1, new.growth = 0, migr.mat = nomig_migmat)   #Commenting this out, unnecessary, allow migration among refugial cells - JDR 6/12/2020
  #hist_ev <- fscHistEv(num.gen = allevents$time, source.deme = allevents$source, sink.deme = allevents$sink, prop.migrants = allevents$prop_m, new.sink.size = allevents$new_size, new.sink.growth = allevents$new_growth, new.mig.mat = allevents$new_mig)
  #hist_ev1 <- fscHistEv(num.gen = ngens, source.deme = newID[which(simhist[,1] > 0)], sink.deme = newID[which(simhist[,1] > 0)], prop.migrants = 0, new.sink.size = 1, new.sink.growth = 0, new.mig.mat = "nomig")
  
  if(length(newID[which(simhist[,1] > 0)]) == 1) {
    eventlist[[length(eventlist)+1]] <- fscEvent(event.time = preLGMparms$preLGM_t[1], source = newID[which(simhist[,1] > 0)], sink = newID[which(simhist[,1] > 0)], prop.migrants = 0, new.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.growth = 0, migr.mat = nomig_migmat)
    #hist_ev2 <- fscHistEv(num.gen = preLGMparms$preLGM_t[1], newID[which(simhist[,1] > 0)], sink.deme = newID[which(simhist[,1] > 0)], prop.migrants = 0, new.sink.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.sink.growth = 0, new.mig.mat = "nomig")
    #hist_ev <- rbind(hist_ev, hist_ev1, hist_ev2)
  } else {
    arb.preLGM.pop <- newID[which(simhist[,1] > 0)][1]
    #for(rem in 2:(length(which(simhist[,1] > 0))-1) {
    for(rem in 2:length(which(simhist[,1] > 0))) {
      eventlist[[length(eventlist)+1]] <- fscEvent(event.time = preLGMparms$preLGM_t[1]-1, source = newID[which(simhist[,1] >0)][rem], sink = arb.preLGM.pop, prop.migrants = 1, new.size = 1, new.growth = 0, migr.mat = nomig_migmat)
    }
    eventlist[[length(eventlist)+1]] <- fscEvent(event.time = preLGMparms$preLGM_t[1], source = arb.preLGM.pop, sink = arb.preLGM.pop, prop.migrants = 0, new.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.growth = 0, migr.mat = nomig_migmat)

    #hist_ev2 <- fscHistEv(num.gen = preLGMparms$preLGM_t[1]-1, source.deme = newID[which(simhist[,1] > 0)][-1], sink.deme = arb.preLGM.pop, prop.migrants = 1, new.sink.size = 1, new.sink.growth = 0, new.mig.mat = "nomig")
    #hist_ev3 <- fscHistEv(num.gen = preLGMparms$preLGM_t[1], source.deme = arb.preLGM.pop, sink.deme = arb.preLGM.pop, prop.migrants = 0, new.sink.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.sink.growth = 0, new.mig.mat = "nomig")
    #hist_ev <- rbind(hist_ev, hist_ev1, hist_ev2, hist_ev3)
  }
  events <- do.call(fscSettingsEvents, eventlist)

  #fscout <- fastsimcoal(label = label, pop.info = pop_info, locus.params = locus_params, mig.rates = migmat, hist.ev = hist_ev, num.cores = num_cores, delete.files = delete_files, exec = exec)

  ###!### This is new!!
  varSNPs <- 0
  while(varSNPs < loc_parms$nloci) {
    p <- fscWrite(demes = demes, genetics = genetics, events = events, migration = migration, label = label, use.wd = TRUE)
    p <- fscRun(p, all.sites = FALSE, inf.sites = FALSE, no.arl.output = FALSE, dna.to.snp = TRUE, quiet = FALSE, num.cores = num_cores, exec = exec)
    out <- fscReadArp(p)
    message(paste("Coalescent simulation with",attr(genetics,"num.chrom"), "loci resulted in", (ncol(out)-2)/2, "polymorphic sites"), appendLF=TRUE)
    
    fscout <- sampleOnePerLocus(mat = out, MAF = MAF)
    tmp_gtype <- df2gtypes(fscout, ploidy = 2)
    varSNPs <- sum(numAlleles(tmp_gtype)$num.alleles == 2)
    message(paste("Subsampling to one SNP per locus...", varSNPs, "loci have at least one polymorphic sites with MAF >", MAF), appendLF=TRUE)

    if(varSNPs < loc_parms$nloci) {
      if(attr(genetics, "num.chrom") == maxloc) {
        stop("Too few SNPs pass the MAF!  Redrawing parameter values for this replicate!  Setting maxSNPstried to a larger value will lead to longer simulations")
      }
      scaleSNP <- 1.1*(attr(genetics, "num.chrom")/varSNPs)
      newSNPnum <- round(scaleSNP*attr(genetics, "num.chrom"),0)
      if(newSNPnum > maxloc) {
        newSNPnum <- maxloc
        message(paste("Coalescent simulation with", attr(genetics, "num.chrom"), "loci resulted in", varSNPs, "variable markers. Trying again with the maximum number of loci allowable -",newSNPnum, "!"), appendLF=TRUE)
      } else {
        message(paste("Coalescent simulation with", attr(genetics, "num.chrom"), "loci resulted in", varSNPs, "variable markers. Trying again with",newSNPnum, "loci!"), appendLF=TRUE)
      }
      attr(genetics, "num.chrom") <- newSNPnum
      rm(scaleSNP, newSNPnum)
    } else if(varSNPs > loc_parms$nloci) {
      #SAMPLE loc_parms$nloci markers here!!!
      #####
      nalleles <- numAlleles(tmp_gtype)
      varchromnames <-  nalleles$locus[nalleles$num.alleles == 2]
      keep.loci <- sample(varchromnames, loc_parms$nloci, replace = FALSE)
      keep.columns <- c(1,2)
      for(locus in keep.loci) {
        keep.columns <- c(keep.columns, grep(locus, colnames(fscout)))
      }
      fscout <- fscout[,keep.columns]
    }
  }
  
  #Need to get the names working as we want them to work...
  popid_inds <- sapply(strsplit(out$id, split = "_"), FUN = function(x) as.numeric(x[1]))
  fscout$deme <- demes$deme.name[popid_inds]
  
  ###!### JDR 3/25/2020
  #I don't think we really want this formatted as a gtypes object anymore :(
  ###!###

  #Naming strata prior to output - this is no longer needed with new strataG!!  - JDR 4/3/2020
  #FSCid <- sapply(strsplit(fscout@data$ids,"_"), function(x){x[1]})
  #FSCgridid <- plyr::mapvalues(FSCid,newID,oldID,warn_missing=FALSE)
  #Then map these old ID's back to the abbreviation of the sampled population
  #FSCabbrev <- plyr::mapvalues(FSCgridid,l$sampdf$cell,l$sampdf$abbrev)
  #fscout@data$strata <- FSCabbrev
  #fscout@data$ids <- paste0(fscout@data$strata, "_", c(1:length(fscout@data[,1])))
  fscout <- fscout[order(fscout$deme),]
  
  message(paste("Coalescent simulation complete -", (ncol(fscout)-2)/2, "SNPs generated"), appendLF=TRUE)
  fscout
  
  
}
