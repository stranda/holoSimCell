#' Run fastsimcoal 
#'
#' Takes a population history from the forward demographic simulation and parameterizes a coalescent simulation in fastsimcoal (Excoffier et al. 2013)
#'
#' @param ph A pophist object, output by the forward demographic simulation function: \code{getpophist2.cells}.  Includes slots: \code{pophist}, \code{Nvecs}, \code{tmat}, \code{struct}, \code{hab_suit}, \code{coalhist}, \code{popslst}, \code{old_hab_suit}, \code{old_tmat}, \code{gmap}.  Cells in typical landscapes (>400 cells) should be aggregated after the forward-time simulation using the functions \code{make.gmap} and \code{pophist.aggregate} to improve computational efficiency. 
#' @param l An object that defines the simulation landscape from the forward demographic simulation, includes slots: \code{details}, \code{occupied}, \code{empty}, \code{sampled}, \code{hab_suit}, \code{sumrast}, \code{samplocsrast}, \code{samplocs}, \code{sampdf}, \code{NAmask}.
#' @param sample_n The number of sampled (diploid) individuals per population.  Assumed to be equal across all sampled populations.
#' @param preLGMparms A one row data frame of parameters related to the size and history of populations in glacial refugia: \code{preLGM_t} (time in generations at which refugial populations are assumed to diverge, e.g., previous interglacial), \code{preLGM_Ne} (effective population size in the ancestral population), \code{ref_Ne} (effective size of populations in refugial cells).
#' @param label A text string specifying a filename prefix for output from fastsimcoal.  Note that hyphens (-) should be avoided; they are converted to periods (.) in files output by fastsimcoal.  In order to delete files with \code{strataG::fscCleanup}, avoid hyphens in file names.
#' @param delete_files A logical (TRUE/FALSE) indicating whether temporary simulation files and directories from fastsimcoal should be deleted after genetic data are read into R.  If TRUE, \code{strataG::fscCleanup} is used to remove files output during the coalescent simulation.
#' @param num_cores Number of processors to use for fastsimcoal simulations.
#' @param exec The filename for the fastsimcoal executable.
#' @param loc_parms A one row data frame of parameters providing details on the simulated genetic markers.  Contains the following named elements: \code{marker} (the type of marker to simulate, e.g., "snp"), \code{nloci} (the final number of loci included in the output dataset), \code{seq_length} (the assumed sequence length of simulated loci, if \code{marker} = "snp"), \code{mu} (the mutation rate for the simulated marker).
#' @param found_Ne A parameter describing the effective population size associated with newly colonized populations in the coalescent simulation.
#' @param gmap A map specifying a scheme for aggregating cells from the finer grained forward demographic simulation for more efficient coalescent simulation.
#' @param MAF A numeric value specifying a minor allele frequency filter.  SNP loci with minor alleles below this frequency are excluded from output.
#' @param maxloc The maximum number of loci to attempt to simulate in fastsimcoal.  If too few variable sites are generated, additional markers are simulated up to this value.
#'
#' @details 
#' Provides a wrapper for several helpful functions from the \code{strataG} R package.  Input parameters include a pophist object (from \code{getpophist2.cells}), a landscape object, and parameter values associated with refugial populations and genetic markers. Population size changes are limited to step-changes from founding population size (\code{found_Ne}) to the population size at the end of the forward simulation (from elements of \code{ph$Nvecs}). To improve computational efficiency, coalescent simulations are conducted for aggregated groups of cells, with an aggregation scheme defined by the function \code{make.gmap}. SNP data output by this function are subsampled to retain only one variable position from each simulated DNA segment (locus).  \emph{Note that the fastsimcoal executable must be installed and in a directory in the user's $PATH for coalescent simulations.}
#'
#' @return
#' Returns a data frame of individual SNP genotypes from fastsimcoal.  The first column of this data frame provides an individual ID, the second specifies the population or deme associated with the individual, and columns that follow provide diploid genotypes for simulated individuals (in two-column format).
#'
#' @examples
#' library(holoSimCell)
#' parms <- drawParms(control = system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell"))
#' modchoice <- 1
#' load(file=paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",pollenPulls[[modchoice]]$file))
#' refpops <- pollenPulls[[modchoice]]$refs
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
#' gmap=make.gmap(ph$pophist,
#'                xnum=2, #number of cells to aggregate in x-direction
#'                ynum=2) #number of aggregate in the y-direction
#' 
#' ph2 <- pophist.aggregate(ph,gmap=gmap)
#'
#' loc_parms <- data.frame(marker = "snp",
#'                         nloci = parms$nloci,           
#'                         seq_length = parms$seq_length,
#'                         mu = parms$mu)
#'   
#' preLGMparms <- data.frame(preLGM_t = parms$preLGM_t/parms$G,   
#'                           preLGM_Ne = parms$preLGM_Ne,
#'                          ref_Ne = parms$ref_Ne)
#' 
#' out <- runFSC_step_agg3(ph = ph2,
#'                         l = landscape,
#'                         sample_n = 14,
#'                         preLGMparms = preLGMparms,
#'                         label = "test",
#'                         delete_files = TRUE,
#'                         num_cores = 1,
#'                         exec = "fsc26",
#'                         loc_parms = loc_parms,
#'                         found_Ne = parms$found_Ne,
#'                         gmap = gmap,
#'                         MAF = 0.01,
#'                         maxloc = 50000)
#'
#' @seealso \code{\link[strataG]{fscTutorial}}, \code{\link[strataG]{fscWrite}}, \code{\link[strataG]{fscRun}}, \code{\link[strataG]{fscReadArp}}, \code{\link{getpophist2.cells}}, \code{\link{pophist.aggregate}}, \code{\link{make.gmap}}, \code{\link{holoStats}}, \url{http://cmpg.unibe.ch/software/fastsimcoal26/man/fastsimcoal26.pdf}
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
  pops <- ph$pophist
  simhist <- ph$Nvecs
  tmat <- ph$tmat
  struct <- ph$struct
  coalhist <- ph$coalhist
  
  migmat <- vector("list", 1)
  migmat[[1]] <- tmat
  diag(migmat[[1]]) <- 0
  
  allevents <- data.frame(time=NA, source = NA, sink = NA, prop_m = NA, new_size = NA, new_growth=NA, new_mig = NA)
  
  ############################################################
  
  
  ############################################################
  #Deal with empty populations 
  empty <- rep(FALSE, nrow(gmap))
  empty[l$empty] <- TRUE
  empty_gpops <- sort(unique(c(                                    
    gmap$gpop[!gmap$gpop %in% gmap$gpop[empty == FALSE]],
    which(rowSums(ph$Nvecs) == 0),
    gmap$gpop[!gmap$gpop %in% coalhist$src & !gmap$gpop %in% coalhist$snk])))

  
  if(length(empty_gpops) > 0) {
    simhist <- simhist[-empty_gpops,]
    migmat[[1]] <- migmat[[1]][-empty_gpops,]
    migmat[[1]] <- migmat[[1]][,-empty_gpops]
  } 
  ############################################################
  
  
  ############################################################
  #Build objects with corrected IDs, backward in time arrivals, etc.
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
  h <- max(gmap$gpop)
    
  pop_size <- simhist[,length(simhist[1,])]
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

  sample_size <- rep(0,h)
  growth_rate <- rep(0,h)
  
  if(length(empty_gpops > 0)) {
    sample_size <- sample_size[-empty_gpops]
    growth_rate <- growth_rate[-empty_gpops]
  }
  
  sample_ids <- plyr::mapvalues(l$sampled, gmap$pop, gmap$gpop, warn_missing=FALSE)
  if(sum(sample_ids %in% empty_gpops) > 0) {
    stop("Some locations where genetic samples have been collected were not colonized during the forward simulation.  Try other parameter values and rerun getpophist!")
  }
  sample_pops <- plyr::mapvalues(sample_ids, oldID, newID, warn_missing=FALSE)
  sample_size[sample_pops+1] <- sample_n
 
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
  
  
  ##################################################
  # PARAMETERS FOR SIMULATED LOCI
  ##################################################
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
    tmp_leaving <- c()
    
    if(length(tmp_pops[,1]) > 0) {
      for(fp in 1:length(tmp_pops[,1])) {
        if(!is.na(tmp_pops$snk[fp])) {
          if(tmp_pops$prop[fp] == 1) {
            allevents[ev,] <- c(G-1, tmp_pops$src[fp], tmp_pops$src[fp], 0, as.numeric(found_Ne)/as.numeric(pop_size[tmp_pops$src[fp]+1]), 0, usemigmat)
            ev <- ev+1
            tmp_leaving <- c(tmp_leaving, tmp_pops$src[fp]+1)
          }
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
  
  eventlist <- vector("list", nrow(allevents))
  for(ev in 1:nrow(allevents)) {
    eventlist[[ev]] <- fscEvent(event.time = allevents$time[ev], source = allevents$source[ev], sink = allevents$sink[ev], prop.migrants = allevents$prop_m[ev], new.size = allevents$new_size[ev], new.growth = allevents$new_growth[ev], migr.mat = allevents$new_mig[ev])
  }

  if(length(newID[which(simhist[,1] > 0)]) == 1) {
    eventlist[[length(eventlist)+1]] <- fscEvent(event.time = preLGMparms$preLGM_t[1], source = newID[which(simhist[,1] > 0)], sink = newID[which(simhist[,1] > 0)], prop.migrants = 0, new.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.growth = 0, migr.mat = nomig_migmat)
  } else {
    arb.preLGM.pop <- newID[which(simhist[,1] > 0)][1]
    for(rem in 2:length(which(simhist[,1] > 0))) {
      eventlist[[length(eventlist)+1]] <- fscEvent(event.time = preLGMparms$preLGM_t[1]-1, source = newID[which(simhist[,1] >0)][rem], sink = arb.preLGM.pop, prop.migrants = 1, new.size = 1, new.growth = 0, migr.mat = nomig_migmat)
    }
    eventlist[[length(eventlist)+1]] <- fscEvent(event.time = preLGMparms$preLGM_t[1], source = arb.preLGM.pop, sink = arb.preLGM.pop, prop.migrants = 0, new.size = round(preLGMparms$preLGM_Ne/preLGMparms$ref_Ne,5), new.growth = 0, migr.mat = nomig_migmat)
  }
  events <- do.call(fscSettingsEvents, eventlist)

  varSNPs <- 0
  while(varSNPs < loc_parms$nloci) {
    p <- fscWrite(demes = demes, genetics = genetics, events = events, migration = migration, label = label, use.wd = TRUE)
    p <- fscRun(p, all.sites = FALSE, inf.sites = FALSE, no.arl.output = FALSE, dna.to.snp = TRUE, quiet = FALSE, num.cores = num_cores, exec = exec)
    out <- fscReadArp(p)
    
    if (delete_files==TRUE)
    {
        print(date())
        print(paste("cleaning up fsc files:",label))
        fscCleanup(label)
    }
    
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
  
  popid_inds <- sapply(strsplit(out$id, split = "_"), FUN = function(x) as.numeric(x[1]))
  fscout$deme <- demes$deme.name[popid_inds]

  fscout <- fscout[order(fscout$deme),]
  
  message(paste("Coalescent simulation complete -", (ncol(fscout)-2)/2, "SNPs generated"), appendLF=TRUE)
  fscout
  
  
}
