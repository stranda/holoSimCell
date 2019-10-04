#' Calculate summary stats for holosim
#'
#' @param out genetic data
#' @param popDF population history
#' @param extent size of landscape
#' @return dataframe with stats
#' @export
#' 
#holostats function
holoStats = function(out, popDF, extent, cores=1) {
  pops.xy = popDF[,c("id","col","row")]
  names(pops.xy)[1] <- "pop"

  allMAF <- mafreq(out)
  totalHe <- 2*allMAF*(1-allMAF)
  popid <- strataNames(out)
  
  SNPs <- sum(allMAF<1)
  names(SNPs) <- "tot_SNPs"
  split_out <- strataSplit(out) #list of strataG objects for each pop
  ###    locMAF = sapply(split_out,function(o){mafreq(o)})
#!#  locMAF <- do.call(cbind,mclapply(split_out,mc.cores=cores,function(o){mafreq(o)}))

  #!# Error here in locMAF -- just because it's the local minor allele doesn't mean it's the global minor allele
  #!# Exploit the names of allMAF to ID the global minor allele
  minor <- sapply(names(allMAF), FUN=function(x){strsplit(x, split = "[.]")[[1]][2]})
  #!# A new function for the locMAF -- this is a big change... does it work with all things that depend on locMAF??
  loc.mafreq <- function(split_out = NULL, minor = NULL) {
    out <- matrix(data = NA, nrow = length(grep("Locus",names(split_out[[1]]@data))), ncol = length(split_out))
    rownames(out) <- names(split_out[[1]]@data)[grep("Locus",names(split_out[[1]]@data))]
    colnames(out) <- names(split_out)

    for(pop in 1:length(split_out)) {
      
      geno <- split_out[[pop]]@data[,-c(1,2)]
      geno <- apply(geno,2,as.character)
      geno <- apply(geno,2,as.numeric)

      for(loc in 1:length(out[,1])) {
        out[loc,pop] <- length(which(geno[,loc] == minor[loc]))/length(geno[,loc])
      }

      rm(geno)
    }

    out
  }

  locMAF <- loc.mafreq(split_out, minor)


#  colnames(locMAF) <- names(split_out)
  locHe <- colMeans(2*locMAF*(1-locMAF))
  varlocHe <- apply(2*locMAF*(1-locMAF),2,var)
  
  locN <- sapply(split_out,function(o){length(o@data$ids)})
#  names(locN) <- popid
  #!# Change here for new locMAF
  #localSNP <- apply(locMAF,2,function(x){sum(x<1)})
  localSNP <- apply(locMAF,2,function(x){sum(x<1 & x>0)})

#  names(localSNP) <- paste0("S.", popid)
  names(localSNP) <- paste0("S.", names(localSNP))
  ##not sure we need this, does not seem to be variable and costs some time
  privateSNP <- colSums(privateAlleles(out))
  #privateSNP = privateSNP[sample.order]
#  names(privateSNP) <- paste0("pS.", popid)
  names(privateSNP) <- paste0("pS.", names(privateSNP))

  total_priv = sum(privateSNP)
  names(total_priv) <- "tot_priv"
  
  pwhet <- pwise.het(locMAF,locN,cores)
  FstMat.loc <- as.dist(pwise.fst.loc(locMAF,allMAF,locN,pwhet))
  neimat <- pwise.nei(locMAF,cores)
  
  pairnei <- as.vector(neimat)
  neinames <- c()
  #This will be different...
  #for(pid in 1:(length(popDF$id)-1)) {
  #  neinames <- c(neinames, paste0("Nei_", popDF$id[pid],".", popDF$id[(pid+1):length(popDF$id)]))
  #}
  neimat_names <- attr(neimat, "Labels")
  for(pid in 1:(length(neimat_names)-1)) {
    neinames <- c(neinames, paste0("Nei_", neimat_names[pid],".", neimat_names[(pid+1):length(neimat_names)]))
  }
  
  #!# Naming issue here - popid.popid is also used for Fst - but merge requires same names
  names(pairnei) = neinames
  
  #pairFst.loc = as.vector(FstMat.loc)
  #Fstnames.loc = c()
  #This will be different...
  #for(pid in 1:(length(popDF$id)-1)) {
  #  Fstnames.loc = c(Fstnames.loc, paste0("Fst_", popDF$id[pid],".", popDF$id[(pid+1):length(popDF$id)]))
  #}
  attr(FstMat.loc,"Labels") <- attr(neimat,"Labels")  
  FstMat_names <- attr(FstMat.loc,"Labels")
  pairFst.loc = as.vector(FstMat.loc)
  Fstnames.loc = c()
  #This will be different...
  for(pid in 1:(length(FstMat_names)-1)) {
    Fstnames.loc = c(Fstnames.loc, paste0("Fst_", FstMat_names[pid],".", FstMat_names[(pid+1):length(FstMat_names)]))
  }	
  #!# Naming issue here - popid.popid is also used for Nei - but merge requires same names
  names(pairFst.loc) = Fstnames.loc
  
  tot_Fst = as.numeric(statFst(out)$result[1])
  
  names(tot_Fst) = "tot_Fst"
  
  eucdist = dist(t(locMAF))
  paireuc = as.vector(eucdist)
  eucnames = c()
  dnames = attr(eucdist,"Labels")
  #for(pid in 1:(length(popDF$id)-1)) {
  #  eucnames = c(eucnames, paste0(popDF$id[pid],".", popDF$id[(pid+1):length(popDF$id)]))
  #}
  for(pid in 1:(length(popDF$id)-1)) {
    eucnames = c(eucnames, paste0(dnames[pid],".", dnames[(pid+1):length(dnames)]))
  }
	
  names(paireuc) <- eucnames
  ########### some spatially focused stats
  ###     get pop coords:
  #pops1=t(matrix(1:prod(extent),ncol=extent[1]))
  #popcrd=data.frame(t(sapply(popDF$grid.cell,function(i) {which(pops1==i,arr.ind=T)})))
  #names(popcrd)=c("y","x")
  #popDF=cbind(popDF,popcrd)
  popDF$y <- popDF$row
  popDF$x <- popDF$col
	
  ###     get pairwise geo distances
  pdist=as.matrix(dist(popDF[,c("x","y")]))
  colnames(pdist) <- popDF$id
  rownames(pdist) <- colnames(pdist)
  diag(pdist) <- NA
  pdist[upper.tri(pdist)] <- NA
  dsts = as.data.frame(as.table(pdist))
  dsts = dsts[complete.cases(dsts),]
  names(dsts) <- c("to","from","d")
  dsts$from <- as.character(dsts$from)
  dsts$to <- as.character(dsts$to)
  dsts <- dsts[order(dsts$to,dsts$from),]
  
  print("this far")
  
  #sum_stats_gi<-summary(genind_out)
  #numAll=sum_stats_gi$pop.n.all
  #numAll=data.frame(id=names(numAll),na=numAll,stringsAsFactors=F)
  
  #        he_by_pop<-as.vector(as.numeric(lapply(lapply(seppop(genind_out),summary),
  #                                               function(x) mean(x$Hexp))))
  he_by_pop <- colMeans(2*locMAF*(1-locMAF))
  
#  hedf <- data.frame(he=he_by_pop, id=unique(out@data$strata),stringsAsFactors=F)
  hedf <- data.frame(he=he_by_pop, id=names(he_by_pop),stringsAsFactors=F)

  pr=prcomp(t(locMAF))
  pcadf <- data.frame(id=rownames(predict(pr)),pc1=predict(pr)[,1],pc2=predict(pr)[,2],pc3=predict(pr)[,3])
  
  popDF <- merge(pcadf,merge(popDF,hedf))
  
  print("this far2")
  polyfit <- function(p,resp="he",ind="y",ord=1)
  {
    fit <- lm(p[,resp]~poly(p[,ind],ord))
    c(coef(fit))
  }
  
  he.lat.stats <- polyfit(popDF,"he","y",ord=2)
  names(he.lat.stats) <- paste0("helat.",c("int","frst","scnd"))
  he.long.stats <- polyfit(popDF,"he","x",ord=2)
  names(he.long.stats) <- paste0("helong.",c("int","frst","scnd"))
  pc1.lat.stats <- polyfit(popDF,"pc1","y",ord=2)
  pc1.long.stats <- polyfit(popDF,"pc1","x",ord=2)
  names(pc1.lat.stats) <- paste0("pc1lat.",c("int","frst","scnd"))
  names(pc1.long.stats) <- paste0("pc1long.",c("int","frst","scnd"))
  pc2.lat.stats <- polyfit(popDF,"pc2","y",ord=2)
  pc2.long.stats <- polyfit(popDF,"pc2","x",ord=2)
  names(pc2.lat.stats) <- paste0("pc2lat.",c("int","frst","scnd"))
  names(pc2.long.stats) <- paste0("pc2long.",c("int","frst","scnd"))
  
  pc3.lat.stats <- polyfit(popDF,"pc3","y",ord=2)
  pc3.long.stats <- polyfit(popDF,"pc3","x",ord=2)
  names(pc3.lat.stats) <- paste0("pc3lat.",c("int","frst","scnd"))
  names(pc3.long.stats) <- paste0("pc3long.",c("int","frst","scnd"))
  
  fsts = data.frame(fst=pairFst.loc)
  fsts = cbind(fsts,data.frame(t(sapply(strsplit(names(pairFst.loc),"\\."),function(nms){c(from=strsplit(nms[1],"_")[[1]][2],to=nms[2])})),stringsAsFactors=F))
  fsts = fsts[order(fsts$to,fsts$from),]
  
  neis = data.frame(nei=pairnei)
  neis = cbind(neis,data.frame(t(sapply(strsplit(names(pairnei),"\\."),function(nms){c(from=strsplit(nms[1],"_")[[1]][2],to=nms[2])})),stringsAsFactors=F))
  neis = neis[order(neis$to,neis$from),]
  
  
  
  
  edist = data.frame(edist=paireuc)
  edist = cbind(edist,data.frame(t(sapply(strsplit(names(paireuc),"\\."),function(nms){c(from=nms[1],to=nms[2])})),stringsAsFactors=F))
  edist = edist[order(edist$to,edist$from),]
  
  dsts <- merge(merge(merge(dsts,fsts),edist),neis)
  
  print("this far3")
  
  IBDfst <- lm(log(fst+1)~log(d),dsts)
  ibdfst.slope <- c(coef(IBDfst)[2])
  ibdfst.int <- c(coef(IBDfst)[1])
  bsfst <- segmentGLM(c(dsts$d),log(c(dsts$fst+1)))
  
  IBDnei <- lm(log(nei+1)~log(d),dsts)
  ibdnei.slope <- c(coef(IBDnei)[2])
  ibdnei.int <- c(coef(IBDnei)[1])
  bsnei <- segmentGLM(c(dsts$d),log(c(dsts$nei+1)))
  
  
  IBDedist <- lm(log(edist+1)~log(d),dsts)
  ibdedist.slope <- c(coef(IBDedist)[2])
  ibdedist.int <- c(coef(IBDedist)[1])
  bsedist <- segmentGLM(c(dsts$d),log(c(dsts$edist+1)))
  
	
  #Turning off spatial stats for now... some issues (JDR 7/26/19)
  if(FALSE) {
  ##Spatial Statistics: written by Ellie Weise
  ##incorporated Dec 19, 2018
  #getting data into a usable format for the stats to calculate
  df <- out@data
  
  #function used to ensure that 0 represents the major allele and 1 is the minor allele in our data set
  #!#Moved to additional_stats_functions.R
  #maj_min_fix <- function(x){
  #  if(mean(x) > 0.5){
  #    x <- -1*(x-1)
  #  }
  #  x
  #}
  
  #making a data file, pops object for calculating statistics from Alvarado-Serrano scripts
  data <- data_reformat(raw_data = df)
  #pops.xy <- pops[samp_pops,c(2,3)]
  npops <- nrow(pops.xy)
  
  #!#avoid  hard-coding variables like this
  #!#need a way to get nsamples from the dataset 
  #nsamples <- as.numeric(table(data$st))/2
  #Some of these stats may not work for unequal sample size
  #nsamples <- 50
  nsamples <- as.numeric(table(data$st)[1])
  
  #Calculating the AFS function in the Alvarado-Serrano script
  #!# Naming issue here - P1-P25 instead of grid cell popid
  NSS_as <- stats.AFS(data = data,sstype = "nss",nsamples = nsamples,pops.xy = pops.xy)
  SSS_as <- stats.AFS(data = data,sstype = "sss",nsamples = nsamples,pops.xy = pops.xy)
  #message("NSS and SSS done")

  #GSSA Calculation from Alvarado-Serrano paper
  #function in separate script
  geodist <- as.matrix(dist(data.frame(popDF$x, popDF$y)))
  gssa <- gssa_raggedness(out = out, dist_mat = geodist)
    gssa <- t(unname(gssa))
  gssa <- as.vector(gssa)
  names(gssa) = paste0("HRi_",pops.xy$pop)
  #message("GSSA done")

  #Calculating spatial PCA from the Alvarado-Serrano script + our plotting technique
  #changed the plot from the Alvarado-Serrano code, but it needs to be wrapped in the function to work
  sPCA_sum_stats <- sPCA.dist(data, pops.xy, nsamples, cpos=2, cneg=2, plot=F)
  #message("SPCA done")

  #Bray-Curtis Distance
  #using vegan package for this as well
  euc_bray <- as.matrix(eucdist)
  bray_curtis <- vegdist(euc_bray)
  #naming stats and turning it into a vector
  pairbc.loc = as.vector(bray_curtis)
  bcnames.loc = c()
  #This will be different...
  for(pid in 1:(length(popDF$id)-1)) {
    bcnames.loc = c(bcnames.loc, paste0("BrayCurt_",popDF$id[pid],".", popDF$id[(pid+1):length(popDF$id)]))
  }
  #!# Naming issue here: popid.popid also used for Nei and Fst
  names(pairbc.loc) <- bcnames.loc
  
  #Moran's I estimate:
  #getting the distance and fst matrices into moran format
  euc_moran <- as.matrix(1/dist(pops.xy[,2:3]))
  #!# Slight change here, adding names to euc_moran for bookkeeping
  attr(euc_moran, "dimnames") <- list(pops.xy$pop, pops.xy$pop)
  moran <- Moran.I(x = localSNP,weight = euc_moran)
  moran <- t(moran)
  moran.names <- paste0("Moran.",colnames(moran))
  moran <- unlist(c(moran))
  names(moran) <- moran.names
  #message("Moran's I done")

  #calculating boundaries for monmonier's algorithm
  pops_xy1 <- pops.xy[,2:3]
#!# Fix here, these column names need to be capitalized for the monmonier fxn!!
  #colnames(pops_xy1) <- c("x","y")
  colnames(pops_xy1) <- c("X","Y")
  monmonier <- monmonier_SSS(NSS=FstMat.loc, xy=pops_xy1,boundLength = 2)
  monmonier <- t(monmonier)
  monmonier.names <- colnames(monmonier)
  monmonier <- unlist(c(monmonier))
  names(monmonier) <- monmonier.names
  #message("Monmonier done")

  #Mantel Test:
  #mantel test with Alvarado-Serrano function (vegan mantel test + summary stats)
#!# Slight change here to avoid warnings, XY argument to mantel.correlog should just be x and y positions
  #mantel_sum <- mantel_SSS(FstMat.loc,pops.xy)
  mantel_sum <- mantel_SSS(FstMat.loc, pops.xy[,-1])
  #message("Mantel done")

  #semi-variogram stat
  #!# Slight change here, pops.xy needs to only be the row/col location
  #semivar <- semivar_SSS(NSS=localSNP,xy=pops.xy)
  semivar <- semivar_SSS(NSS=localSNP,xy=pops_xy1)
  semivar <- t(semivar)
  #message("Semivariogram done")

  ##calculating LD per population
  LD_pop <- NULL
  i <- 1

  #calculate mean LD for each population
  for(i in 1:npops){
    #!# Change here: names aren't 1:npops
    #loci_pop <- subset(data,data$st == i)
    loci_pop <- subset(data, data$st == pops.xy$pop[i])
    snps <- loci_pop[,-c(1,2)]
    pair_LD <- LD.Measures(snps,data = "H")
    mean_LD_pop <- mean(pair_LD$r2)
    ifelse(i==1,
           LD_pop <- as.data.frame(mean_LD_pop,row.names = paste0("LD",pops.xy$pop[i])),
           LD_pop <- rbind(LD_pop,as.data.frame(mean_LD_pop,row.names = paste0("LD",pops.xy$pop[i])))
    )
  }
  LD_pop2 <- t(unname(LD_pop))
  colnames(LD_pop2) <- paste0("LD_",pops.xy$pop)
  #apply the polyfit function to LD
  #!# Changed... risky to point back at the out object.  LD was made going through pops.xy, stick with that
  #LD_pop1 <- data.frame(id = unique(out@data$strata), ld = LD_pop$mean_LD_pop)
  LD_pop1 <- data.frame(id = pops.xy$pop, ld = LD_pop$mean_LD_pop)
  popDF <- merge(popDF,LD_pop1)
  
  ld.lat.stats <- polyfit(popDF,"ld","y",ord=2)
  ld.long.stats <- polyfit(popDF,"ld","x",ord=2)
  names(ld.lat.stats) <- paste0("ldlat.",c("int","frst","scnd"))
  names(ld.long.stats) <- paste0("ldlong.",c("int","frst","scnd"))
  ld.lat.stats <- t(data.frame(ld.lat.stats))
  ld.long.stats <- t(data.frame(ld.long.stats))
  
  ld_stats <- cbind(LD_pop2,ld.lat.stats,ld.long.stats)
  #message("LD Stats done")

  #psi with Fst matrix (allele frequency clines)
  psi<- psiCalc(locMAF, samplen=nsamples)
  #message("Directionality Index done")

  #calculating graph theory statistics
  gt_sum <- graph_theory(data_frame = FstMat.loc, pops = pops.xy)
  #message("Graph done")

  #naming stats to be put into the big stats list
  pop_num <- gt_sum$node_stats[,1]
  stats <- t(gt_sum$node_stats[,-1])
  colnames(stats) <- pop_num
  statnames <- rownames(stats)
  node_stats <- unmatrix(stats)
  
  stats <- gt_sum$edge_stats
  stats1 <- data.frame(edge_between_igraph = stats$edge_between_igraph)
  #rownames(stats1) <- paste0("node",stats$node1,"_",stats$node2)
  #edgenames = c()
  #for(pid in 1:length(popDF$id)) {
  #  edgenames = c(edgenames, paste0(popDF$id[pid],".", popDF$id[-pid]))
  #}
  
  tmpedge <- data.frame(node1 = rep(NA, length(stats$node1)), node2 = rep(NA, length(stats$node2)))	
  for(x in 1:length(popDF$id)) {
	  tmpedge$node1[stats$node1 == x] <- levels(popDF$id)[x]
  	tmpedge$node2[stats$node2 == x] <- levels(popDF$id)[x]
  }
  edgenames <- paste0(tmpedge$node1,".",tmpedge$node2)
  
  rownames(stats1) <- edgenames
  #!# Naming issue here: Using P1-P25 instead of grid cell ID
  edge_stats <- unmatrix(t(stats1))
  }
  ################################################################
  
  stats = c(SNPs, localSNP, privateSNP, total_priv, pairFst.loc, pairnei,tot_Fst,
            ibdfst.slope=ibdfst.slope,ibdfst.int=ibdfst.int,bsfst.break=bsfst[1],bsfst.ll=bsfst[2],
            ibdedist.slope=ibdedist.slope,ibdedist.int=ibdedist.int,bsedist.break=bsedist[1],bsedist.ll=bsedist[2],
            ibdnei.slope=ibdnei.slope,ibdnei.int=ibdnei.int,bsnei.break=bsnei[1],bsnei.ll=bsnei[2],
            
            he.lat.stats, he.long.stats,
            pc1.lat.stats, pc1.long.stats,
            pc2.lat.stats, pc2.long.stats,
            pc3.lat.stats, pc3.long.stats)
            #NSS_as,SSS_as,gssa,sPCA_sum_stats,pairbc.loc,moran,monmonier,semivar[1,],ld_stats[1,],psi,node_stats,edge_stats)
  
  stats1 = matrix(data=stats, nrow = 1)
  colnames(stats1) = names(stats)
  stats = as.data.frame(stats1)
  message("all stats calculated - returning vector of summary stats")

  stats
}

#holoStats <- cmpfun(holoStats.r)  #requires compiler
# not needed for packages, byte compile happens at install
