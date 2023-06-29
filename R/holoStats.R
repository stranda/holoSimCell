#' Calculate population genetic summary statistics
#'
#' Calculates a variety of summary statistics from simulated SNP data.
#'
#' @param out SNP dataset produced by coalescent simulation (using \code{strataG} and \code{fastsimcoal}) - the first column of this data frame provides an individual ID, the second specifies the population or deme associated with the individual, and columns that follow provide diploid genotypes for simulated individuals (in two-column format)
#' @param popDF a three-column data frame produced by \code{makePopdf()} with population IDs (pop) and column (x) / row (y) locations on the simulated landscape for each population
#' @param cores the number of CPUs to use in parallelized portions of the function.
#'
#' @details This function calculates a wide variety of summary statistics from simulated or observed SNP datasets.  Summary statistics can broadly be grouped into statistics: i) related to levels of genetic diversity within populations, ii) related to levels of differentiation between populations, and iii) summarizing spatial patterns in the dataset.
#' \itemize{
#' Statistics measuring levels of genetic diversity within populations:
#' \item{Number of segregating sites (polymorphic SNPs; S) per population}
#' \item{Number of private segregating sites (pS) per popoulation and total across populations}
#' \item{"Frequency down-weighted marker values", calculated as the sum (across loci) of location-specific minor allele frequencies dividied by global minor allele frequency (see Schonswetter & Tribsch 2005)}
#' \item{Mean and standard deviation of minor allele counts per locus}
#' \item{Mean and standard deviations of Euclidean distances in spatial PCA space, considering within population comparisons (see Alvarado-Serrano & Hickerson 2016)}
#' \item{Levels of linkage disequilibrium among loci within populations (using correlation in allele frequencies within individuals, r^2; Hill 1981)}
#' }
#'
#' \itemize{
#' Statistics measuring levels of genetic differentiation between populations:
#' \item{Pairwise Fst (Wright 1949) from local and combined expected heterozygosity (i.e., Hs and Ht)}
#' \item{Pairwise Nei's genetic distance (Nei 1973)}
#' \item{Mean and standard deviation of pairwise differences in minor allele counts per locus}
#' \item{Mean and standard deviations of Euclidean distances in spatial PCA space, considering between population comparisons (see Alvarado-Serrano & Hickerson 2016)}
#' \item{Pairwise Bray-Curtis dissimilarity calculated from Euclidean distance in allele frequencies between populations}
#' \item{Conditional genetic distance (see Dyer et al. 2010) between populations on a population graph calculated from SNP data}
#' }
#'
#' \itemize{
#' Statistics summarizing spatial patterns in the dataset:
#' \item{Summaries of regressions between geographic distance and genetic distance (Fst, Nei's distance, etc.) from Isolation By Distance analyses}
#' \item{Summaries of polynomial regressions between latitude / longitude and measures of diversity (expected heterozygosity, principal component scores on the first 3 axes, LD)}
#' \item{Summaries of the Geographic Spectrum of Shared Alleles (Harpending's raggedness index, GSSA mean, and GSSA variance per population; see Alvarado-Serrano & Hickerson 2018)}
#' \item{Summaries of spatial autocorrelation in the site frequency spectrum (see Smouse & Peakall 1999 and Alvarado-Serrano & Hickerson 2016)}
#' \item{Measures of spatial autocorrelation in genetic data,  including Moran's I (see Moran 1950) and the beta, sill, nugget, and range of the variogram}
#' \item{Summaries of Monmonier's algorithm - mean and standard deviation of path length, x and y positions of path vertices (see Alvarado-Serrano & Hickerson 2016)}
#' \item{Pairwise directionality index (see Peter & Slatkin 2013)}
#' \item{Per population betweenness and closeness centralities (Freeman 1979) from a population graph calculated from SNP genotype data}
#' }
#'
#' @return 
#' \itemize{
#' Returns a one-row dataframe with variable numbers of columns, representing summary statistics calculated from SNP data. Summaries in the output are named as follows:
#' \item{\code{tot_SNPs}} {The total number of polymorphic sites in the dataset, this should equal the \code{nloci} argument to \code{runFSC_step_agg3()}, and is not used in subsequent analyses}
#' \item{\code{S.<POP_ID>}} {The number of polymorphic SNP loci in each population.}
#' \item{\code{pS.<POP_ID>}} {The number of private polymorphic SNP loci (SNPs that are only variable in the focal population) per population.}
#' \item{\code{DW_<POP_ID>}} {Frequency-weighted marker values (see Schonswetter & Tribsch 2005) for each population.}
#' \item{\code{tot_priv}} {The total number of private SNP loci across populations.}
#' \item{\code{Fst_<POP_ID1>.<POP_ID2>}} {Pairwise Fst (=1-Hs/Ht) between populations.}
#' \item{\code{Nei_<POP_ID1>.<POP_ID2>}} {Pairwise Nei's genetic distance between populations.}
#' \item{\code{ibdfst.*}} {Summaries of linear relationships between genetic distance and geographic distance (slope and intercept of IBD relationship)}
#' \item{\code{bsfst.*}} {Summaries of broken-stick relationships between genetic distance and geographic distance (breakpoint, difference in log-likelihood compared to IBD model)}
#' \item{\code{ibdedist.*}} {Similar to \code{ibdfst.*}, but using Euclidean distance rather than Fst to measure genetic differentiation.}
#' \item{\code{bsedist.*}} {Similar to \code{bsfst.*}, but using Euclidean distance rather than Fst to measure genetic differentiation.}
#' \item{\code{ibdnei.*}} {Similar to \code{ibdfst.*}, but using Nei's genetic distance rather than Fst to measure genetic differentiation.}
#' \item{\code{bsnei.*}} {Similar to \code{bsfst.*}, but using Nei's genetic distance rather than Fst to measure genetic differentiation.}
#' \item{\code{helat.*}} {Summaries of a polynomial model (intercept, first, and second coefficients) relating heterozygosity to latitude}
#' \item{\code{helong.*}} {Summaries of a polynomial model (intercept, first, and second coefficients) relating heterozygosity to longitude}
#' \item{\code{pc*lat.*}} {Summaries of a polynomial model (intercept, first, and second coefficients) relating position on principal component (axis 1, 2, or 3) to latitude}
#' \item{\code{pc*long.*}} {Summaries of a polynomial model (intercept, first, and second coefficients) relating position on principal component (axis 1, 2, or 3) to longitude}
#' \item{\code{W.Mean:<POP_ID>}} {Mean minor allele count (from the 1-D SFS) across loci per population.}
#' \item{\code{W.SD:<POP_ID>}} {Standard deviation of minor allele count (from the 1-D SFS) across loci per population.}
#' \item{\code{W.Mean.Diff:<POP_ID1>_<POP_ID2>}} {Mean pairwise difference in allele counts between two populations (from the 2-D SFS).}
#' \item{\code{W.SD.Diff:<POP_ID1>_<POP_ID2>}} {Standard deviation of pairwise differences in allele counts between two populations (from the 2-D SFS).}
#' \item{\code{r.*}} {Summaries of spatial autocorrelation in the site frequency spectrum - mean and standard deviation, maximum correlation coefficient, lag associated with max correlation.}
#' \item{\code{HRi_<POP_ID>}} {Harpending's raggedness index calculated from the Geographic Spectrum of Shared Alleles for each population.}
#' \item{\code{gssa_mean_<POP_ID>}} {Mean distance of allele sharing for each population, calculated from the Geographic Spectrum of Shared Alleles.}
#' \item{\code{gssa_var_<POP_ID>}} {Variance in the distribution of allele sharing distances for each population, calculated from the Geographic Spectrum of Shared Alleles.}
#' \item{\code{Spca.Dmean_<POP_ID>}} {The mean inter-individual distance in PCA space among individuals within a population, calculated from a spatial PCA analysis.}
#' \item{\code{Spca.Dsd_<POP_ID>}} {The standard deviation in inter-individual distance in PCA space among individuals within a population, calculated from a spatial PCA analysis.}
#' \item{\code{Spca.Dmean_<POP_ID1>_<POP_ID2>}} {The mean inter-individual distance in PCA space between individuals sampled from two populations, calculated from a spatial PCA analysis.}
#' \item{\code{Spca.Dsd__<POP_ID1>_<POP_ID2>}} {The standard deviation in inter-individual distance in PCA space between individuals sampled from two populations, calculated from a spatial PCA analysis.}
#' \item{\code{BrayCurt_<POP_ID1>.<POP_ID2>}} {Bray-Curtis dissimilarity calculated from Euclidean distance in allele frequencies between populations.}
#' \item{\code{Moran.Beta}} {Estimate of Moran's I measuring spatial autocorrelation in genetic data.}
#' \item{\code{Mon.*}} {Summaries of results of Monmonier's algorithm applied to the genetic dataset, mean and standard deviation of path length, x and y positions of path vertices}
#' \item{\code{Var.*}} {Summaries of the variogram measuring spatial autocorrelation in the genetic data - beta, sill, nuggget, and range of the variogram.}
#' \item{\code{LD_<POP_ID>}} {Estimate of linkage disequilibrium within a population, using the correlation in allele frequencies within individuals (r^2).}
#' \item{\code{ldlat.*}} {Summaries of a polynomial model (intercept, first, and second coefficients) relating estiamted LD to latitude.}
#' \item{\code{ldlong.*}} {Summaries of a polynomial model (intercept, first, and second coefficients) relating estiamted LD to longitude.}
#' \item{\code{Psi_<POP_ID1>.<POP_ID2>}} {Peter & Slatkin's (2013) directionality index between a pair of populations.}
#' \item{\code{cGD-<POP_ID1>.<POP_ID2>}} {Conditional genetic distance between populations, from a population graph calculated from SNP genotype data.}
#' \item{\code{bwness-<POP_ID>}} {Betweenness centrality for a population, from a population graph calculated from SNP genotype data.}
#' \item{\code{cness-<POP_ID>}} {Closeness centrality for a population, from a population graph calculated from SNP genotype data.}
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
#' popDF <- makePopdf(landscape,"cell")
#'
#' stats <- holoStats(out, popDF, cores = 1)
#'
#' @seealso \code{\link{makePopdf}}, \code{\link{run_FSC_step_agg3}}, \code{\link[strataG]{privateAlleles}}, \code{\link[stats]{prcomp}}, \code{\link[adegenet]{spca}}, \code{\link[LDcorSV]{LD.Measures}}, \code{\link[igraph]{betweenness}}, \code{\link[igraph]{closeness}}, \code{\link[popgraph]{popgraph}}, \code{\link[vegan]{vegdist}}, \url{https://onlinelibrary.wiley.com/doi/full/10.1111/evo.12202}, \url{https://www.nature.com/articles/6885180}, \url{https://www.biorxiv.org/content/10.1101/457556v1}, \url{https://doi.org/10.1017/S0016672300020553}, \url{https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1469-1809.1949.tb02451.x}, \url{https://www.pnas.org/doi/10.1073/pnas.70.12.3321}, \url{https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12489}, \url{https://onlinelibrary.wiley.com/doi/abs/10.2307/25065429}, \url{https://doi.org/10.1016/0378-8733(78)90021-7}, \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-294X.2010.04748.x}
#' @export
#' 
#holostats function
#holoStats = function(out, popDF, extent, cores=1) {
holoStats = function(out, popDF, cores=1) {
  options(mc.cores = cores)
  pops.xy = popDF[,c("pop","x","y")]
  names(pops.xy) <- c("pop", "col", "row")

  allMAF <- mafreq(out)
  totalHe <- 2*allMAF*(1-allMAF)
  popid <- unique(out$deme)
  
  SNPs <- sum(allMAF<1)
  names(SNPs) <- "tot_SNPs"
  split_out <- vector("list", length(unique(out$deme)))
  for(p in 1:length(unique(out$deme))) {
    split_out[[p]] <- out[out$deme == unique(out$deme)[p],]
    names(split_out)[p] <- unique(split_out[[p]]$deme)
  }

  minor <- sapply(names(allMAF), FUN=function(x){strsplit(x, split = "[.]")[[1]][2]})

  locMAF <- loc.mafreq(split_out, minor)
  DW <- colSums(locMAF/allMAF)
  names(DW) <- paste0("DW_",names(DW))

  locHe <- colMeans(2*locMAF*(1-locMAF))
  varlocHe <- apply(2*locMAF*(1-locMAF),2,var)
  
  locN <- sapply(split_out,function(o){nrow(o)})
  localSNP <- apply(locMAF,2,function(x){sum(x<1 & x>0)})

  names(localSNP) <- paste0("S.", names(localSNP))

  privateSNP <- colSums(privateAlleles(df2gtypes(out, ploidy = 2)))
  privateSNP <- privateSNP[match(unique(out$deme), names(privateSNP))]

  names(privateSNP) <- paste0("pS.", names(privateSNP))

  total_priv = sum(privateSNP)
  names(total_priv) <- "tot_priv"
  
  pwhet <- pwise.het(locMAF,locN,cores)
  FstMat.loc <- as.dist(pwise.fst.loc(locMAF,allMAF,locN,pwhet))
  neimat <- pwise.nei(locMAF,cores)
  
  pairnei <- as.vector(neimat)
  neinames <- c()

  neimat_names <- attr(neimat, "Labels")
  for(pid in 1:(length(neimat_names)-1)) {
    neinames <- c(neinames, paste0("Nei_", neimat_names[pid],".", neimat_names[(pid+1):length(neimat_names)]))
  }
  
  names(pairnei) = neinames
  
  attr(FstMat.loc,"Labels") <- attr(neimat,"Labels")  
  FstMat_names <- attr(FstMat.loc,"Labels")
  pairFst.loc = as.vector(FstMat.loc)
  Fstnames.loc = c()

  for(pid in 1:(length(FstMat_names)-1)) {
    Fstnames.loc = c(Fstnames.loc, paste0("Fst_", FstMat_names[pid],".", FstMat_names[(pid+1):length(FstMat_names)]))
  }	

  names(pairFst.loc) = Fstnames.loc
  
  eucdist = dist(t(locMAF))
  paireuc = as.vector(eucdist)
  eucnames = c()
  dnames = attr(eucdist,"Labels")

  for(pid in 1:(length(popDF$id)-1)) {
    eucnames = c(eucnames, paste0(dnames[pid],".", dnames[(pid+1):length(dnames)]))
  }
	
  names(paireuc) <- eucnames

  pdist=as.matrix(dist(popDF[,c("x","y")]))

  colnames(pdist) <- popDF$pop
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
  
  he_by_pop <- colMeans(2*locMAF*(1-locMAF))
  hedf <- data.frame(he=he_by_pop, pop=names(he_by_pop),stringsAsFactors=F)

  pr=prcomp(t(locMAF))
  pcadf <- data.frame(pop=rownames(predict(pr)),pc1=predict(pr)[,1],pc2=predict(pr)[,2],pc3=predict(pr)[,3])
  
  popDF <- merge(pcadf,merge(popDF,hedf))
  
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
  
  IBDfst <- lm((fst/(1-fst))~log(d),dsts) #changed from the previous line to implement Rousset's (1997) version 
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
  
	
  chrom.names <- unique(regmatches(colnames(out[,-c(1,2)]), regexpr("^C[[:digit:]]+", colnames(out[,-c(1,2)]))))
  df <- matrix(data = NA, nrow = 2*nrow(out), ncol = length(chrom.names)+2)
  df <- as.data.frame(df)
  colnames(df) <- c("id","deme",chrom.names)
  df$id <- c(out$id, out$id)
  df$deme <- c(out$deme, out$deme)
  for(loc in chrom.names) {
    tmp.columns <- grep(loc, colnames(out))
    if(length(tmp.columns) == 2) {
      tmp <- c(out[,tmp.columns[1]], out[,tmp.columns[2]])
      df[,loc] <- tmp
      rm(tmp.columns, tmp)
    } else {
      stop(paste("too many columns for chromosome", loc, "subsample to one per locus?"))
    }
  }
  myorder <- c()
  for(pop in unique(out$deme)) {
    tmp <- which(df$deme == pop)
    tmp <- tmp[order(df$id[tmp])]
    myorder <- c(myorder, tmp)
    rm(tmp)
  }
  df <- df[myorder,]
  
  data <- data_reformat(raw_data = df)

  npops <- nrow(pops.xy)
  
  nsamples <- as.numeric(table(data$st)[1])
  
  NSS_as <- stats.AFS(data = data,sstype = "nss",nsamples = nsamples,pops.xy = pops.xy)
  SSS_as <- stats.AFS(data = data,sstype = "sss",nsamples = nsamples,pops.xy = pops.xy)

  geodist <- as.matrix(dist(data.frame(popDF$x, popDF$y)))
  gssa <- gssa_raggedness(data = data, dist_mat = geodist)
  gssa_HRi <- t(unname(gssa$HRi))     
  gssa_HRi <- as.vector(gssa_HRi)
  names(gssa_HRi) = paste0("HRi_",pops.xy$pop)
  gssa_mean <- t(unname(gssa$gssa_mean))
  gssa_mean <- as.vector(gssa_mean)
  names(gssa_mean) = paste0("gssa_mean_", pops.xy$pop)
  gssa_var <- t(unname(gssa$gssa_var))
  gssa_var <- as.vector(gssa_var)
  names(gssa_var) = paste0("gssa_var_", pops.xy$pop)
	  

  sPCA_sum_stats <- sPCA.dist(data, pops.xy, nsamples, cpos=2, cneg=0, plot=F)



  euc_bray <- as.matrix(eucdist)
  bray_curtis <- vegdist(euc_bray)

  pairbc.loc = as.vector(bray_curtis)
  bcnames.loc = c()

  for(pid in 1:(length(popDF$pop)-1)) {
    bcnames.loc = c(bcnames.loc, paste0("BrayCurt_",popDF$pop[pid],".", popDF$pop[(pid+1):length(popDF$pop)]))
  }

  names(pairbc.loc) <- bcnames.loc
  
  moran <- moran_SSS(NSS = localSNP, xy = pops.xy[,c(2:3)])
		     

  pops_xy1 <- pops.xy[,2:3]

  colnames(pops_xy1) <- c("X","Y")
  monmonier <- monmonier_SSS(NSS=FstMat.loc, xy=pops_xy1,boundLength = 2)
  monmonier <- t(monmonier)
  monmonier.names <- colnames(monmonier)
  monmonier <- unlist(c(monmonier))
  names(monmonier) <- monmonier.names

  mantel_sum <- mantel_SSS(FstMat.loc, pops.xy[,-1])


  semivar <- semivar_SSS(NSS=localSNP,xy=pops_xy1)
  semivar <- t(semivar)



  LD_pop <- NULL
  i <- 1


  SNPsamp <- sample(c(1:(ncol(data)-2)), min(500, (ncol(data)-2)), replace = FALSE)
	
  df2 <- matrix(data = NA, nrow = nrow(data)/2, ncol = ncol(data))
  df2 <- as.data.frame(df2)
  colnames(df2) <- colnames(data)
  allele1seq <- seq(1,nrow(data),2)
  rw <- 1
  for(x in allele1seq){
    tmp <- data[c(x,x+1),]
    df2$st[rw] <- tmp$st[1]
    df2$ones[rw] <- tmp$ones[1]
    df2[rw,-c(1,2)] <- colSums(tmp[,-c(1,2)])
    rm(tmp)
    rw <- rw+1
  }  
  rm(rw)

  for(i in 1:npops){
    loci_pop <- subset(df2, df2$st == pops.xy$pop[i])
    snps <- loci_pop[,-c(1,2)]
    snps <- snps[,SNPsamp]
    pair_LD <- LD.Measures(snps,data = "G")
    mean_LD_pop <- mean(pair_LD$r2)
    ifelse(i==1,
           LD_pop <- as.data.frame(mean_LD_pop,row.names = paste0("LD",pops.xy$pop[i])),
           LD_pop <- rbind(LD_pop,as.data.frame(mean_LD_pop,row.names = paste0("LD",pops.xy$pop[i])))
    )
  }
  LD_pop2 <- t(unname(LD_pop))
  colnames(LD_pop2) <- paste0("LD_",pops.xy$pop)

  LD_pop1 <- data.frame(pop = pops.xy$pop, ld = LD_pop$mean_LD_pop)
  popDF <- merge(popDF,LD_pop1)
  
  ld.lat.stats <- polyfit(popDF,"ld","y",ord=2)
  ld.long.stats <- polyfit(popDF,"ld","x",ord=2)
  names(ld.lat.stats) <- paste0("ldlat.",c("int","frst","scnd"))
  names(ld.long.stats) <- paste0("ldlong.",c("int","frst","scnd"))
  ld.lat.stats <- t(data.frame(ld.lat.stats))
  ld.long.stats <- t(data.frame(ld.long.stats))
  
  ld_stats <- cbind(LD_pop2,ld.lat.stats,ld.long.stats)

  psi<- psiCalc(locMAF, samplen=nsamples)

  graphstats <- graph_theory(data = data, stats = c("cGD", "betweenness", "closeness"), plot = FALSE)

 
  stats = c(SNPs, localSNP, privateSNP, DW, total_priv, pairFst.loc, pairnei,
            ibdfst.slope=ibdfst.slope,ibdfst.int=ibdfst.int,bsfst.break=bsfst[1],bsfst.ll=bsfst[2],
            ibdedist.slope=ibdedist.slope,ibdedist.int=ibdedist.int,bsedist.break=bsedist[1],bsedist.ll=bsedist[2],
            ibdnei.slope=ibdnei.slope,ibdnei.int=ibdnei.int,bsnei.break=bsnei[1],bsnei.ll=bsnei[2],
            he.lat.stats, he.long.stats,
            pc1.lat.stats, pc1.long.stats,
            pc2.lat.stats, pc2.long.stats,
            pc3.lat.stats, pc3.long.stats,
            NSS_as,SSS_as,gssa_HRi,gssa_mean,gssa_var,
	    sPCA_sum_stats,pairbc.loc,moran,monmonier,
	    semivar[1,],ld_stats[1,],psi,graphstats)
  
  stats1 = matrix(data=stats, nrow = 1)
  colnames(stats1) = names(stats)
  stats = as.data.frame(stats1)
  message("all stats calculated - returning vector of summary stats")

  stats
}


