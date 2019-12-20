#Additional functions to calculate spatial statistics for holostats
#functions all run in Holostats R script

#Moved here, wasn't working inside holoStats fxn
  maj_min_fix <- function(x){
    if(mean(x) > 0.5){
      x <- -1*(x-1)
    }
    x
  }


#'Rename strata in strataG object
#'
#' renames the strata to match the forward time simulation (I think)
#'
#' @export
nameStrat = function(fscout, pops, sample_pops, sample_n) {
  strat_id = c()
  for(popu in 1:length(sample_pops)) {
    if(popu == 1) {
      strat_id = rep(pops$pophist$pop[sample_pops[popu]], 2*sample_n[popu])
    } else {
      strat_id = append(strat_id, rep(pops$pophist$pop[sample_pops[popu]], 2*sample_n[popu]))
    }
  }
  strat_id <- as.numeric(strat_id)
  str.length <- nchar(strat_id)
  maxlen = max(str.length)
  for(x in 1:(maxlen-1)){
    lead0 <- maxlen-x
    parta <- paste0(rep(0,lead0),collapse ="")
    strat_id[str.length == x] = paste0(parta,strat_id[str.length == x])
  }
  fscout@data$strata <- paste0("pop-",strat_id)
  fscout
}

#function to change SNP table into Alvarado-Serrano format
data_reformat <- function(raw_data){
  #make strata into single number numeric
  #!# This is new, sorts raw_data by strata, as in split_out
  #!# This should help with consistency of the ordering of populations
  raw_data <- raw_data[order(raw_data$strata),]
  #!#
  #CHANGE HERE -- REDOING HOW STRATA ARE NAMED -- JDR May 2019
  #st <- unlist(strsplit(raw_data$strata," "))
  #st <- as.numeric(st[seq(2,length(st),2)])
  st <- raw_data$strata
  #!#

  #make a repeat list of ones
  ones <- rep(1:1,each=length(st))
  #turn SNP data into numeric values
##    act_snps <- raw_data[,!c(1,2)] #AES: I dont think this construct works
    act_snps <- raw_data[,-c(1,2)] #AES: this uses the correct construct (I think); #JDR: Changed to what I normally use -c(1,2)
    act_snps <- apply(act_snps,2,as.character)
  act_snps <- apply(act_snps,2,as.numeric)
  act_snps <- apply(act_snps,2,maj_min_fix)
  #cbind all components into one reformatted matrix
  SNPs_reformat <- as.data.frame(cbind(st,ones,act_snps))
  
  #!# Also added this down here, these were still factors after cbinding
  SNPs_reformat[,-1] <- apply(SNPs_reformat[,-1],2,as.character)
  SNPs_reformat[,-1] <- apply(SNPs_reformat[,-1],2,as.numeric)
  SNPs_reformat
  #!#

}

#GSSA calculation and raggedness index calculation
old_gssa_raggedness <- function(out=NULL,dist_mat=NULL){
  #require(strataG)
  ##working with the gtypes object to get in into the correct format
  full <- out
  split <- strataSplit(out)
  
  #First, need to ID the minor alleles for each locus
  maID <- function(x) {
    as.numeric(names(which(table(x) == min(table(x)))))
  }
  
  ma_vec <- apply(full@data[,-c(1,2)], 2, maID)
  
  #Next make a table with npops columns and nloci rows that will give the frequencies of minor alleles in each population at each locus
  getMAFreq <- function(x=NULL, ma=NULL) {
    length(which(x == ma))
  }
  
  ma_freq_mat = matrix(data = NA, nrow = length(ma_vec), ncol = length(split))
  for(pop in 1:length(split)) {
    ma_freq_mat[,pop] = mapply(getMAFreq, split[[pop]]@data[,-c(1,2)], ma_vec)
  }
  
  #create an empty list of vectors to fill later
  gssa_vecs = vector("list",length(split))
  
  #filling fssa_vecs with the distances for each allele
  for(pop in 1:length(split)) {    
    #First add the 0's for all minor alleles with freq > 1 in the focal population
    #For each of these there are (n(n-1))/2 0's added to the GSSA vector
    
    tmp = rep(0, sum(ma_freq_mat[ma_freq_mat[,pop] > 1,pop]))   #Logical to exclude singletons!
    
    gssa_vecs[[pop]] = append(gssa_vecs[[pop]],tmp)
    
    #Then loop through each of the OTHER pops and add to the GSSA vector
    other = c(1:length(gssa_vecs))
    other = other[-pop]
    for(op in other) {
      times = sum(ma_freq_mat[,pop]*ma_freq_mat[,op])
      gssa_vecs[[pop]] = append(gssa_vecs[[pop]], rep(dist_mat[pop,op],times))
      rm(times)	
    }
  }
  
  #making the historgram for each population and calculate the raggedness index
  index <- as.data.frame(rep(0,nrow(dist_mat)))
  colnames(index) <- "stat"
  
  #calculate the number of bins for the population
  bins <- nclass.Sturges(dist_mat[i,])
  
  #loop to make a histogram and calculate raggedness index for each population
  for(i in 1:nrow(dist_mat)){
    #make a histogram with the binning scheme
    hist_gssa <- hist(gssa_vecs[[i]],breaks = bins)
    #calculate the raggedness index 
    sfs <- hist_gssa$counts
    comb <- combn(length(sfs),2)
    paired <- combn(sfs,2)
    paired <- paired[,which(colDiffs(comb)==1)]
    diffs2 <- (colDiffs(paired))^2
    ragged <- sum(diffs2)/(bins-1)
    index$stat[i] <- ragged
  }
  index
}
gssa_raggedness <- function(out=NULL,dist_mat=NULL){
  #require(strataG)
  #require(matrixStats)
  ##working with the gtypes object to get in into the correct format
  split <- strataSplit(out)
  
  #First, need to ID the minor alleles for each locus
  maID <- function(x) {
    as.numeric(names(which(table(x) == min(table(x)))))[1]
  }
  
  ma_vec <- apply(out@data[,-c(1,2)], 2, maID)
  
  #Next make a table with npops columns and nloci rows that will give the frequencies of minor alleles in each population at each locus
  getMAFreq <- function(x=NULL, ma=NULL) {
    length(which(x == ma))
  }
  
  ma_freq_mat = matrix(data = NA, nrow = length(ma_vec), ncol = length(split))
  for(pop in 1:length(split)) {
    ma_freq_mat[,pop] = mapply(getMAFreq, split[[pop]]@data[,-c(1,2)], ma_vec)
  }
  
  #create an empty list of vectors to fill later
  gssa_vecs = vector("list",length(split))
  
  #filling fssa_vecs with the distances for each allele
  for(pop in 1:length(split)) {    
    #First add the 0's for all minor alleles with freq > 1 in the focal population
    #For each of these there are (n(n-1))/2 0's added to the GSSA vector
################################
# This doesn't do what it says above... and I'm not sure if it should... IT SHOULD!  Pairwise comparisons!
# Should this be allele copy number of 0's?
# In Fig 1 example it's not n choose 2, it's just n (not true, it is pairwise comparisons)
# Maybe this is correct as is, but it would be good to verify this (it was wrong, fixed below now)
################################
    #tmp = rep(0, sum(ma_freq_mat[ma_freq_mat[,pop] > 1,pop]))   #Logical to exclude singletons!
    tmp = rep(0, sum(sapply(ma_freq_mat[ma_freq_mat[,pop] > 1,pop],choose,k=2)))

    gssa_vecs[[pop]] = append(gssa_vecs[[pop]],tmp)
    
    #Then loop through each of the OTHER pops and add to the GSSA vector
    other = c(1:length(gssa_vecs))
    other = other[-pop]
    for(op in other) {
      times = sum(ma_freq_mat[,pop]*ma_freq_mat[,op])
      gssa_vecs[[pop]] = append(gssa_vecs[[pop]], rep(dist_mat[pop,op],times))
      rm(times)	
    }
  }
  
  
  
  #calculate the number of bins for the population
  breakpts = seq(range(dist_mat)[1],range(dist_mat)[2],length=nclass.Sturges(dist_mat))

  #mean, variance, skew and kurtosis(measures spread of a distribution) per population
  
  
  ##loop to make a histogram and calculate raggedness index for each population
  ##also calculates mean, variance, skew and kurtosis(measures spread of a distribution) for each population histogram
  ##output is a grid of all 5 statistics
  #making a blank data frame to fill in the 
  index <- as.data.frame(rep(0,nrow(dist_mat)))
  colnames(index) <- "HRi"
  i <- 1
  for(i in 1:nrow(dist_mat)){
    #make a histogram with the binning scheme
    #!# Change to breakpts from bins, bins is a suggestion, breakpts uses breaks
    #hist_snps <- hist(gssa_vecs[[i]],breaks = bins,plot = FALSE)
    hist_snps <- hist(gssa_vecs[[i]],breaks = breakpts,plot = FALSE)
    snps_bins <- hist_snps$counts
    #bin the distance matrix to determine a null distribution
    dist_mat_i <- dist_mat[i,]
    #!# Change to breakpts from bins, bins is a suggestion, breakpts uses breaks
    #hist_dist <- hist(dist_mat_i,breaks = bins,plot = FALSE)
    hist_dist <- hist(dist_mat_i,breaks = breakpts,plot = FALSE)
    dist_bins <- hist_dist$counts
    
    #bias correction using distance binn scheme
    #takes a linear regression of the snp bin scheme compared to the distance binning scheme
    gen_bin_corr <- abs(lm(as.numeric(snps_bins)~as.numeric(dist_bins))$residuals)
    
    #RESUME REVIEW DOWN HERE!!
    #actually calculate the raggedness index for one population
    non_zero_bins <- as.numeric(which(dist_bins!=0.0))
    snps_bin1 <- gen_bin_corr[non_zero_bins]
    comb <- combn(length(snps_bin1),2)
    paired <- combn(snps_bin1,2)
    paired <- paired[,which(colDiffs(comb)==1)]
    diffs2 <- (colDiffs(paired))^2
    ragged <- sum(diffs2)/(length(non_zero_bins)-1)
    index$HRi[i] <- ragged
  }
  index
}

#psi function from John - May 13, 2019 *JDR*
#!# Changing this to "frequency" - not proportion
#!# Updates - back to frequency in the proportion sense of the word.  I think I misread the paper! - JDR 8/20/19
psiCalc = function(data, samplen=20) {
  combos <- combn(colnames(data),2)
  
  combo_names <- apply(combos,2,FUN=function(x){paste0(x[1],".",x[2])})
 
  psi = c()
  pair <- 1
  for(pair in 1:length(combos[1,])) {
    tmp_locMAF <- data[,c(combos[1,pair], combos[2,pair])]
  #!# This won't work anymore with the changes to locMAF
  #varloci <- which(rowSums(tmp_locMAF == 1) == 0)
    varloci <- which(rowSums(tmp_locMAF == 1) == 0 & rowSums(tmp_locMAF == 0) == 0)
    if(length(varloci) > 0) {
      tmp_locMAF <- tmp_locMAF[varloci,]
      if(length(varloci) == 1) {
        #!# I think I had these wrong... check eqn 1a from Peter & Slatkin 2013
        #psi[pair] = (1/samplen)*((samplen*mean(tmp_locMAF[1]))-(samplen*mean(tmp_locMAF[2])))
        psi[pair] = (1/samplen)*(mean(tmp_locMAF[1])-mean(tmp_locMAF[2]))
      } else {
        #psi[pair] = (1/samplen)*((samplen*mean(tmp_locMAF[,1]))-(samplen*mean(tmp_locMAF[,2])))
        psi[pair] = (1/samplen)*(mean(tmp_locMAF[,1])-mean(tmp_locMAF[,2]))
      }
    } else {
      psi[pair] = NA
    }
  }
  #!# Unnecessary, defined above
  #names(psi) <- apply(combn(colnames(locMAF),2), 2, FUN=function(x) {paste0(x[1],".",x[2])})
  names(psi) <- paste0("Psi_",combo_names)
  psi
}


#OLD AND INCORRECT psi function from John
psiCalc_OLD = function(data, samplen=20) {
  combos <- combn(colnames(locMAF),2)
  
  combo_names <- apply(combos,2,FUN=function(x){paste0(x[1],".",x[2])})
  
  #for(pair in 1:length(combos[1,])) {
  #	tmp_locMAF <- locMAF[,c(combos[1,pair], combos[2,pair])]
  #}
  
  combos <- combn(c(1:length(locMAF[1,])),2)
  psi = c()
  pair <- 1
  for(pair in 1:length(combos[1,])) {
    tmp_locMAF <- locMAF[,c(combos[1,pair], combos[2,pair])]
    #this line didn't do anything results wise so I commented it out fn
    #!# Uncommented this line, it keeps loci that are not variable in either pop out of the mean calculation
    tmp_locMAF <- tmp_locMAF[!rowSums(tmp_locMAF) == 2,]
    psi[pair] = (1/samplen)*(mean(tmp_locMAF[,1])-mean(tmp_locMAF[,2]))
  }
  #!# Unnecessary, defined above
  #names(psi) <- apply(combn(colnames(locMAF),2), 2, FUN=function(x) {paste0(x[1],".",x[2])})
  names(psi) <- paste0("Psi_",combo_names)
  psi
}

#calculating graph theory summary
graph_theory <- function(data_frame,pops){
  #require(igraph)
  #require(sna)
  #making sure there's no negative in the data frame - if there are, turn to zero
  #!# Negatives may be in Fst matrix, but shouldn't be in the df_adj version... would have to have Fst > 1 for that to happen
  #!# JDR - 8/20/19
  data_frame[data_frame < 0] <- 0   #!# There are sometimes negative Fst values in the matrix!
  df_adj <- 1-data_frame
  #!# Turning this off due to above
  #for(i in 1:length(df_adj)){
  #  if(df_adj[i] <= 0){
  #    df_adj[i] = 0
  #  }
  #}
  
  #turn the adjusted data frame into a graph type object
  graph <- graph.adjacency(as.matrix(df_adj),weighted = T)
  
  #making a matrix of by population summary statistics
  gt_sum_pop <- pops[1]
  colnames(gt_sum_pop) <- "pop"
  #not variable in test matrix
  #calculating betweenness of graph nodes
  gt_sum_pop$between_igraph <- betweenness.estimate(graph,cutoff = 0)
  #variable in test matrix
  #calculating betweenness and closeness of data frame
  gt_sum_pop$close_sna <- closeness(as.matrix(graph))
  gt_sum_pop$between_sna <- betweenness(as.matrix(graph))
  #the degree of centrality of each population
  gt_sum_pop$degree_sna <- degree(as.matrix(graph))
  #the eigenvector centrality of each population
  gt_sum_pop$evcent_sna <- evcent(as.matrix(graph))$vector
  #calculating closeness of graph nodes
  gt_sum_pop$close_igraph <- closeness.estimate(graph,cutoff = 0)
  
  #adding in edge statistics
  #calculating edge betweenness for the graph edges (pairwise)
  pairs <- as.data.frame(as_edgelist(graph))
  colnames(pairs) <- c("node1","node2")
  pairs$edge_between_igraph <- edge.betweenness(graph)
  
  #make a summary list object, with one component being the by population node statistics
  #the other component is the pairwise edge statistics
  gf_sum <- list(node_stats = gt_sum_pop,edge_stats = pairs)
}





