#Additional functions to calculate spatial statistics for holostats
#functions all run in Holostats R script

#Fix major (0) and minor (1) allele designations from simulated data
  maj_min_fix <- function(x){
    if(mean(x) > 0.5){
      x <- -1*(x-1)
    }
    x
  }


#Change SNP table into Alvarado-Serrano format
data_reformat <- function(raw_data){
  st <- raw_data$deme
  ones <- rep(1:1,each=length(st))
  act_snps <- raw_data[,-c(1,2)] 
  act_snps <- apply(act_snps,2,as.character)
  act_snps <- apply(act_snps,2,as.numeric)
  act_snps <- apply(act_snps,2,maj_min_fix)
  SNPs_reformat <- as.data.frame(cbind(st,ones,act_snps))
  
  SNPs_reformat[,-1] <- apply(SNPs_reformat[,-1],2,as.character)
  SNPs_reformat[,-1] <- apply(SNPs_reformat[,-1],2,as.numeric)
  SNPs_reformat

}

#Calculate GSSA and raggedness index (Alvarado-Serrano & Hickerson 2018)
gssa_raggedness <- function(data=NULL,dist_mat=NULL){
  ##working with the gtypes object to get in into the correct format
  split <- vector("list", length(unique(data$st)))
  for(p in 1:length(unique(data$st))) {
    split[[p]] <- data[data$st == unique(data$st)[p],]
    names(split)[p] <- unique(split[[p]]$st)
  }
  #First, need to ID the minor alleles for each locus
  maID <- function(x) {
    as.numeric(names(which(table(x) == min(table(x)))))[1]
  }
  ma_vec <- apply(data[,-c(1,2)], 2, maID)
  
  #Next make a table with npops columns and nloci rows to give MAFs per pop and locus
  getMAFreq <- function(x=NULL, ma=NULL) {
    length(which(x == ma))
  }
  
  ma_freq_mat = matrix(data = NA, nrow = length(ma_vec), ncol = length(split))
  for(pop in 1:length(split)) {
    ma_freq_mat[,pop] = mapply(getMAFreq, split[[pop]][,-c(1,2)], ma_vec)
  }
  
  #create an empty list of vectors to fill later
  gssa_vecs = vector("list",length(split))
  
  #filling gssa_vecs with the distances for each allele
  for(pop in 1:length(split)) {    
    #First add the 0's for all minor alleles with freq > 1 in the focal population
    #For each of these there are (n(n-1))/2 0's added to the GSSA vector
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

 
  
  ##loop to make a histogram and calculate raggedness index for each population
  ##also calculates mean, variance, skew and kurtosis
  index <- as.data.frame(rep(0,nrow(dist_mat)))
  colnames(index) <- "HRi"  
  index$gssa_mean <- NA
  index$gssa_var <- NA 
  i <- 1
  for(i in 1:nrow(dist_mat)){
    #make a histogram with the binning scheme
    hist_snps <- hist(gssa_vecs[[i]],breaks = breakpts,plot = FALSE)
    snps_bins <- hist_snps$counts
    #bin the distance matrix to determine a null distribution
    dist_mat_i <- dist_mat[i,]
    hist_dist <- hist(dist_mat_i,breaks = breakpts,plot = FALSE)
    dist_bins <- hist_dist$counts
    
    #bias correction using distance binn scheme
    #takes a linear regression of the snp bin scheme compared to the distance binning scheme
    gen_bin_corr <- abs(lm(as.numeric(snps_bins)~as.numeric(dist_bins))$residuals)
    
    #actually calculate the raggedness index for one population
    non_zero_bins <- as.numeric(which(dist_bins!=0.0))
    snps_bin1 <- gen_bin_corr[non_zero_bins]
    snps_bin1 <- snps_bin1/sum(snps_bin1)    
    comb <- combn(length(snps_bin1),2)
    paired <- combn(snps_bin1,2)
    paired <- paired[,which(colDiffs(comb)==1)]
    diffs2 <- (colDiffs(paired))^2
    ragged <- sum(diffs2)
    index$HRi[i] <- ragged
    
    #Adding calculations of moments off the vector of distances over which alleles are shared
    index$gssa_mean[i] <- mean(gssa_vecs[[i]])
    index$gssa_var[i] <- var(gssa_vecs[[i]])
  }
  index
}

#Function to calculate directionality index - psi - Peter & Slatkin (2013) 
psiCalc = function(data, samplen=20) {
  combos <- combn(colnames(data),2)
  
  combo_names <- apply(combos,2,FUN=function(x){paste0(x[1],".",x[2])})
 
  psi = c()
  pair <- 1
  for(pair in 1:length(combos[1,])) {
    tmp_locMAF <- data[,c(combos[1,pair], combos[2,pair])]
    varloci <- which(rowSums(tmp_locMAF == 1) == 0 & rowSums(tmp_locMAF == 0) == 0)
    if(length(varloci) > 0) {
      tmp_locMAF <- tmp_locMAF[varloci,]
      if(length(varloci) == 1) {
        psi[pair] = tmp_locMAF[1]-tmp_locMAF[2]
      } else {
        psi[pair] = mean(tmp_locMAF[,1])-mean(tmp_locMAF[,2])
      }
    } else {
      psi[pair] = NA
    }
  }
  names(psi) <- paste0("Psi_",combo_names)
  psi
}


#Use popgraph package to calculate conditional genetic distance, betweenness centrality, closeness centrality
graph_theory <- function(data = data, stats = c("cGD", "betweenness", "closeness"), plot = FALSE) {
	indat <- apply(data[,-c(1,2)], 2, as.numeric)
	indat2 <- matrix(data = NA, nrow = nrow(indat)/2, ncol = ncol(indat))
	allele1row <- seq(1,nrow(indat),2)
	for(r in 1:nrow(indat2)){
		indat2[r,] <- indat[allele1row[r],]+indat[(allele1row[r]+1),]
	}
	pops <- as.character(data$st[allele1row])
	lat <- rep(NA, nrow(indat2))
	long <- rep(NA, nrow(indat2))
	colnames(indat2) <- paste0("loc",c(1:ncol(indat2)))
	writeme <- data.frame(Population = pops, Latitude = lat, Longitude = long, indat2)
	if (FALSE) 
	{
	  tmpfilename <- paste0("tmpsnp-",round(runif(1),5),".csv")
	  write.csv(writeme, tmpfilename, quote = FALSE, row.names = FALSE)
	  indat3 <- read_population(tmpfilename, type = "snp", locus.columns = c(4:ncol(writeme)))
	  file.remove(tmpfilename)
	} else {# use text connection
	  vec=  c(paste0(paste(names(writeme),collapse=","),"\n"),sapply(1:nrow(writeme),function(l){paste0(paste(writeme[l,],collapse=", "),"\n")}))
	  indat3 <- read_population(textConnection(vec), type = "snp", locus.columns = c(4:ncol(writeme)))
	}
	
  indat3$Population <- writeme$Population
	indat4 <- to_mv(indat3)
	pops <- indat3$Population
	

	myg <- popgraph(x = indat4, groups = pops)
		
	if(plot == TRUE) {
		plot(myg)
	}

	tmpout <- c()
	if("cGD" %in% stats) {
		cGD <- to_matrix(myg, mode = "shortest path")
		combo_names <- combn(rownames(cGD),2)
		combo_names <- apply(combo_names, 2, FUN = function(x){paste0("cGD-",x[1],".",x[2])})
		cGD_out <- cGD[lower.tri(cGD) == TRUE]
		cGD_out <- as.vector(cGD_out)
		names(cGD_out) <- combo_names
		tmpout <- append(tmpout, cGD_out)
	}

	if("betweenness" %in% stats) {
		bwness <- betweenness(myg)
		names(bwness) <- paste0("bwness-", names(bwness))
		tmpout <- append(tmpout, bwness)
	}	

	if("closeness" %in% stats) {
		cness <- closeness(myg)
		names(cness) <- paste0("cness-", names(cness))
		tmpout <- append(tmpout, cness)
	} 

	tmpout
}



