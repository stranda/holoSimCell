#script to calculate AFS statistics based on arlequin result files from SPLATCHE
#Diego F. Alvarado-S.
#March 4, 2014
#Modified: Aug 25, 2014

#########################################################################################################################
#function to catch errors
tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}
#########################################################################################################################



#########################################################################################################################
stats.AFS = function(data, sstype, nsamples, pops.xy){
  
  #function to calculate AFS summary statisitics
  # data: table of SNP sequences generated with my python script called arfile_SNPs_extract.py
  # sstype: type of summary statistics wanted (either 'nss' --> non-spatial ss or 'sss' --> spatial ss)
  # nsamples: number of individual seqs per population (assumes same number of seqs in all populations)
  # pops.xy: coordinates of population samples
  #suppressMessages(require(gdata))
  
  #calculating AFS:    
  AFS.fol = matrix(rep(NA, (ceiling(nsamples%/%2)+1)*nrow(pops.xy)), nrow=nrow(pops.xy))      #the ceiling is required to deal with uneven numbers of samples
  rownames(AFS.fol) = as.character(unique(data$st))
  for (I in 1:nrow(AFS.fol)){
    freqs = colSums(data[which(data[,1]==unique(data[,1])[I]),3:ncol(data)]) / nsamples
    freq.counts = table(freqs)
    head.freqs = as.numeric(names(freq.counts))
    head.AFS = seq(from = 0, to = nsamples, by = 1)/nsamples
    AFS.unfol = rep(0,nsamples+1)
    names(AFS.unfol) = head.AFS
    for (J in 1:length(head.freqs)){
      AFS.unfol[as.character(head.freqs[J])] = as.numeric(freq.counts[J])
    }
    colnames(AFS.fol) = head.AFS[1:((nsamples/2)+1)]
    for (L in 1:(ncol(AFS.fol)-1)){
      AFS.fol[I,L] = AFS.unfol[L] + AFS.unfol[length(AFS.unfol)-L+1]
    }
    AFS.fol[I,ncol(AFS.fol)] = AFS.unfol[length(AFS.fol[1,])]
  }
  
  #univariate AFS stats:
  AFS.unistats = matrix(rep(NA,(nrow(AFS.fol)*2)),ncol=2)
  colnames(AFS.unistats) = c('W.Mean', 'W.SD')
  rownames(AFS.unistats) = rownames(AFS.fol)
  for (x in 1:nrow(AFS.unistats)){
    allele.count <- seq(ncol(AFS.fol))-1 
    faketmpdat <- rep(allele.count, AFS.fol[x,])
    AFS.unistats[x,1] = mean(faketmpdat)
    AFS.unistats[x,2] = sd(faketmpdat)
    rm(faketmpdat)
  }
  AFS.uni = unmatrix(t(AFS.unistats))
  
  #pairwise AFS stats:
  comp = combn(pops.xy$pop,2)
  AFS.pairstats = matrix(rep(NA,(2*ncol(comp))),ncol=2)
  colnames(AFS.pairstats) = c('W.Mean.Diff', 'W.SD.Diff')
  rownames(AFS.pairstats) = seq(nrow(AFS.pairstats))
  AFS.pairdiffs = NULL
  for (z in 1:ncol(comp)){
    diff = abs(AFS.fol[comp[1,z],] - AFS.fol[comp[2,z],])
    AFS.pairdiffs = rbind(AFS.pairdiffs, diff)
    allele.count <- seq(ncol(AFS.fol))-1 
    faketmpdat <- rep(allele.count, diff)
    AFS.pairstats[z,1] = mean(faketmpdat)
    AFS.pairstats[z,2] = sd(faketmpdat)
    rownames(AFS.pairstats)[z] = paste(comp[1,z],comp[2,z],sep="_")
    rm(faketmpdat)
  }
  AFS.pair = unmatrix(t(AFS.pairstats))
  
  AFS.stats = c(AFS.uni,AFS.pair)
  if (sstype == 'nss'){
    return(AFS.stats)
  }
  
  
  #Modification of Smouse and Peakall 2009 spatial autocorrelation method
  rownames(AFS.pairdiffs) = rownames(AFS.pairstats)
  Matrix = matrix(rep(NA,(nrow(pops.xy)^2)),ncol=nrow(pops.xy))
  rownames(Matrix) = rownames(AFS.unistats); colnames(Matrix) = rownames(AFS.unistats)
  Matrix[upper.tri(Matrix, diag=T)] = rep(0,(nrow(AFS.pairdiffs)+nrow(pops.xy)))
  Matrix[lower.tri(Matrix)] = rowSums(AFS.pairdiffs)
  Matrix = Matrix + t(Matrix)
  
  Cov.mat = Matrix - Matrix

  for (c in 1:ncol(Matrix)){
    for (r in 1:nrow(Matrix)){
      Cov.mat[r,c] = (-Matrix[r,c] + (sum(Matrix[,c]) + sum(Matrix[r,]))/nrow(Matrix)
                      - sum(Matrix)/(nrow(Matrix)^2))
      
    }
  }
  
  geo.matrix = as.matrix(dist(cbind(as.numeric(pops.xy[,2]),as.numeric(pops.xy[,3])),diag=T,upper=T))

  geo.breaks = unique(as.integer(as.numeric(attributes(table(geo.matrix))$dimnames$geo.matrix)[-1]))
  r.coeffs = rep(0,length(geo.breaks))
  names(r.coeffs) = geo.breaks
  for (z in 1:length(geo.breaks)){
    X.matrix = Matrix - Matrix
    for (c in 1:ncol(geo.matrix)){
      for (r in 1:nrow(geo.matrix)){
        if (as.integer(geo.matrix[r,c])==geo.breaks[z]){
          X.matrix[r,c] = 1
        }}}
    diag(X.matrix) = rowSums(X.matrix)
    X_C.matrix = X.matrix * Cov.mat
    r.coeffs[z] = (sum(X_C.matrix) - sum(diag(X_C.matrix))) / sum(diag(X_C.matrix))
  }
  
  #Spatial correlation stats:
  r.mean = mean(r.coeffs)	
  r.SD = sd(r.coeffs)	 
  r.max = max(r.coeffs)		
  r.lagmax = min(as.numeric(which(r.coeffs==max(r.coeffs))))	
  
  AFS.spatial.stats = c(r.mean, r.SD, r.max, r.lagmax)
  names(AFS.spatial.stats) = c('r.mean', 'r.SD', 'r.max', 'r.lagmax')
  
  if (sstype == 'sss'){
    return(AFS.spatial.stats)
  }
}
#########################################################################################################################





#########################################################################################################################
sPCA.dist = function(data, pops.xy, nsamples, cpos=2, cneg=0, plot=T){
  
  #function to calculate the eculidean distance in sPCA (Jombart et al. 2008) space between individual SNP sequences
  # data: table of SNP sequences generated with my python script called XXX
  # pops.xy: coordinates of population coordinates
  # nsamples: number of individual seqs per population (assumes same number of seqs in all populations)
  # cpos: number of sPCA global components to retain
  # cneg: number of sPCA local components to retain
  # plot: defines whether sPCA plot for global and local (if significant) axes should be plotted
  
  #getting the data into a gneind format 
  data2 = data[,3:ncol(data)]
  
  names = rep(NA,nrow(data2))
  names <- seq(1,nrow(data2))
  row.names(data2) = names
  pops1 = data[,1]
  DATA = genind(data2,pop=pops1,ploidy=1, type="PA")
  
  XY <- data.frame(pop = data[,1],y = NA, x = NA)
  for(p in 1:length(pops.xy[,1])) {
    XY$y[XY$pop == pops.xy$pop[p]] <- pops.xy$row[p]
    XY$x[XY$pop == pops.xy$pop[p]] <- pops.xy$col[p]
  }

  #running definitive SPCA
  mySpca = spca(obj = DATA, xy = XY[,-1], type = 7, a=2, dmin=0.1, scannf = FALSE, plot.nb = FALSE, nfposi = cpos, nfnega = cneg)

  if (as.logical(plot)){
    plot(mySpca$li[,1],mySpca$li[,2])
  }

  Spca.scores = cbind(XY, mySpca$li)
  Spca.scores

  #within distances
  spca.Wit.mean = rep(NA,nrow(pops.xy))
  spca.Wit.dsd = spca.Wit.mean
  for (j in 1:nrow(pops.xy)){
    comp.set = Spca.scores[Spca.scores[,1] == pops.xy$pop[j],4:ncol(Spca.scores)]
    spca.Wit.mean[j] = mean(dist(comp.set))
    names(spca.Wit.mean)[j] = paste('Spca.Dmean',pops.xy$pop[j], sep='_')
    spca.Wit.dsd[j] = sd(dist(comp.set))
    names(spca.Wit.dsd)[j] = paste('Spca.Dsd',pops.xy$pop[j], sep='_')
  }

  #between distances
  comp = combn(pops.xy$pop,2)
  spca.Bet.mean = rep(NA,ncol(comp))
  names(spca.Bet.mean) = paste(comp[1,],comp[2,],sep='_')
  spca.Bet.dsd = spca.Bet.mean
  for (k in 1:ncol(comp)){
    comp.set = rbind(Spca.scores[Spca.scores[,1]==toString(comp[1,k]),4:ncol(Spca.scores)], Spca.scores[Spca.scores[,1]==toString(comp[2,k]),4:ncol(Spca.scores)])
    spca.Bet.dist = as.matrix(dist(comp.set, upper=T))[1:(nrow(comp.set)/2),((nrow(comp.set)/2)+1):nrow(comp.set)]
    spca.Bet.mean[k] = mean(spca.Bet.dist)
    names(spca.Bet.mean)[k] = paste('Spca.Dmean',comp[1,k],comp[2,k],sep='_')
    spca.Bet.dsd[k] = sd(spca.Bet.dist)
    names(spca.Bet.dsd)[k] = paste('SPCA.Dsd',comp[1,k],comp[2,k],sep='_')
  }

  spca.stats = c(spca.Wit.mean, spca.Wit.dsd, spca.Bet.mean, spca.Bet.dsd)
  return(spca.stats)
}
#########################################################################################################################




#########################################################################################################################
moran_SSS = function(NSS,xy){  
  #function to calculate Moran's I coefficient as spatial summary statistic
  #NSS: vector of non-spatial summary stats used as input
  #xy: matrix of coordinates of sampled populations
  
  xy.dist.inv = as.matrix(1/dist(xy))
  moran.I.TC = tryCatch.W.E(Moran.I(NSS,xy.dist.inv))
  if (!is.null(moran.I.TC$value$message)){
    print(moran.I.TC$warning)
    moran = NaN
  } else{    
    moran = moran.I.TC$value$observed	        
  }
  names(moran) = 'Moran.Beta'
  return(moran)
}
#########################################################################################################################





#########################################################################################################################
semivar_SSS = function(NSS, xy, pdf='F'){
  #function to calculate spatial summary statistics based on the semivariogram
  #NSS: vector of non-spatial summary stats used as input
  #xy: matrix of coordinates of sampled populations
  #pdf: whether to generate or not a pdf
  
  do.pdf = as.logical(pdf)
  
  xy.dist = dist(xy)
  bins = seq(0,max(xy.dist),by=max(xy.dist)/round(sqrt(length(xy.dist))))
  
  vario.TC = tryCatch.W.E(variog(coord=xy, data=NSS, breaks=bins))
  vario.fit = tryCatch.W.E(variofit(vario.TC$value,fix.nugget=F))
  if (!is.null(vario.TC$value$message) || is.null(summary(vario.fit$value)$spatial.component[1])){
    print(vario.TC$warning)
    vario.B = NaN
    vario.N = NaN
    vario.S = NaN
    vario.R = NaN
  } else {
    vario.B = vario.TC$value$beta.ols	
    vario.N = summary(vario.fit$value)$nugget.component   
    vario.S = summary(vario.fit$value)$spatial.component[1] 
    vario.R = summary(vario.fit$value)$spatial.component[2]   
    if (do.pdf){
      pdf('Variogram_plot.pdf')
      plot(vario.TC$value)
      lines(vario.fit)
      dev.off()
    }}
  vario.stats = list('Var.Beta'=vario.B, 'Var.Nugget'=vario.N, 'Var.Sill'=vario.S, 'Var.Range'=vario.R)
  return(vario.stats)
}
#########################################################################################################################




#########################################################################################################################
monmonier_SSS = function(NSS, xy, boundLength, pdf='F'){
  #function to calculate spatial summary statistics based on a monmonier function
  #NSS: vector of non-spatial summary stats used as input
  #xy: matrix of coordinates of sampled populations
  #boundLength: vector of bins' upper limits for semivariogram
  #pdf: whether to generate or not a pdf   
  
  do.pdf = as.logical(pdf)
  
  if (sd(xy[,1])!=0 & sd(xy[,2])!=0){
    trimesh.TC = tryCatch.W.E(summary(tri.mesh(xy[,'X'],xy[,'Y']))$na)
    if (!is.null(trimesh.TC$message) | length(grep('Error',trimesh.TC$values))==0){
      n.tries = nrow(xy)-1
    } else {
      n.tries = summary(tri.mesh(xy[,'X'],xy[,'Y']))$na	
    }}
  if (sd(xy[,1]) == 0){
    xy[,1] = xy[,1] + seq(0.001,(0.001*nrow(xy)), by=0.001)
    n.tries = nrow(xy)-1
  }
  if (sd(xy[,2]) == 0){
    xy[,2] = xy[,2] + seq(0.001,(0.001*nrow(xy)), by=0.001)
    n.tries = nrow(xy)-1
  }
  CN = chooseCN(xy, type=1, plot.nb=F)	
  
  monmonier.TC = tryCatch.W.E(optimize.monmonier(xy, NSS, CN, threshold=0, ntry=n.tries, display.graph=F, bd.length=boundLength))	#for more details on the function check its help site   ###NOTE THAT I TAKE THIS OUT: bd.length=10 of the arguments!!!!
  monmonier = monmonier.TC$value
  if (!is.null(monmonier$message)){
    print(monmonier.TC$warning)
    mon.steps = NaN
    mon.Tdist = NaN
    mon.AVGdist = NaN
    mon.SDdist = NaN
    mon.AVGx = NaN
    mon.SDx = NaN
    mon.AVGy = NaN
    mon.SDy = NaN
  } else {
    mon.lengths = c(monmonier$run1$dir1$values,monmonier$run1$dir2$values)	
    mon.xy = rbind(monmonier$run1$dir1$path,monmonier$run1$dir2$path)	
    mon.steps = length(mon.lengths)
    mon.Tdist = sum(mon.lengths)
    mon.AVGdist = mean(mon.lengths)
    mon.SDdist = sd(mon.lengths)
    mon.AVGx = mean(mon.xy[,'x'])
    mon.SDx = sd(mon.xy[,'x'])
    mon.AVGy = mean(mon.xy[,'y'])
    mon.SDy = sd(mon.xy[,'y'])
    if (do.pdf == 'yes'){
      pdf('Monmonier_plot.pdf')
      plot(monmonier)
      dev.off()
    }}
  mon.stats = list('Mon.Nsteps'=mon.steps, 'Mon.Totdist'=mon.Tdist, 'Mon.Avgdist'=mon.AVGdist, 'Mon.Sddist'=mon.SDdist, 'Mon.AvgX'=mon.AVGx, 'Mon.SdX'=mon.SDx, 'Mon.AvgY'=mon.AVGy, 'Mon.SdY'=mon.SDy)
  return(mon.stats)
}
#########################################################################################################################




#########################################################################################################################
mantel_SSS = function(NSS,xy){
  #function to calculate Mantel correlogram values as spatial summary statistics
  #NSS: vector of non-spatial summary stats used as input
  #xy: matrix of coordinates of sampled populations
  
  Mantel.TC = tryCatch.W.E(mantel.correlog(NSS,XY=xy))
  Mantel = Mantel.TC$value
  if (!is.null(Mantel$message)){
    print(Mantel.TC$warning)
    MantCor.x = NaN
    MantCor.sd = NaN
    MantCor.max = NaN
    MantCor.ml = NaN
  } else {
    MantCor.x = mean(Mantel$mantel.res[,'Mantel.cor'], na.rm=T)
    MantCor.sd = sd(Mantel$mantel.res[,'Mantel.cor'], na.rm=T)
    MantCor.max = max(Mantel$mantel.res[,'Mantel.cor'], na.rm=T)
    MantCor.ml = Mantel$mantel.res[as.numeric(which(Mantel$mantel.res[,'Mantel.cor']==max(Mantel$mantel.res[,'Mantel.cor'],na.rm=T))[1]),'class.index']
  }
  mant.stats = list('Man.AvgCor'=MantCor.x, 'Man.SdCor'=MantCor.sd, 'Man.MaxCor'=MantCor.max, 'Man.MaxLag'=MantCor.ml)
  return(mant.stats)
}
#########################################################################################################################
