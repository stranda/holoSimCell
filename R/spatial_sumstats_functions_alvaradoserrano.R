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
  suppressMessages(require(gdata))
  
  #calculating AFS:    
  AFS.fol = matrix(rep(NA, (ceiling(nsamples%/%2)+1)*nrow(pops.xy)), nrow=nrow(pops.xy))      #the ceiling is required to deal with uneven numbers of samples
  rownames(AFS.fol) = seq(nrow(pops.xy))
  for (I in 1:nrow(AFS.fol)){
    #!# Slight change here, pop names aren't 1-25
    #freqs = colSums(data[which(data[,1]==I),3:ncol(data)]) / nsamples
    freqs = colSums(data[which(data[,1]==unique(data[,1])[I]),3:ncol(data)]) / nsamples
    freq.counts = table(freqs)
    head.freqs = as.numeric(names(freq.counts))
    head.AFS = c(0,hist(seq(1:nsamples)/nsamples, breaks=nsamples, plot=F)$breaks)
    AFS.unfol = rep(0,nsamples+1)
    names(AFS.unfol) = head.AFS
    for (J in 1:length(head.freqs)){
      AFS.unfol[as.character(head.freqs[J])] = as.numeric(freq.counts[J])
    }
    colnames(AFS.fol) = head.AFS[1:((nsamples/2)+1)]
    for (L in 1:(ncol(AFS.fol)-1)){
      AFS.fol[I,L] = AFS.unfol[L] + AFS.unfol[length(AFS.unfol)-L]
    }
    AFS.fol[I,ncol(AFS.fol)] = AFS.unfol[length(AFS.unfol)]
    rownames(AFS.fol)[I] = paste('P',I,sep='')
  }
  
  #univariate AFS stats:
  AFS.unistats = matrix(rep(NA,(nrow(AFS.fol)*2)),ncol=2)
  colnames(AFS.unistats) = c('W.Mean', 'W.SD')
  rownames(AFS.unistats) = rownames(AFS.fol)
  for (x in 1:nrow(AFS.unistats)){
    AFS.unistats[x,1] = mean(AFS.fol[x,] * seq(ncol(AFS.fol)))
    AFS.unistats[x,2] = sd(AFS.fol[x,] * seq(ncol(AFS.fol)))
  }
  AFS.uni = unmatrix(t(AFS.unistats))
  
  #pairwise AFS stats:
  comp = combn(nrow(pops.xy),2)
  AFS.pairstats = matrix(rep(NA,(2*ncol(comp))),ncol=2)
  colnames(AFS.pairstats) = c('W.Mean.Diff', 'W.SD.Diff')
  rownames(AFS.pairstats) = seq(nrow(AFS.pairstats))
  AFS.pairdiffs = NULL
  for (z in 1:ncol(comp)){
    diff = abs(AFS.fol[paste('P',comp[1,z],sep=''),] - AFS.fol[paste('P',comp[2,z],sep=''),])
    AFS.pairdiffs = rbind(AFS.pairdiffs, diff)
    AFS.pairstats[z,1] = mean(diff * seq(ncol(AFS.fol)))
    AFS.pairstats[z,2] = sd(diff * seq(ncol(AFS.fol)))
    rownames(AFS.pairstats)[z] = paste(rownames(AFS.fol)[comp[1,z]],rownames(AFS.fol)[comp[2,z]],sep='_')
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
      #using formula (13) from Smouse and Peakall 2009
      Cov.mat[r,c] = (-Matrix[r,c] + (sum(Matrix[,c]) + sum(Matrix[r,]))/length(Matrix)
                      - sum(Matrix)/(length(Matrix)^2))
    }}
  
  #!# CHANGE HERE -- pops.xy first col is pop ID, y and x locations are #2 and #3
  #geo.matrix = as.matrix(dist(cbind(as.numeric(pops.xy[,1]),as.numeric(pops.xy[,2])),diag=T,upper=T))
  geo.matrix = as.matrix(dist(cbind(as.numeric(pops.xy[,2]),as.numeric(pops.xy[,3])),diag=T,upper=T))
  #geo.matrix = as.matrix(dist(pops.xy[,2:3], diag= T, upper = T))   #!# This would be equivalent, somewhat easier to read

  geo.breaks = unique(as.integer(as.numeric(attributes(table(geo.matrix))$dimnames$geo.matrix)[-1]))
  r.coeffs = rep(0,length(geo.breaks))
  names(r.coeffs) = geo.breaks
  for (z in 1:length(geo.breaks)){
    X.matrix = Matrix - Matrix
    for (c in 1:ncol(geo.matrix)){
      for (r in 1:nrow(geo.matrix)){
        #following analogous weighting scheme Fig.2 from Smouse and Peakall 2009
        if (as.integer(geo.matrix[r,c])==geo.breaks[z]){
          X.matrix[r,c] = 1
        }}}
    diag(X.matrix) = rowSums(X.matrix)
    X_C.matrix = X.matrix * Cov.mat
    #using formula (13) from Smouse and Peakall 2009
    r.coeffs[z] = (sum(X_C.matrix) - sum(diag(X_C.matrix))) / sum(diag(X_C.matrix))
  }
  
  #Spatial correlation stats:
  r.mean = mean(r.coeffs)	#standard deviation correlation
  r.SD = sd(r.coeffs)	 
  r.max = max(r.coeffs)		#max abs correlation 
  r.lagmax = min(as.numeric(which(r.coeffs==max(r.coeffs))))	#spatial lag at which max abs correlation occur
  
  AFS.spatial.stats = c(r.mean, r.SD, r.max, r.lagmax)
  names(AFS.spatial.stats) = c('r.mean', 'r.SD', 'r.max', 'r.lagmax')
  
  if (sstype == 'sss'){
    return(AFS.spatial.stats)
  }
}
#########################################################################################################################





#########################################################################################################################
sPCA.dist = function(data, pops.xy, nsamples, cpos=2, cneg=2, plot=T){
  
  #function to calculate the eculidean distance in sPCA (Jombart et al. 2008) space between individual SNP sequences
  # data: table of SNP sequences generated with my python script called XXX
  # pops.xy: coordinates of population coordinates
  # nsamples: number of individual seqs per population (assumes same number of seqs in all populations)
  # cpos: number of sPCA global components to retain
  # cneg: number of sPCA local components to retain
  # plot: defines whether sPCA plot for global and local (if significant) axes should be plotted
  
  suppressMessages(require(adegenet))
  
  #getting the data into a gneind format 
  data2 = data[,3:ncol(data)]
  ###used an alternative way to do names because our data set was incompatible with this naming system
  names = rep(NA,nrow(data2))
  #for (n in 1:nrow(data2)){
  #  if (n/2 > n%/%2){
  #    names[n] = paste(data[n,2],'a',sep='')
  #  } else {
  #    names[n] = paste(data[n,2],'b',sep='')
  #  }}
  names <- seq(1,nrow(data2))
  row.names(data2) = names
  pops1 = data[,1]
  DATA = genind(data2,pop=pops1,ploidy=1, type="PA")
  
  #!# CHANGE HERE: adding pop column to XY object for use later in the b/w and w/i distance calcs
  XY = data.frame(pop=rep(as.character(unique(data[,1])),nsamples),y=rep(pops.xy[,2],nsamples), x= rep(pops.xy[,3],nsamples))
  XY = XY[order(as.numeric(as.character(XY$pop))),]
  #for (i in 1:nrow(pops.xy)){
  #  row = cbind(as.numeric(pops.xy[i,2]),as.numeric(pops.xy[i,3]))
  #  for (j in 1:nsamples){
  #    XY = rbind(XY,row)
  #  }
  #}

  mySpca.prelim = spca(obj = DATA, xy = XY[,-1], type = 7, a=2, dmin=0.1, scannf = FALSE, plot.nb = FALSE)
  mySpca.prelim
  #barplot(mySpca.prelim$eig, main = "Eigenvalues of sPCA")
  
  #testing significance of global and local components
  myGtest <- global.rtest(DATA$tab, mySpca.prelim$lw, nperm = 99)
  if (myGtest$pvalue > 0.05){
    cpos = 2        #NOTE here that even if global structures are not significant, at least the first two components are kept!!!
  }
  #cat(paste('\nGlobal test p-value:', myGtest$pvalue, ' --> ', cpos, 'global components retained\n'))
  #plot(myGtest)
  
  myLtest <- local.rtest(DATA$tab, mySpca.prelim$lw, nperm = 99)
  if (myLtest$pvalue > 0.05){
    cneg = 0
  }
  #cat(paste('Local test p-value:', myLtest$pvalue, ' --> ', cneg, 'local components retained\n'))
  #plot(myLtest)
  
  #running definitive SPCA
  mySpca = spca(obj = DATA, xy = XY[,-1], type = 7, a=2, dmin=0.1, scannf = FALSE, plot.nb = FALSE, nfposi = cpos, nfnega = cneg)
  #cat(paste('Components retained:', mySpca$nfposi, 'pos. and', mySpca$nfnega, 'neg.\n'))
  if (as.logical(plot)){
    plot(mySpca$li[,1],mySpca$li[,2])
  }
    #source(paste(scriptPath,'bounding_ellipse.R',sep='/'))
    #samples = row.names(mySpca$li)
    #color.range = colorRampPalette(rainbow(7))(nrow(pops.xy))
    #colors = rep(NA,length(samples))
    #for (s in 1:length(samples)){
      #samples[s] = strsplit(samples[s], '_')[[1]][1]
      #colors[s] = color.range[as.numeric(strsplit(samples[s], '_')[[1]][1])]
    #}
    #par(mfrow=c(2,1))
    #plot(mySpca$li[,1:2],pch=samples, col=colors, main='Global components scatterplot')
    #for (P in 1:nrow(pops.xy)){
      #mvee(mySpca$li[which(samples==as.character(P)),1:2], color=color.range[P],do.points=F,do.add=T,marker=P)
    #}
    #if (mySpca$nfnega == 2){
      #plot(mySpca$li[,3:4],pch=samples, col=colors, main='Local components scatterplot')
      #for (P in 1:nrow(pops.xy)){
        #mvee(mySpca$li[which(samples==as.character(P)),3:4], color=color.range[P],do.points=F,do.add=T,marker=P)
      #}
    #}
  #}
  
  #obtaining median distance between and within pops
  #!# CHANGE HERE: Need to bind XY, not pops.xy
  #Spca.scores = cbind(pops.xy,mySpca$li)
  Spca.scores = cbind(XY, mySpca$li)
  Spca.scores
  #within distances
  spca.Wit.mean = rep(NA,nrow(pops.xy))
  spca.Wit.dsd = spca.Wit.mean
  for (j in 1:nrow(pops.xy)){
    #!# Names aren't 1-npops, columns 1 through 3 are not PCA scores
    #comp.set = Spca.scores[Spca.scores[,1]==toString(j),2:ncol(Spca.scores)]
    comp.set = Spca.scores[Spca.scores[,1] == pops.xy$pop[j],4:ncol(Spca.scores)]
    spca.Wit.mean[j] = mean(dist(comp.set))
    #!# Slight change to naming stats, instead of (j-1) append pop ID
    names(spca.Wit.mean)[j] = paste('Spca.Dmean',pops.xy$pop[j], sep='_')
    spca.Wit.dsd[j] = sd(dist(comp.set))
    names(spca.Wit.dsd)[j] = paste('Spca.Dsd',pops.xy$pop[j], sep='_')
  }
  #between distances
  #!# Slight chane, combinations of the pop ID
  #comp = combn(nrow(pops.xy),2)
  comp = combn(pops.xy$pop,2)
  spca.Bet.mean = rep(NA,ncol(comp))
  names(spca.Bet.mean) = paste(comp[1,],comp[2,],sep='_')
  spca.Bet.dsd = spca.Bet.mean
  for (k in 1:ncol(comp)){
    #!# As above, slight change here to 4:ncol(Spca.scores) instead of 2:ncol(Spca.scores)
    comp.set = rbind(Spca.scores[Spca.scores[,1]==toString(comp[1,k]),4:ncol(Spca.scores)], Spca.scores[Spca.scores[,1]==toString(comp[2,k]),4:ncol(Spca.scores)])
    spca.Bet.dist = as.matrix(dist(comp.set, upper=T))[1:(nrow(comp.set)/2),((nrow(comp.set)/2)+1):nrow(comp.set)]
    spca.Bet.mean[k] = mean(spca.Bet.dist)
    #!# Slight change in naming here too... same as abovce
    names(spca.Bet.mean)[k] = paste('Spca.Dmean',comp[1,k],comp[2,k],sep='_')
    spca.Bet.dsd[k] = sd(spca.Bet.dist)
    names(spca.Bet.dsd)[k] = paste('SPCA.Dsd',comp[1,k],comp[2,k],sep='_')
  }
  spca.stats = c(spca.Wit.mean, spca.Wit.dsd, spca.Bet.mean, spca.Bet.dsd, myGtest$pvalue, myLtest$pvalue)
  names(spca.stats)[(length(spca.stats)-1):length(spca.stats)] = c('Gtest.pv','Ltest.pv')
  return(spca.stats)
}
#########################################################################################################################




#########################################################################################################################
moran_SSS = function(NSS,xy){  
  #function to calculate Moran's I coefficient as spatial summary statistic
  #NSS: vector of non-spatial summary stats used as input
  #xy: matrix of coordinates of sampled populations
  
  suppressMessages(require(ape)); suppressMessages(require(gdata))
  
  xy.dist.inv = as.matrix(1/dist(xy))
  moran.I.TC = tryCatch.W.E(Moran.I(NSS,xy.dist.inv))
  if (!is.null(moran.I.TC$value$message)){
    print(moran.I.TC$warning)
    moran = NaN
  } else{    
    moran = moran.I.TC$value$observed	        #for details on the function check its help site
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
  
  suppressMessages(require(geoR))
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
    vario.B = vario.TC$value$beta.ols	#for details on the function check its help site
    vario.N = summary(vario.fit$value)$nugget.component   #tausq = nugget
    vario.S = summary(vario.fit$value)$spatial.component[1] #sigmasq = sill
    vario.R = summary(vario.fit$value)$spatial.component[2]   #phi = range
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
  
  suppressMessages(require(adegenet)); suppressMessages(require(tripack))
  do.pdf = as.logical(pdf)
  
  #this step introduces minor displacements to sample coordinates which are required when populations are lined-up in 1-dimension because there is no triangulation possible in that case. The minor distortion do not significantly change the neighbor calculation
  if (sd(xy[,1])!=0 & sd(xy[,2])!=0){
    trimesh.TC = tryCatch.W.E(summary(tri.mesh(xy[,'X'],xy[,'Y']))$na)
    if (!is.null(trimesh.TC$message) | length(grep('Error',trimesh.TC$values))==0){
      n.tries = nrow(xy)-1
    } else {
      n.tries = summary(tri.mesh(xy[,'X'],xy[,'Y']))$na	#provides the number of total possible arcs in the Delauny triangulation to inform the monomier function of all possible start points to try)
    }}
  if (sd(xy[,1]) == 0){
    xy[,1] = xy[,1] + seq(0.001,(0.001*nrow(xy)), by=0.001)
    n.tries = nrow(xy)-1
  }
  if (sd(xy[,2]) == 0){
    xy[,2] = xy[,2] + seq(0.001,(0.001*nrow(xy)), by=0.001)
    n.tries = nrow(xy)-1
  }
  CN = chooseCN(xy, type=1, plot.nb=F)	#This function builds a connection network using in this case a Delauny triangulation
  
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
    mon.lengths = c(monmonier$run1$dir1$values,monmonier$run1$dir2$values)	#stores the lengths of the vectorsof the best break path		
    mon.xy = rbind(monmonier$run1$dir1$path,monmonier$run1$dir2$path)	#stores the lengths of the vectorsof the best break path
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
  
  suppressMessages(require(vegan))
  
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
    MantCor.ml = as.numeric(which(Mantel$mantel.res[,'Mantel.cor']==max(Mantel$mantel.res[,'Mantel.cor'],na.rm=T))[1])
  }
  mant.stats = list('Man.AvgCor'=MantCor.x, 'Man.SdCor'=MantCor.sd, 'Man.MaxCor'=MantCor.max, 'Man.MaxLag'=MantCor.ml)
  return(mant.stats)
}
#########################################################################################################################