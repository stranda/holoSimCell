##
## functions that support holostat, primarily
##

mafreq <- function(out)
{
    #locs <- names(out@data)[-1:-2]
    locs <- regmatches(colnames(out)[-c(1:2)], regexpr("^C[[:digit:]]+", colnames(out)[-c(1:2)]))
    locs <- unique(locs)
    sapply(locs,function(l)
    {
        #tbl=table(out@data[[l]])
        tbl <- table(c(out[,grep(l,colnames(out))[1]], out[,grep(l,colnames(out))[2]]))
        res= (sort(tbl)/sum(tbl))[1]
        res
    })
}

loc.mafreq <- function(split_out = NULL, minor = NULL) {
  outm <- matrix(data = NA,
                 nrow = length(grep("^C[[:digit:]]+",names(split_out[[1]])))/2,
                 ncol = length(split_out))
  rownames(outm) <- unique(gsub(".1$|.2$","",names(split_out[[1]])[grep("^C[[:digit:]]+",names(split_out[[1]]))]))
  colnames(outm) <- names(split_out)
  
  for(pop in names(split_out)) {
    
    geno <- split_out[[pop]][,-c(1,2)]
    #seq(1,ncol(geno)/2,by=1)
    geno <- apply(geno,2,as.character)
    geno <- apply(geno,2,as.numeric)
    
    for(mloc in gsub(".1$|.0$","",names(minor))) { #cruise through the minor alleles
      outm[grep(mloc,rownames(outm)),colnames(outm)==pop] <-
        length(which(geno[,grep(mloc,colnames(geno))] == minor[grep(mloc,names(minor))]))/(nrow(geno)*2)
    }
  }
  outm
}


pwise.fst.glob <- function(locMAF,allMAF,locN)  #uses the Nei formulation parameterized in adegenet docs
{
    ahet = mean(2*allMAF*(1-allMAF))
    phet = colMeans(2*locMAF*(1-locMAF))
    d = matrix(NA,nrow=length(phet),ncol=length(phet))
    for (row in 2:length(phet))
        for (col in 1:row)
        {
            d[row,col] <- (ahet-(locN[row]*phet[row] + locN[col]*phet[col])/(locN[row]+locN[col])) / ahet
        }
    diag(d)=0
    d
}

pwise.fst.loc <- function(locMAF,allMAF,locN,pHetTot)  #uses the actual calc in adegenet pairwise.fst function
{

    phet = colMeans(2*locMAF*(1-locMAF))
    d = matrix(NA,nrow=length(phet),ncol=length(phet))
    for (row in 2:length(phet))
        for (col in 1:row)
        {
             ahet=pHetTot[row,col]
            d[row,col] <- (ahet-(locN[row]*phet[row] + locN[col]*phet[col])/(locN[row]+locN[col])) / ahet
        }
    diag(d)=0
    d
}


pwise.het.slowest <- function(out) #pooled het among pairs of pops
{
    ret = matrix(0,nrow=length(unique(out@data$strata)),
                 ncol=length(unique(out@data$strata)))
    for (i in 2:dim(ret)[1])
        for (j in 1:(i-1))
        {
#            print(c(i,j))
            frq=mafreq(out[,,c(i,j)])
            ret[i,j]=mean(2*frq*(1-frq))
        }
    ret+t(ret)
}

pwise.het.slow <- function(out) #pooled het among pairs of pops
{
    strata=unique(out@data$strata)
    
    ret = matrix(0,nrow=length(strata),
                 ncol=length(strata))
    
    outdf <- do.call(cbind,out@data)
    
    for (i in 2:dim(ret)[1])
        for (j in 1:(i-1))
        {
###            print(c(i,j))
            tdf=outdf[outdf[,2] %in% c(strata[i],strata[j]),]
            maf=apply(tdf,2,function(x){tbl=table(x);(tbl/sum(tbl))[1]})
            ret[i,j] <-  mean(2*maf*(1-maf))
        }
    ret+t(ret)
}

#pooled het among pairs of pops
pwise.het <- function(locMAF,locN,cores=1)
{
    strata=colnames(locMAF)
    df=expand.grid(x=strata,y=strata)
    df$x=as.character(df$x)
    df$y=as.character(df$y)
    
    name_map <- unique(c(df[,1],df[,2]))
    df$xn=as.numeric(sapply(df$x,function(n){which(name_map==n)}))
    df$yn=as.numeric(sapply(df$y,function(n){which(name_map==n)}))
    
    alldf=df
    df=df[df$x!=df$y,]
    
    reduc=data.frame(unique(t(apply(df[,1:2],1,sort))))
    names(reduc) <- c("x","y")
    df=merge(reduc,df)
    
    df$x=as.character(df$x)
    df$y=as.character(df$y)
    
    

    df$het <- unlist(mclapply(1:dim(df)[1],mc.cores=cores,function(i)
    {
        tlm = locMAF[,colnames(locMAF) %in% c(df$x[i],df$y[i])]
        ns = locN[colnames(locMAF) %in% c(df$x[i],df$y[i])]
        jaf = (tlm[,1]*ns[1]+tlm[,2]*ns[2]) / sum(ns)
        mean(2*jaf*(1-jaf))
    }))

    df <- merge(df,alldf,all.y=T)
    df=df[order(df$y,df$x),]
    dfmat = matrix(df$het,nrow=sqrt(dim(df)[1]))
    rownames(dfmat) <- df$x[1:dim(dfmat)[1]]
    colnames(dfmat) <- rownames(dfmat)
    dfmat[is.na(dfmat)]=0
    dfmat+t(dfmat)
}

##' process genotype data
##' @param out genotype data

get.gSum = function(out)
{
    nall=sapply(out@data[names(out@data)[grep("Locus",names(out@data))]],
                function(x) {length(levels(x))})
    nvar=sum(nall>1)
    
    list(nall=nall,nvar=nvar)
}

##' get number of snps
##' @param out genotype data

get.nSNP = function(out) {
    get.gSum(out)$nvar
}

fsc.rename = function(out,popDF,dip=F)
{
    strat_id = c()
    mult = as.numeric(dip)+1
    for(popu in 1:length(popDF$id)) {
        if(popu == 1) {
            strat_id = rep(as.character(popDF$id[popu]), popDF$sample.size[popu]*mult)
        } else {
            strat_id = append(strat_id, rep(as.character(popDF$id[popu]), popDF$sample.size[popu]*mult))
        }
    }
    out@data$strata <- strat_id
    out
}


pwise.nei <- function(locMAF,cores=1)
{
    df=expand.grid(x=colnames(locMAF),y=colnames(locMAF))
    df$x=as.character(df$x)
    df$y=as.character(df$y)

    name_map <- unique(c(df[,1],df[,2]))
    df$xn=as.numeric(sapply(df$x,function(n){which(name_map==n)}))
    df$yn=as.numeric(sapply(df$y,function(n){which(name_map==n)}))

    alldf=df
    df=df[df$x!=df$y,]
 
    reduc=data.frame(unique(t(apply(df[,1:2],1,sort))))
    names(reduc) <- c("x","y")
    df=merge(reduc,df)

     df$x=as.character(df$x)
    df$y=as.character(df$y)
  
    L=dim(locMAF)[1]
    df$nei=unlist(mclapply(1:dim(df)[1],mc.cores=cores,function(i)
    {
        tlm = locMAF[,colnames(locMAF) %in% c(df$x[i],df$y[i])]
        Jx = sum(tlm[,1]^2+(1-tlm[,1])^2)/L
        Jy = sum(tlm[,2]^2+(1-tlm[,2])^2)/L
        Jxy= sum((tlm[,1]*tlm[,2])+((1-tlm[,1])*(1-tlm[,2])))/L
        -log(Jxy/((Jx*Jy)^0.5))
    }
    ))
    df <- merge(df,alldf,all.y=T)
    df=df[order(df$y,df$x),]
    dfmat = matrix(df$nei,nrow=sqrt(dim(df)[1]))
    rownames(dfmat) <- df$x[1:dim(dfmat)[1]]
    colnames(dfmat) <- rownames(dfmat)
    d <- as.dist(t(dfmat))
    d
}

#Modified function from Eric Archer's strataG fscTutorial() markdown
#Adding pair.per.loc, simulating diploid data, need to retain 2 columns for each locus
sampleOnePerLocus <- function(mat, MAF = NULL) {
  # Extract the SNP names from the matrix column names
  snp.name <- colnames(mat[, -(1:2)])
  snpcol1 <- snp.name[grep(pattern = ".1", snp.name, fixed = TRUE)]
  # Extract the chromosome name (starts with "C" and is followed by numbers) 
  #   from the SNP names
  chrom.names <- regmatches(snpcol1, regexpr("^C[[:digit:]]+", snpcol1))
  
  if(!is.null(MAF)) {
    minorfreq <- c()
    for(col in 1:length(snpcol1)) {
      tmp_locuscalls <- c(mat[,snpcol1[col]],mat[,paste0(strsplit(snpcol1[col], split = ".1", fixed = TRUE), ".2")])
      tmp_locustab <- table(tmp_locuscalls)
      minorfreq[col] <- min(tmp_locustab/sum(tmp_locustab))
      rm(tmp_locuscalls, tmp_locustab)
    }
    filtered <- snpcol1[minorfreq >= MAF & minorfreq <= 0.5]
    filt.chrom.names <- regmatches(filtered, regexpr("^C[[:digit:]]+", filtered))
    one.per.loc <- tapply(filtered, filt.chrom.names, sample, size = 1)
  } else {
    one.per.loc <- tapply(snpcol1, chrom.names, sample, size = 1)  
  }
  # Choose one SNP per chromosome
  
  for(l in 1:length(one.per.loc)) {
    if(l == 1) {
      pair.per.loc <- c(one.per.loc[l], paste0(strsplit(one.per.loc[l], split = ".1", fixed = TRUE), ".2"))
    } else {
      pair.per.loc <- c(pair.per.loc, one.per.loc[l], paste0(strsplit(one.per.loc[l], split = ".1", fixed = TRUE), ".2"))
    }
  }
  # Return matrix of 
  mat[, c("id", "deme", pair.per.loc)]
}
