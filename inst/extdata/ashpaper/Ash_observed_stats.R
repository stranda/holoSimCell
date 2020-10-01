#Script to calculate observed data for Fraxinus pennsylvanica samples...
#Updated 10/1/2020 - JDR

library(holoSimCell)

rownames(popmap) <- popmap[,1]
table(popmap[gsub("fp","",names(imputed)),2]) #printout popnames and samples

#imputed.pruned holds our observed data
#Drop the same populations as we do during the simulation
#Remove samples down to n = 14 per population
imputed.pruned=imputed[,-which(gsub("fp","",names(imputed))%in%popmap[popmap$abbrev=="Michigan","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="UNK","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="MO1","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="ON1","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="VA1","id"])]
imputed.pruned=imputed.pruned[,-which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev=="MB1","id"])]
removes <- c()
popids <- popmap[gsub("fp","",names(imputed.pruned)),2]
table(popids)
for (a in unique(popids))
{
  if (sum(popids==a)>14)
  {
    removes <- c(removes,sample(which(popids==a),1))
  }
}
imputed.pruned <- imputed.pruned[,-1*removes]

poptbl <- table(popmap[gsub("fp","",names(imputed.pruned)),2])
poptbl

samppts <- pts[pts$abbrev %in% names(poptbl),]
samppts

#Need to build the landscape as we would for a forward simulation
#Some of the stats require information from the structure of the landscape
rs <- brick(paste0(system.file("extdata","rasters",package="holoSimCell"),"/","study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif"))

newrs <- newLandscapeDim(rs,0.45)

icelakesland <- def_grid_pred2(pred=1-newrs,
                                     samps=transSampLoc(samppts,
                                                          range.epsg=4326,
                                                          raster.proj=crs(rs)@projargs),
                                     raster.proj=crs(rs)@projargs
                                     )

landscape <- icelakesland

#Create vectors of population names and cell IDs that will be incorporated into the genotypes object
strata_abbrev <- rep(NA, length(imputed.pruned[1,]))
strata_cells <- strata_abbrev

for(x in 1:length(landscape$sampdf$abbrev)) {
  this_pop <- landscape$sampdf$abbrev[x]
  this_cell <- landscape$sampdf$cell[x]
  strata_abbrev[which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev==this_pop,"id"])] <- this_pop
  strata_cells[which(gsub("fp","",names(imputed.pruned))%in%popmap[popmap$abbrev==this_pop,"id"])] <- this_cell
}

#Now work with the data, get it into strataG
allele1 <- seq(1,((2*dim(imputed.pruned)[1])-1),2)
allele2 <- seq(2,(2*dim(imputed.pruned)[1]),2)
loc_name <- rownames(imputed.pruned)


tobe_gtypes <- matrix(data = NA, nrow = sum(poptbl), ncol = 2*dim(imputed.pruned)[1])
for(ind in 1:sum(poptbl)){
  for(loc in 1:length(imputed.pruned[,1])) {
    tmpcall <- strsplit(as.character(imputed.pruned[loc,ind]),"")[[1]]
    tobe_gtypes[ind,allele1[loc]] <- tmpcall[1]
    tobe_gtypes[ind,allele2[loc]] <- tmpcall[2]
    rm(tmpcall)
    }
}


#RECODE TO 0/1 for maj/min allele
tobe_gtypes2 <- tobe_gtypes
for(loc in 1:length(imputed.pruned[,1])) {
  loc_vars <- table(c(tobe_gtypes2[,allele1[loc]], tobe_gtypes2[,allele2[loc]]))
  if(length(loc_vars)==1) {
    tobe_gtypes2[,allele1[loc]] <- 0
    tobe_gtypes2[,allele2[loc]] <- 0
  } else {
    maj <- names(which(loc_vars == max(loc_vars)))
    tmp1 <- rep(0,length(tobe_gtypes2[,allele1[loc]]))
    tmp2 <- rep(0,length(tobe_gtypes2[,allele1[loc]]))
    tmp1[tobe_gtypes2[,allele1[loc]] != maj] <- 1
    tmp2[tobe_gtypes2[,allele2[loc]] != maj] <- 1
    tobe_gtypes2[,allele1[loc]] <- tmp1
    tobe_gtypes2[,allele2[loc]] <- tmp2
    rm(tmp1, tmp2)
  }
  rm(loc_vars)
}

#Minor allele frequency - only want markers with MAF > 0.01
minor_freq <- c()
for(x in 1:nrow(imputed.pruned)) {
	tmpdat <- c(tobe_gtypes2[,allele1[x]], tobe_gtypes2[,allele2[x]])
	minor_freq[x] <- min(table(tmpdat))/sum(table(tmpdat))
	rm(tmpdat)
}

locsamp <- sample(which(minor_freq >= 0.01), 1000, replace = FALSE)
keepcols <- sort(c(allele1[locsamp], allele2[locsamp]))
keepcols <- c(1,2,keepcols+2)

obs_df <- data.frame(id = paste0(strata_abbrev, "_", colnames(imputed.pruned)),
                           deme = strata_abbrev,
                           tobe_gtypes2)
obs_df <- obs_df[order(obs_df$deme),]
obs_df <- obs_df[,keepcols]
dim(obs_df)

for(x in 1:dim(obs_df)[2]) {
  if(x > 2) {
    obs_df[,x] <- as.numeric(as.character(obs_df[,x]))
  } else {
    obs_df[,x] <- as.character(obs_df[,x])
  }
}

columns <- sort(rep(c(1:1000),2))
columns <- formatC(columns, width = 4, format = "d", flag = "0")
columns <- paste0("C", columns, "_SNP_L01")
columns <- paste0(columns, c(".1",".2"))
columns <- c("id", "deme", columns)
colnames(obs_df) <- columns

popDF <- makePopdf(landscape, "cell")

stats <- holoStats(obs_df, popDF, cores = 1)

#Check that there are no NAs in the output
sum(is.na(stats))

#Write output file of summary stats...
write.table(c(stats), "Ash_obs_1Oct20.csv", sep = ",", quote = FALSE, row.names = FALSE)

save.image(file = "Ash_observed_calcs_1Oct20.Rws")

stats[grep("DW",names(stats))]



