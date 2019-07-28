#simrep.R
#Runs the holoSimCell model for Ash using an ENM to define habitat suitability

#Set the "node" you're simulating 1-500 possible
#Set the number of replicates to simulate
#Add the refuge location (grid cell ID)
#Add the user's initials to keep this file separate from others...
#Usage: Rscript simrep.R 1 10 5 JDR


library(holoSimCell)
source("../../../R/plothist.R")

args <- commandArgs(TRUE)
i <- as.numeric(args[1])
nreps <- as.numeric(args[2]) 
who <- as.character(args[3])  
#refs <- as.character(args[3])  #!!# If we want to have the refuge location passed as part of the command line 
#who <- as.character(args[4])

if(length(args) == 0){
    i <- 1
    nreps <- 1
    who <- "JDR"
    #refs <- 5
}

simdir <- system("echo $TMPDIR", intern = TRUE)
outdir <- "/mnt/research/TIMBER/Ash/OUT" 

library(holoSimCell)


### imputed, popmap (individualID->pop mapping), pts (sample locations) and ashpred
### are now built into holoSimCell
### as built in dataframes (in data/ directory)

rownames(popmap) <- popmap[,1]
table(popmap[gsub("fp","",names(imputed)),2])
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


    landscape <- ashpollen
    
    #Draw parms
    parms <- drawParms(control = system.file("extdata/csv","priors.csv",package="holoSimCell"))
    #parms$refs <- refs   #!!# If we want to have the refuge location passed as part of the command line
    parms$mu <- 1e-6  #!!!# Bumping up mutation rate!!  Does this help with speeds?

    if (FALSE) #don't run any without suitability
    {
                                        #Logical parameter of teh simulation, use hab_suit or not...
        if(parms$use.hab_suit == 0) {
                                        #landscape = NULL   #Don't do it this way, entire matrix is habitable
            landscape$hab_suit[!is.na(landscape$hab_suit)] <- 1  #This way ignores the glacier
                                        #landscape$hab_suit[landscape$hab_suit > 0] <- 1   #This way maintains the 0 suitability for glaciated cells
        }
    } else { #but set the model number to 2 for use.hab_suit
            parms$use.hab_suit = 2
    }

    
ph = getpophist.cells(#hab_suit=NULL,
                      hab_suit=landscape,
                      samptime=1,refs=parms$refs,refsz=parms$ref_Ne,
                      mix=parms$mix,
                      shortscale=parms$shortscale,shortshape=parms$shortshape,
                      longmean=parms$longmean,
                      sz=parms$sz,
                      K=parms$Ne,maxtime=701)

        plothist(ph,maxtime=NULL)
