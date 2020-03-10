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

samppts <- pts[pts$abbrev %in% names(poptbl),]



##this should produce a landscape with (x,y) _square_ cells that also have
##21empirical samples in separate grid cells (otherwise need to figure out something else)
## ashland is a stored R object as well
if (!exists("ashland"))
{
    ashland <- def_grid_pred(pred=ashpred[,,701:1],samppts=samppts,init.ext=c(40,36),keep.thresh=0.05)
}

landscape <- ashland

###seed is based on time in seconds and the number of characters in the library path
###
###
sec=as.numeric(Sys.time())-1500000000
lp= as.numeric(as.character(nchar(paste(.libPaths(), collapse = " "))))
slp <- as.integer(floor(sec*lp))

set.seed(as.integer(slp))

ph = getpophist2.cells(hab_suit=landscape,
                       refs=(5),
                       refsz=100,
                       mix=0.000005,  #note how small.
                       shortscale=6,  # scale parameter of weibull with shape below
                       shortshape=1, #weibull shape
                       longmean=75,  # mean of normal with sd = longmean
                       sz=150) #size of a cell (same units as longmean and shortscale)

gmap=make.gmap(ph$pophist,
               xnum=2,
               ynum=2)

ph2 <- pophist.aggregate(ph,gmap=gmap)


pdf("subset_example.pdf")
plothist(ph)
plothist(ph2)
dev.off()
