library(holoSimCell)
library(ggplot2)

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



treats <- list()
treats[[1]] <- c(70,100,6)
treats[[2]] <- c(90,134,12)
treats[[3]] <- c(50,85,3)
treats[[4]] <- c(60,90,3)
treats[[5]] <- c(20,33,2)
treats[[6]] <- c(40,65,3)
treats[[7]] <- c(30,55,2)
treats[[8]] <- c(45,70,3)

res <- mclapply(treats,mc.cores=4,mc.preschedule=FALSE,function(trt)
#res <- lapply(treats,function(trt)
    {

        rs <- stack("study_region_dalton_ice_mask_lakes_masked_linearInterpolation.tif")
        e <- extent(rs)
        corners <- (matrix(c( e[1], e[4],
                             e[1], e[3],
                             e[2], e[3],
                             e[2],e[4]),ncol=2,byrow=T))
        colnames(corners) <- c("x","y")
        rownames(corners) <- c("ul","ll","lr","ur")
        
        icelakes <- as.array(1-rs) #really slow! thats one reason these layers get stored
        print(trt)
        icelakesland <- def_grid_pred(pred=icelakes,samppts=samppts,
                                      init.ext=c(trt[1],trt[2]),
                                      keep.thresh=0.05,corners=corners,base.ext="col"
                                      )
        
        landscape <- icelakesland 
        print(landscape$details)
###seed is based on time in seconds and the number of characters in the library path
###
###
        sec=as.numeric(Sys.time())-1500000000
        lp= as.numeric(as.character(nchar(paste(.libPaths(), collapse = " "))))
        slp <- as.integer(floor(sec*lp))
        
        set.seed(as.integer(sec))
        rows=landscape$details$x.dim
        cols=landscape$details$y.dim
        rf = round(rows*cols-cols*trt[3])
strt <- Sys.time()
        ph = getpophist2.cells(hab_suit=landscape,
                               refs=rf,
                               refsz=100,
                               mix=0.00,  #note how small.
                               shortscale=0.04,  # scale parameter of weibull with shape below
                               shortshape=1, #weibull shape
                               longmean=0.1,  # mean of normal with sd = longmean
                               sz=1) #size of a cell (same units as longmean and shortscale)

        c(dur=as.numeric(Sys.time())-as.numeric(strt),rows=rows,cols=cols,npops=landscape$details$npops,
          nempty=landscape$details$nempty,ref=rf)
    })
resdf <- as.data.frame(do.call(rbind,res))

save(file="resdata.rda",res,resdf)

resdf <- resdf[order(resdf$rows),]

ggplot(resdf,aes(y=dur,x=rows))+
    geom_point() + geom_line() +
    scale_y_continuous("sec", trans="log10")+
    xlab("Number of rows in demographic object")
