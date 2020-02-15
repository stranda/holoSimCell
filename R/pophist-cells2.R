###
### dispersal kernel
###
#' simple distance density function for short and ldd
#'
#' returns a density for a vector of x's
#'
#' @export
distancePDF <- function(x, ssh=1,ssc=1,lmn=100,lsd=100,mix=0)
    {
        (1-mix)*dweibull(x,shape=ssh,scale=ssc)+mix*dnorm(x,mean=lmn,sd=lsd)
    }


###
### Cell-based version of getpophist
###  exploring the cost (loss of some realism) vs benefit (speed and temporal and spatial granularity)
###

#' cell-based getpophist
#'
#' creates a population history on a landscape
#'
#' @export

getpophist.cells <- function(h=225, #humber of habitats (populations)
                             xdim=15, #x -extent
                             ydim=15, #y-extent
                             maxtime=1000, 
                             distance.fun=distancePDF,   # unlike previous versions, needs the dispersal
                             shortscale=0.35,longmean=3, # parameters along with landscape
                             shortshape=1,mix=0.01,        #weibull shape for short distance, mix
                             lambda=1.01,  #underlying population growth rate
                             deltLambda=rep(0,h), #adjustment to pop growth per population
                             CVn=NULL, #!# coefficient of variation in population size (demographic stochasticity) 
                             pois.var = FALSE, #!# if TRUE, population size each generation is drawn from a Poisson
                             extFUN=NULL, #!# Function passed to extirpate() to implement local extinction due to demographic stochasticity - 2 arguments
                             ##hab_suit = NULL, #!# matrix of habitat suitability (0-1), npops rows and maxtime/samptime columns
                             hab_suit = NULL, #object returned from def_grid_pred
                             K=10000, #underlying carrying capacity per pop
                             deltK = rep(1,h), #!# per population adjustment to K
                             refs=c(2), #population number of the refuges
                             refsz=c(1000),#the size of each refuge
                             popDispInfl=function(x){return (x)}, #function that maps pop size to col ability
                             sz=1,  #size of each cell
                             enmstep=NULL, #vector of the number of time clicks for each habitat
                             samptime=0 #if zero runs through.  If >zero reports every samptime
                             
                             )
{

    if (FALSE)
    {
        
        h=540;
        xdim=20;
        ydim=27;
        maxtime=1000
        distance.fun=distancePDF;   # unlike previous versions, needs the dispersal
        shortscale=50;longmean=3; # parameters along with landscape
        shortshape=1;mix=0.001;        #weibull shape for short distance, mix
        lambda=1.01;  #underlying population growth rate
        deltLambda=rep(0,h); #adjustment to pop growth per population
        K=10000; #underlying carrying capacity per pop
        deltK = rep(1,h) #!# per population adjustment to K
        refs=c(1,2) #population number of the refuges
        refsz=c(1000,1000)#the size of each refuge
        sz=150
        samptime=5 #!# census the population every 5 generations
        popDispInfl=function(x){return (x)}
        pois.var = FALSE
        enmstep=NULL
        CVn=NULL #!# coefficient of variation in population size (demographic stochasticity)
        extFUN=NULL #!# Function passed to extirpate() to implement local extinction due to demographic stochasticity - 2 arguments

    }

    if (!is.null(hab_suit))
    {
        if ((ydim*xdim)!=dim(hab_suit$hab_suit)[2])
        {
            warning("x and y dimensions must be the same as the undelying hab_suit; changing them")
            
        }
        if (h!=dim(hab_suit$hab_suit)[2])
        {
            warning("h must be the same as the length of underlying hab_suit; changing")
        }
        xdim=hab_suit$details$x.dim[1]
        ydim=hab_suit$details$y.dim[1]
        h <- hab_suit$details$ncells[1]
        deltLambda=rep(0,h); #adjustment to pop growth per population
        deltK = rep(1,h) #!# per population adjustment to K
        if (is.null(enmstep)) {enmstep <- 1:dim(hab_suit$hab_suit)[1]; maxtime=length(enmstep)}
    }
    
    print(paste("starting...",date()))
    struct <- c(
        xdim=xdim,
        ydim=ydim,
        maxtime=maxtime,
        shortscale=shortscale,
        longmean=longmean,
        shortshape=shortshape,
        mix=mix,
                                        #deltLambda=deltLambda,
        K=K,
                                        #deltK=deltK,
        refs=refs,
        refsz=refsz,
        sz=sz,
        samptime=samptime,
        pois.var=pois.var)

    pops <- cbind(data.frame(pop=1:(xdim*ydim)),expand.grid(col=1:xdim,row=1:ydim),
                  arrive = NA, source = NA)  #!# adding in the arrive and source here, list element 1

    
    
    ##Set initial popsizevec
    Nvec <- rep(0,dim(pops)[1])
    npop <- length(Nvec)
        #if (length(refs)>0) Nvec[refs] <- refsz   #!# changing this to account for mismatches in length of refs and refsz
    
    #!# Account for refsz and refs not matching in length, fill in arrive times for refs
    if(length(refsz) == length(refs)) {
        Nvec[refs] <- refsz
        pops$arrive[refs] <- 0
    } else {
        Nvec[refs] <- rep(refsz[1],length(refs))    
        pops$arrive[refs] <- 0
    }

    #!# If samptime is switched on, we record the history of population size every 'samptime' generations
    if(samptime > 0) {
        Nvec_hist <- matrix(data = NA, nrow = h, ncol = length(seq(0,maxtime,samptime)))
        colnames(Nvec_hist) <- paste0("gen",seq(0,maxtime,samptime))
        Nvec_hist[,"gen0"] <- Nvec
    } else {
        Nvec_hist <- NULL
    }       


###create an object that can hold the entire pophistory
###has to be a flexible object, choosing 3d matrix
    colhist <- vector("list",maxtime)
#    colhist <- array(NA,dim=c(nrow(pops),2,maxtime))
#    dimnames(colhist) <- list(1:dim(colhist)[1],c("from","to"),1:dim(colhist)[3])
###this is supposed to be the among-population migration matrix
###should represent the prob of movement from the center of pop i to the center of pop j
    tmat <- integratedMigMat(landx=xdim,landy=ydim,xnum=5,ynum=5,ysz=sz,xsz=sz,
                             sshp=shortshape,ssc=shortscale,mix=mix,nmean=longmean,nvar=longmean^2)
        
    #pops$arrive <- NA  #!# both of these are now done in the definition of pops
    #pops$source <- NA

    if (!is.null(hab_suit)) enmcnt=0  #which enm hindcast are we using each time step? If any
    
    print(paste("setup done; simulating...",date()))

####added to record complete histories

    for (gen in 1:maxtime)
    {
#        print(paste("gen:",gen))
        src <- getsrcC(tmat,popDispInfl(Nvec))
        if (sum(src!=1:npop)!=0)  #some cells had movement
        {
            dispersed=T
            disp <- data.frame(src=src[src!=(1:npop)],snk=c(1:npop)[src!=(1:npop)])
            disp <- disp[Nvec[disp$snk]==0,]  #sink empty
            disp <- disp[Nvec[disp$src]>0,]   #source had some individuals
            
                                        #                dn <- disp$src
                                        #                names(dn) <- disp$snk
            
            
            pops$arrive[as.numeric(disp$snk)] <- gen
            pops$source[as.numeric(disp$snk)] <- disp$src
            Nvec[disp$snk] <- 1

#            print(paste("len disp$snk",length(disp$snk)))

            
            colhist[[gen]] <- pops[as.numeric(disp$snk),]
            
        } else {#end if cells had movement
            dispersed <- F
        }
        
                                        #if ((samptime>0)&((gen==1) |(!(gen %% samptime)))) #do something every 'samptime' time units
                                        #!# changed above logical slightly... fill in a column of Nvec_hist every samptime generations
        if(samptime > 0) {
            if(!(gen %% samptime)) 
            {
                                        #print(gen)
                                        #print(Nvec)
                Nvec_hist[,paste0("gen",gen)] <- Nvec
            }
        }
        
###!# Implementing habitat suitability, if it's passed as an argument to the function
        if(!is.null(hab_suit)) {
            if (gen %in% enmstep) {enmcnt <- enmcnt+1}
            dk <- hab_suit$hab_suit[enmcnt,]/max(hab_suit$hab_suit,na.rm=T)
            Nvec <- growpops(Nvec, lambda = lambda, deltLambda = deltLambda, K=K, deltK=dk,
                             CVn=CVn, pois.var = pois.var)   #!# If CVn is not null, implements a stochastic model of population growth (see below)
        } else {
            Nvec <- growpops(Nvec, lambda = lambda, deltLambda = deltLambda, K=K, deltK=deltK, CVn=CVn, pois.var = pois.var)
        }
        
                                        #!# Use 'extirpate' function to model extinction due to demographic stochasticity at small population sizes
        if(!is.null(extFUN)) {
            Nvec <- extirpate(Nvec, extFUN=extFUN)
        }
        
                                        #!# Deal with extinctions - only record most recent colonization event in pops
                                        #extinct <- which(Nvec == 0 & !is.na(pops$arrive)) 
        extinct <- which(Nvec == 0)
        
        if(length(extinct) > 0) {
            pops$arrive[extinct] <- NA
#            print(extinct)
#            print(dim(Nvec_hist))
            
            Nvec_hist[extinct,] <- 0   #!# This is new, may need to be eliminated  Do we zero population size after extinction (for previous generations?)
###pops$source[extinct] <- NA   #Turn this off to keep track of sources for populations that went extinct - they may have had a role in colonization prior to extinction!
        }
###Nvec <- growpops(Nvec,lambda=lambda,deltLambda=deltLambda,K=K,deltK=deltK)  #!# see above for growpops line

        if (dispersed)
        {
            p1 <- pops[as.numeric(disp$snk),]
            p2 <- cbind(p1,N=Nvec[p1$pop])#[Nvec[disp$snk]>0,]
            colhist[[gen]] <- p2[!is.na(p2$arrive),]
        } else    colhist[[gen]] <- NULL

    }
    print(paste("simulating done...",date()))
                                        #pops  #!# output is below

    out <- list(pophist = pops, Nvecs = Nvec_hist, tmat = tmat, struct = struct, hab_suit=hab_suit, coalhist=parseColhist(colhist))   #!# output is now a list
    
    return(out)
}



###
### grow a vector of population sizes
###
#!# modified this function to:
#!# i) take habitat suitability
#!# ii) zero change deltK = 1 (instead of 0)
#!# iii) allow for demographic stochasticity (poorly?)
growpops <- function(nv,      
                     lambda=1,  #underlying population growth rate
                     deltLambda=rep(0,length(nv)), #adjustment to pop growth per population
                     K=100, #underlying carrying capacity per pop
                     deltK = rep(1,length(nv)), #per population adjustment to K   #!# new zero-change point is deltK = 1 (for ease of habitat suitability implementation)
                     CVn=NULL,  #!# coefficient of variation in population size, if !NULL then demographic stochasticity is included
                     pois.var = FALSE   #!# if TRUE, expected population size is drawn from a Poisson distribution (for demographic stochasticity)
                     )

{
    rl <- lambda*(1+deltLambda) #make a vector of real lambdas base times the deltas
    #loK <- K*(1+(deltK))
    loK <- K*deltK  #!# slightly different, but eases habitat suitability implementation
    rl <- ((rl-1) * ((loK-nv)/loK))+1 #adjust that lambda based on K to approximate logistic growth
    rl[rl < 0] <- 0         #Trying to fix a problem with -Inf population sizes in unsuitable patches!!  JDR 3/23/19
    rl[is.na(rl)] <- 0
    new_nv <- round(nv*rl,2)

    if(pois.var == TRUE) {
        new_nv = rpois(length(new_nv), new_nv)
    }

    #!# new stuff down here - demographic stochasticity
    #!# modeled as a negative binomial (mean and size parameters)
    #!# mean is the Nvec entry for each grid cell
    #!# size is calculated to give coefficient of variation = CVn (*** some issues here ***)
    if(!is.null(CVn)) {
        size <- new_nv/(new_nv*((CVn/100)^2)-1)  #!# this is sometimes negative?  arbitrarily set to a value if negative
        size[size <= 0] <- 2
        out_nv <- rnbinom(n=length(nv), mu=new_nv, size=size)
    } else {
        out_nv <- new_nv
    }
    out_nv
}



####returns a vector of sources of colonization for
####each population.  Expects tmat to have dimensions: lenth(Nvec)xlength(Nvec)
getsrc <- function(tmat,Nvec)
    {
        sapply(1:dim(tmat)[1],function(i){
                rw = tmat[i,]*(Nvec)
                rw[i]=1-sum(rw[-i])
                rw[i] <- ifelse(rw[i]<0,0,rw[i]) 
                which(rmultinom(1,1,rw)>0)
                })
    }



###
### PARSE THE FULL POPLST TO MAKE pophist and coalhist
###
#' parse the full demographic history to make a pophist and a coalhist object
parsePopslst <- function(popslst,pops)
{
    pophist <- as.data.frame(do.call(rbind,lapply(popslst,function(p)
    {
        if (nrow(p$arrive)>0) ret <- p$arrive else ret <- data.frame(time=NA,source=NA)
        ret$col=p$col
        ret$row=p$row
        ret$pop=p$pop
        names(ret)[1] <- "arrive"
        ret[,c("pop","col","row","arrive","source")]
    })))

    coalhist <- as.data.frame(do.call(rbind,lapply(popslst,function(p)
    {
        if (nrow(p$arrive)>0) ret <- p$arrive else ret <- data.frame(time=NA,source=NA)
        ret$src=p$pop
        names(ret) <- c("time","snk","src")
        ret[!is.na(ret$time),c("time","src","snk")]
    })))
    list(pophist=pophist,coalhist=coalhist)
}



###
### Cell-based version of getpophist (better extinction records than above)
###  exploring the cost (loss of some realism) vs benefit (speed and temporal and spatial granularity)
###

#' cell-based getpophist updated
#'
#' creates a population history on a landscape
#'
#' @export

getpophist2.cells <- function(h=225, #humber of habitats (populations)
                             xdim=15, #x -extent
                             ydim=15, #y-extent
                             maxtime=1000, 
                             distance.fun=distancePDF,   # unlike previous versions, needs the dispersal
                             shortscale=0.35,longmean=3, # parameters along with landscape
                             shortshape=1,mix=0.01,        #weibull shape for short distance, mix
                             lambda=1.01,  #underlying population growth rate
                             deltLambda=rep(0,h), #adjustment to pop growth per population
                             CVn=NULL, #!# coefficient of variation in population size (demographic stochasticity) 
                             pois.var = FALSE, #!# if TRUE, population size each generation is drawn from a Poisson
                             extFUN=NULL, #!# Function passed to extirpate() to implement local extinction due to demographic stochasticity - 2 arguments
                             ##hab_suit = NULL, #!# matrix of habitat suitability (0-1), npops rows and maxtime/samptime columns
                             hab_suit = NULL, #object returned from def_grid_pred
                             K=10000, #underlying carrying capacity per pop
                             deltK = rep(1,h), #!# per population adjustment to K
                             refs=c(2), #population number of the refuges
                             refsz=c(1000),#the size of each refuge
                             popDispInfl=function(x){return (x)}, #function that maps pop size to col ability
                             sz=1,  #size of each cell
                             enmstep=NULL, #vector of the number of time clicks for each habitat
                             samptime=0 #if zero runs through.  If >zero reports every samptime
                             
                             )
{

    if (FALSE)
    {
        h=540;
        xdim=20;
        ydim=27;
        maxtime=1000
        distance.fun=distancePDF;   # unlike previous versions, needs the dispersal
        shortscale=50;longmean=3; # parameters along with landscape
        shortshape=1;mix=0.001;        #weibull shape for short distance, mix
        lambda=1.01;  #underlying population growth rate
        deltLambda=rep(0,h); #adjustment to pop growth per population
        K=10000; #underlying carrying capacity per pop
        deltK = rep(1,h) #!# per population adjustment to K
        refs=c(1,2) #population number of the refuges
        refsz=c(1000,1000)#the size of each refuge
        sz=150
        samptime=1 #!# census the population every 5 generations #DEPRECATED
        popDispInfl=function(x){return (x)}
        pois.var = FALSE
        enmstep=NULL
        CVn=NULL #!# coefficient of variation in population size (demographic stochasticity)
        extFUN=NULL #!# Function passed to extirpate() to implement local extinction due to demographic stochasticity - 2 arguments
    }

############DEPRECATING SAMPTIME
############# Setting Samptime to 1 here for all runs.  This means that census happens every
#####          time click and we have a full record of history
####
    samptime=1
    
    if (!is.null(hab_suit))
    {
        if ((ydim*xdim)!=dim(hab_suit$hab_suit)[2])
        {
            warning("x and y dimensions must be the same as the undelying hab_suit; changing them")
            
        }
        if (h!=dim(hab_suit$hab_suit)[2])
        {
            warning("h must be the same as the length of underlying hab_suit; changing")
        }
        xdim=hab_suit$details$x.dim[1]
        ydim=hab_suit$details$y.dim[1]
        h <- hab_suit$details$ncells[1]
        deltLambda=rep(0,h); #adjustment to pop growth per population
        deltK = rep(1,h) #!# per population adjustment to K
        if (is.null(enmstep)) {enmstep <- 1:dim(hab_suit$hab_suit)[1]; maxtime=length(enmstep)}
    }
    
    print(paste("starting...",date()))
    struct <- c(
        xdim=xdim,
        ydim=ydim,
        maxtime=maxtime,
        shortscale=shortscale,
        longmean=longmean,
        shortshape=shortshape,
        mix=mix,
                                        #deltLambda=deltLambda,
        K=K,
                                        #deltK=deltK,
        refs=refs,
        refsz=refsz,
        sz=sz,
        samptime=samptime,
        pois.var=pois.var)

    pops <- cbind(data.frame(pop=1:(xdim*ydim)),expand.grid(col=1:xdim,row=1:ydim),
                  arrive = NA, source = NA)  #!# adding in the arrive and source here, list element 1

    ##Set initial popsizevec
    Nvec <- rep(0,dim(pops)[1])
    npop <- length(Nvec)

    popslst <- vector("list",length=npop)
    for (i in 1:length(popslst))
    {
        popslst[[i]] <- list(pop=i,col=pops[i,"col"],row=pops[i,"row"],
                             arrive = data.frame(time=NA,source=NA)[-1,],  #times and sources of colonization
                             extinct= c(),                                 #times of extinction
                             N= rep(0,length(seq(0,maxtime,samptime))))    #popsize through time
        names(popslst[[i]]$N) <- paste0("gen", seq(0,maxtime,samptime))
    }
    
###if (length(refs)>0) Nvec[refs] <- refsz   #!# changing this to account for mismatches in length of refs and refsz
###refugia
    #!# Account for refsz and refs not matching in length, fill in arrive times for refs

    if(length(refsz) != length(refs))
    {
        refsz <- rep(refsz[1],length(refs))
    }
    
    Nvec[refs] <- refsz
    pops$arrive[refs] <- 0
    cnt <- 1

    for (l in refs)
    {
        popslst[[l]]$arrive <- rbind(popslst[[l]]$arrive,data.frame(time=0,source=NA))
        popslst[[i]]$N[1] <- refsz[cnt]
        cnt <- cnt+1
    }

    
###basically samptime will always be switched on and set to 1
###!# If samptime is switched on, we record the history of population size every 'samptime' generations
#    if(samptime > 0) {

    Nvec_hist <- matrix(data = NA, nrow = h, ncol = length(seq(0,maxtime,samptime)))
    colnames(Nvec_hist) <- paste0("gen",seq(0,maxtime,samptime))
    
    Nvec_hist[,"gen0"] <- Nvec
    
#    } else {
#        Nvec_hist <- NULL
#    }       


###create an object that can hold the entire pophistory
###has to be a flexible object, choosing 3d matrix
    colhist <- vector("list",maxtime)
#    colhist <- array(NA,dim=c(nrow(pops),2,maxtime))
#    dimnames(colhist) <- list(1:dim(colhist)[1],c("from","to"),1:dim(colhist)[3])
###this is supposed to be the among-population migration matrix
###should represent the prob of movement from the center of pop i to the center of pop j
    tmat <- integratedMigMat(landx=xdim,landy=ydim,xnum=5,ynum=5,ysz=sz,xsz=sz,
                             sshp=shortshape,ssc=shortscale,mix=mix,nmean=longmean,nvar=longmean^2)

    if (!is.null(hab_suit)) enmcnt=0  #which enm hindcast are we using each time step? If any
    
    print(paste("setup done; simulating...",date()))

####added to record complete histories

    for (gen in 1:maxtime)
    {
#        print(paste("gen:",gen))

        
        src <- getsrcC(tmat,popDispInfl(Nvec))

#        print(Nvec)
#        print(src)
        
        if (sum(src!=1:npop)!=0)  #some cells had movement
        {
            dispersed=T
            disp <- data.frame(src=src[src!=(1:npop)],snk=c(1:npop)[src!=(1:npop)])
            disp <- disp[Nvec[disp$snk]==0,]  #sink empty
            disp <- disp[Nvec[disp$src]>0,]   #source had some individuals
            
            pops$arrive[as.numeric(disp$snk)] <- gen
            pops$source[as.numeric(disp$snk)] <- disp$src
            Nvec[disp$snk] <- 1
            colhist[[gen]] <- pops[as.numeric(disp$snk),]
            
        } else {#end if cells had movement
            dispersed <- F
        }
        

        
###!# Implementing habitat suitability, if it's passed as an argument to the function
        if(!is.null(hab_suit)) {
            if (gen %in% enmstep) {enmcnt <- enmcnt+1}
            dk <- hab_suit$hab_suit[enmcnt,]/max(hab_suit$hab_suit,na.rm=T)
            Nvec <- growpops(Nvec, lambda = lambda, deltLambda = deltLambda,
                             K=K, deltK=dk,
                             CVn=CVn, pois.var = pois.var)   #!# If CVn is not null, implements a stochastic model of population growth (see below)
        } else {
            Nvec <- growpops(Nvec, lambda = lambda, deltLambda = deltLambda, K=K, deltK=deltK, CVn=CVn, pois.var = pois.var)
        }


####
#### Extinctions:  There are two sources of extinction possible right now.  Demographic stochasticity and poor habitat suitability.
####        
####   demographic stochasticity here refers to the degree to which probability of extinction in a cell is determined by population size.  The 
####       function 'extirpate' does the calculation but the actual functional relationship between pop size and extintion prob is encoded
####       in the function extFUN, passed as a parameter.  If NULL, there is no demographc stochasticity and only deterministic population growth
####       remains in place
####
####   poor habitat suitability can cause cells to go extinct.  The habitat suitability parameter and the population growth rate parameter determine
####       when and where this happens.

####   The code below, after possibly running extirpate, then checks for extinction by comparing the previous population size in a cell to the current
####       population size.  If the previous was greater than 0 and the current equals 0, and extinction has occurred.
####

        
        
###!# Use 'extirpate' function to model extinction due to demographic stochasticity at small population sizes
###        extFUN is a function passed as a getpophist parameter
        if(!is.null(extFUN)) {
            Nvec <- extirpate(Nvec, extFUN=extFUN)
        }
        
###!# Deal with extinctions - only record most recent colonization event in pops
####extinct <- which(Nvec == 0 & !is.na(pops$arrive)) 

        extinct <- which(Nvec_hist[,paste0("gen",gen-1)]>0 & Nvec==0)  # if it is zero now AND was greater than zero last gen
                
        if(length(extinct) > 0) {
            pops$arrive[extinct] <- NA
            for (pp in extinct) popslst[[pp]]$extinct <- c(popslst[[pp]]$extinct,gen)
#            print(Nvec)
#            print(extinct)
#            print(dim(Nvec_hist))
            
#            Nvec_hist[extinct,] <- 0   #!# This is new, may need to be eliminated  Do we zero population size after extinction (for previous generations?)
###pops$source[extinct] <- NA   #Turn this off to keep track of sources for populations that went extinct - they may have had a role in colonization prior to extinction!
        }
###Nvec <- growpops(Nvec,lambda=lambda,deltLambda=deltLambda,K=K,deltK=deltK)  #!# see above for growpops line

        if (dispersed)
        {
            p1 <- pops[as.numeric(disp$snk),]
            p2 <- cbind(p1,N=Nvec[p1$pop])#[Nvec[disp$snk]>0,]
            colhist[[gen]] <- p2[!is.na(p2$arrive),]

            for (i in 1:nrow(disp))
                if (Nvec[disp$snk[i]]>0)
                {
                    popslst[[disp$snk[i]]]$arrive <- rbind(popslst[[disp$snk[i]]]$arrive,
                                                           data.frame(time=gen,source=disp$src[i]))
                }

        } else    colhist[[gen]] <- NULL

###!# changed  slightly... fill in a column of Nvec_hist every  generation
        Nvec_hist[,paste0("gen",gen)] <- Nvec
        for (pp in 1:length(Nvec)) popslst[[pp]]$N[gen+1] <- Nvec[pp]
        
    } #end of simulation loop
    print(paste("simulating done...",date()))
                                        #pops  #!# output is below
    rl <- parsePopslst(popslst,pops)
    out <- list(pophist = rl$pophist, Nvecs = Nvec_hist, tmat = tmat, struct = struct, hab_suit=hab_suit, coalhist=rl$coalhist, popslst=popslst)   #!# output is now a list
    
    return(out)
}




















####implements a model of stochastic extinction
####functional form should be flexible, defined in terms of population size only
extirpate <- function(Nvec=NULL, extFUN=NULL)
    {
        ext_prob <- extFUN(Nvec)
        is_ext <- runif(length(Nvec)) < ext_prob
        Nvec[is_ext == TRUE] = 0
        Nvec
    }



#Other extFUN tried 
#ext^nvec



    
