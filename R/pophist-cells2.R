#' simple distance density function for short and ldd
#'
#' returns a density for a vector of x's
#'
#' @export
distancePDF <- function(x, ssh=1,ssc=1,lmn=100,lsd=100,mix=0)
    {
        (1-mix)*dweibull(x,shape=ssh,scale=ssc)+mix*dnorm(x,mean=lmn,sd=lsd)
    }


#' getpophist2.cells
#'
#' Run a forward demographic simulation that creates a recolonization history of a species across a grid landscape. 
#'
#' @param h integer defining the total number of cells in the landscape matrix.
#' @param xdim integer indicating the number of cells in the x-axis of the landscape grid.
#' @param ydim integer indicating the number of cells in the y-axis of the landscape grid.
#' @param maxtime the maximum number of generations in the forward demographic simulation.
#' @param distance.fun a function describing the dispersal kernel for the simulation.
#' @param shortscale scale parameter of the Weibull probability density function (i.e., short-distance component of species' dispersal kernel).
#' @param longmean mean of the Normal probability density function (i.e., long-distance component of species' dispersal kernel).
#' @param shortshape shape parameter of the Weibull probability density function (i.e., short-distance component of species' dispersal kernel).
#' @param mix proportion of long-distance dispersal events out of the total number of dispersal events. 
#' @param lambda population growth rate. It can be a fixed value or drawn from the prior distribution of the lambda parameter.
#' @param deltLambda a vector (one entry per cell) that allows cell-specific growth rates to be modeled.  The default (0) implies no variation in growth across populations, while 0.1 increases growth rate by 10 percent and -0.1 reduces growth rate by 10 percent.
#' @param CVn the coefficient of variation in population size, if demographic stochasticity is modeled.
#' @param pois.var if TRUE, population size in each time step is drawn from a Poisson distribution (to incorporate demographic stochasticity).
#' @param extFUN a function that describes the probability of local extinction as a function of cell-specific population size.
#' @param hab_suit list with the following elements:
#' \itemize{
#' \item{\code{details}} {Data frame with the spatial extent of the landscape grid}
#' \item{\code{occupied}} {Integer indicating the number of occupied cells in the landscape grid}
#' \item{\code{empty}} {Integer indicating the number of empty cells in the landscape grid}
#' \item{\code{sampled}} {Integer indicating the spatial location (i.e., cells) of the sampled populations in the landscape grid}
#' \item{\code{hab_suit}} {Matrix with species habitat suitability (ranging from zero to one) through time.}
#' \item{\code{sumrast}} {Raster that illustrates the summed habitat suitability in each cell across all time steps}
#' \item{\code{samplocsrast}} {Raster showing the locations of cells that correspond to sampled populations}
#' \item{\code{samplocs}} {Simple feature encoding spatial vector data related to the sampling populations and a data frame with metadata of the sampling populations (population id, number of individuals, type of spatial data and coordinates)}
#' \item{\code{sampdf}} {Data frame with the spatial location of the sampling populations in the landscape grid}
#' \item{\code{NAmask}} {A matrix used to mask cells that are unsuitable (e.g., covered by glaciers, out of the study region, in large lakes or oceans)}}
#' @param K maximum population size in each grid cell (i.e., carrying capacity) scaled by habitat suitability in the cell.
#' @param deltK per-population adjustment to carrying capacity of a cell.
#' @param refs vector of integers identifying the cells occupied at the beginning of the forward demographic simulations (i.e., species refugia).
#' @param refsz vector of integers indicating the effective population size for each cell occupied at the beginning of the forward demographic simulation. Each value is drawn from the prior distribution of the ancestral effective population size parameter.
#' @param popDispInfl a function that describes the influence of population size on the probability of serving as a source for colonists.
#' @param sz the size of each cell in the landscape.
#' @param ysz height of cells in the landscape grid.
#' @param xsz width of cells in the landscape grid.
#' @param enmstep specifies the number of generations over which each habitat suitability layer should be applied.
#' @param numcolonists the number of individuals that colonize a new population.
#'
#' @details
#' Conducts forward simulation of landscape colonization from a set of refugial populations.  Requires parameters related to refugial locations, maximum possible population size, dispersal, habitat suitability through time, and population growth rate.
#'
#' @return
#' Returns a pophist object
#' \itemize{
#' \item{\code{pophist}} {a data frame describing the population history (source of colonists, time of colonization, etc.).}
#' \item{\code{Nvecs}} {a matrix that records cell-specific population size at each time step in the forward demographic simulation.}
#' \item{\code{tmat}} {migration matrix showing probabilities of moving between populations in the landscape.}
#' \item{\code{strct}} {A numeric vector with information about spatial extent of the landscape grid: number of generations, dispersal parameters, population carrying capacity, refugia location and size.}
#' \item{\code{hab_suit}} {List of ten objects including information on the spatial extent of the landscape grid, number of occupied and empty cells, spatial location and metadata of sampling populations, and habitat suitability across the landscape grid}
#' \item{\code{coalhist}} {a data frame that recodes colonization history in pophist for coalescent simulation (e.g., backward in time).}
#' \item{\code{popslst}} {the colonization history of populations in the landscape in list format.  Records information on timing of extinction events, if extinction is allowed in the simulation.}
#' }
#'
#' @examples
#' library(holoSimCell)
#' parms <- drawParms(control = system.file("extdata/ashpaper","Ash_priors.csv",package="holoSimCell"))
#' modchoice <- 1
#' load(file=paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",pollenPulls[[modchoice]]$file))
#' refpops <- pollenPulls[[modchoice]]$refs
#' avgCellsz <- mean(c(res(landscape$sumrast)))
#'
#' ph = getpophist2.cells(h = landscape$details$ncells, xdim = landscape$details$x.dim, ydim = landscape$details$y.dim,
#'                        hab_suit=landscape,
#'                        refs=refpops,   
#'                        refsz=parms$ref_Ne,
#'                        lambda=parms$lambda,
#'                        mix=parms$mix,  
#'                        shortscale=parms$shortscale*avgCellsz,  
#'                        shortshape=parms$shortshape, 
#'                        longmean=parms$longmean*avgCellsz,  
#'                        ysz=res(landscape$sumrast)[2], 
#'                        xsz=res(landscape$sumrast)[1], 
#'                        K = parms$Ne) 
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
                             ysz=sz,
                             xsz=sz,
                             enmstep=NULL, #vector of the number of time clicks for each habitat
                             numcolonists=1 #number of colonists in a newly colonized cell
                             )
{

    samptime=1

    hab_suit$hab_suit[K*(hab_suit$hab_suit/max(hab_suit$hab_suit, na.rm = TRUE)) < 1] <- 0
    
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
        deltLambda=rep(0,h); 
        deltK = rep(1,h) 
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
        K=K,
        refsz=refsz,
        xsz=xsz,
        ysz=ysz,
        samptime=samptime,
        pois.var=pois.var)

    pops <- cbind(data.frame(pop=1:(xdim*ydim)),expand.grid(col=1:xdim,row=1:ydim),
                  arrive = NA, source = NA)  

    Nvec <- rep(0,dim(pops)[1])
    npop <- length(Nvec)

    popslst <- vector("list",length=npop)
    for (i in 1:length(popslst))
    {
        popslst[[i]] <- list(pop=i,col=pops[i,"col"],row=pops[i,"row"],
                             arrive = data.frame(time=NA,source=NA)[-1,],  
                             extinct= c(),                                 
                             N= rep(0,length(seq(0,maxtime,samptime))))    
        names(popslst[[i]]$N) <- paste0("gen", seq(0,maxtime,samptime))
    }
    
    if(length(refsz) != length(refs))
    {
        refsz <- rep(refsz[1],length(refs))
    }

    Nvec[refs] <- refsz * (hab_suit$hab_suit[1,refs]/max(hab_suit$hab_suit,na.rm=T))
    
    pops$arrive[refs] <- 0
    cnt <- 1

    for (l in refs)
    {
        popslst[[l]]$arrive <- rbind(popslst[[l]]$arrive,data.frame(time=0,source=NA))
        popslst[[i]]$N[1] <- refsz[cnt]
        cnt <- cnt+1
    }

    Nvec_hist <- matrix(data = NA, nrow = h, ncol = length(seq(0,maxtime,samptime)))
    colnames(Nvec_hist) <- paste0("gen",seq(0,maxtime,samptime))
    
    Nvec_hist[,"gen0"] <- Nvec
    
    colhist <- vector("list",maxtime)

    tmat <- integratedMigMat(landx=xdim,landy=ydim,xnum=5,ynum=5,ysz=ysz,xsz=xsz,
                             sshp=shortshape,ssc=shortscale,mix=mix,nmean=longmean,nvar=longmean^2)

    if (!is.null(hab_suit)) enmcnt=0  
    
    print(paste("setup done; simulating...",date()))

    for (gen in 1:maxtime)
    {

        src <- getsrcC(tmat,popDispInfl(Nvec))

        if (sum(src!=(1:npop))!=0)  
        {
            dispersed=T
            disp <- data.frame(src=src[src!=(1:npop)],snk=c(1:npop)[src!=(1:npop)])
            disp <- disp[Nvec[disp$snk]==0,]  
            disp <- disp[Nvec[disp$src]>0,]   
            
            pops$arrive[as.numeric(disp$snk)] <- gen
            pops$source[as.numeric(disp$snk)] <- disp$src
            Nvec[disp$snk] <- numcolonists
            colhist[[gen]] <- pops[as.numeric(disp$snk),]
            
        } else {
            dispersed <- F
        }
        
        if(!is.null(hab_suit)) {
            if (gen %in% enmstep) {enmcnt <- enmcnt+1}
            dk <- hab_suit$hab_suit[enmcnt,]/max(hab_suit$hab_suit,na.rm=T)
            Nvec <- growpops(Nvec, lambda = lambda, deltLambda = deltLambda,
                             K=K, deltK=dk,
                             CVn=CVn, pois.var = pois.var)   
        } else {
            Nvec <- growpops(Nvec, lambda = lambda, deltLambda = deltLambda, K=K, deltK=deltK, CVn=CVn, pois.var = pois.var)
        }

        Nvec[!is.finite(Nvec)] <- 0 

        if (sum(!is.finite(Nvec))>0)
        
        {
            print(paste("gen:",gen))
            print(dk)
            print(Nvec)
        }

        if(!is.null(extFUN)) {
            Nvec <- extirpate(Nvec, extFUN=extFUN)
        }

        extinct <- which(Nvec_hist[,paste0("gen",gen-1)]>0 & Nvec==0)  
                
        if(length(extinct) > 0) {
            pops$arrive[extinct] <- NA
            for (pp in extinct) popslst[[pp]]$extinct <- c(popslst[[pp]]$extinct,gen)
     }


        if (dispersed)
        {
            p1 <- pops[as.numeric(disp$snk),]
            p2 <- cbind(p1,N=Nvec[p1$pop])
            colhist[[gen]] <- p2[!is.na(p2$arrive),]

            if (nrow(disp)>0)
                for (i in 1:nrow(disp))
                {
                    if (Nvec[disp$snk[i]]>0)
                    {
                        popslst[[disp$snk[i]]]$arrive <- rbind(popslst[[disp$snk[i]]]$arrive,
                                                               data.frame(time=gen,source=disp$src[i]))
                    }
                }

        } else    colhist[[gen]] <- NULL

        Nvec_hist[,paste0("gen",gen)] <- Nvec
        for (pp in 1:length(Nvec)) popslst[[pp]]$N[gen+1] <- Nvec[pp]
        
    } 
    print(paste("simulating done...",date()))
                                        
    rl <- parsePopslst(popslst,pops)
    tmpcoord <- coordinates(hab_suit$sumrast)
    popcoords <- data.frame(pop = c(1:length(tmpcoord[,1])), longitude = tmpcoord[,1], latitude = tmpcoord[,2])
    rl$pophist <- merge(rl$pophist, popcoords)
    out <- list(pophist = rl$pophist, Nvecs = Nvec_hist, tmat = tmat,
                struct = struct, hab_suit=hab_suit, coalhist=rl$coalhist,
                popslst=popslst)   
    
    return(out)
}



###
### grow a vector of population sizes
###
growpops <- function(nv,      
                     lambda=1,  #underlying population growth rate
                     deltLambda=rep(0,length(nv)), #adjustment to pop growth per population
                     K=100, #underlying carrying capacity per pop
                     deltK = rep(1,length(nv)), #per population adjustment to K   #!# new zero-change point is deltK = 1 (for ease of habitat suitability implementation)
                     CVn=NULL,  #!# coefficient of variation in population size, if !NULL then demographic stochasticity is included
                     pois.var = FALSE   #!# if TRUE, expected population size is drawn from a Poisson distribution (for demographic stochasticity)
                     )

{
    rl <- lambda*(1+deltLambda) 
    
    loK <- K*deltK  
    rl <- ((rl-1) * ((loK-nv)/loK))+1 
    rl[rl < 0] <- 0        
    rl[is.na(rl)] <- 0
    new_nv <- nv*rl

    if(pois.var == TRUE) {
        new_nv = rpois(length(new_nv), new_nv)
    }

    if(!is.null(CVn)) {
        size <- new_nv/(new_nv*((CVn/100)^2)-1)  
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
### parse the full demographic history to make a pophist and a coalhist object
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



####implements a model of stochastic extinction
####functional form should be flexible, defined in terms of population size only
extirpate <- function(Nvec=NULL, extFUN=NULL)
    {
        ext_prob <- extFUN(Nvec)
        is_ext <- runif(length(Nvec)) < ext_prob
        Nvec[is_ext == TRUE] = 0
        Nvec
    }






    


