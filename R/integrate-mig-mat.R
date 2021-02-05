##
## calculate a "more correct" dispersal rate into rectangular landscape 
##

#' calc dispersal into a cell
#'
#' calculates the amount of dispersal into a cell based on the popsize and distance of other cells
#'
#' returns a matrix with the centermost value corresponding to the source and the surrounding cells represent
#' the probability of dispersal into those cells.  Essentially this produces the dispersal kernel for a single
#' source population.  The migration matrix could be created by moving the center of this functions output to 
#' each population to determine dispersal into the others.
#' 
#' @export
integratedDispMat <- function(xnum=4,ynum=4,xsz=100,ysz=100,sshp=1,ssc=10,mix=0.1,nmean=100,nvar=(nmean))
{

    if (FALSE)
    {
       xnum=10;ynum=10;xsz=100;ysz=100;sshp=1;ssc=60;mix=0.05;nmean=100;nvar=sqrt(nmean )
    }

    topolar <- function(x,y)
    {
        r <- sqrt(x^2+y^2)
        theta <- atan(y/x)
        cbind(theta,r)
    }
    
    
  ynum1 = ynum+1
  xnum1 = xnum+1
  
  qmat = matrix(0,nrow = ynum1, ncol = xnum1)

 
  for (x in 1:xnum1)
  {
      if (x==1) lft = xsz/2 else lft = (x-1)*xsz
      rgt = lft+xsz
      for (y in 1:ynum1)
      {
         if (y==1) bt=ysz/2 else  bt = (y-1)*ysz
          tp = bt+ysz
#          print(paste(lft,rgt))
#          print(paste(bt,tp))
          if (FALSE)
          {####this is wrong, need to work in polar coordinates and integrate under a surface of rotation
              qmat[y,x] <- cubature::adaptIntegrate(
                                         function(l)
                                         {
                                             (1-mix)*dweibull(l[1],shape=sshp,scale=ssc)+mix*dnorm(l[1],mean=nmean,sd=sqrt(nvar))+
                                                 (1-mix)*dweibull(l[2],shape=sshp,scale=ssc)+mix*dnorm(l[2],mean=nmean,sd=sqrt(nvar))
                                         },
                                         c(lft,bt),c(rgt,tp))$integral
          } else {##integrating under a surface of rotation.  effen polar coordinates! This works AES 2/2/21
              f=function(x,y)
              {
                  d=sqrt(x^2+y^2)
                  (1-mix)*dweibull(d,shape=sshp,scale=ssc) + mix*dnorm(d,mean=nmean,sd=sqrt(nvar))
              }
              
              polar <- topolar(x=c(lft,rgt,lft,rgt),y=c(tp,tp,bt,bt))
              qmat[y,x] <- pracma::integral2(f,min(polar[,1],na.rm=T),max(polar[,1],na.rm=T),min(polar[,2],na.rm=T),
                                             max(polar[,2],na.rm=T),sector=T)$Q
          }
      }
  }
    
  qmat[1,1]=4*qmat[1,1]  #becaus we only integrated across 1/4 of central cell
  
  #make the single quadrant reflect across a 2d surface
  fmat <- matrix(0,nrow=1+ynum*2,ncol=1+xnum*2)
  fmat[(ynum)+1,] <-qmat[1,c(ynum1:1,2:ynum1)]
  fmat[,(xnum)+1] <-qmat[c(xnum1:1,2:xnum1),1]
  fmat[1:ynum,1:xnum] = qmat[ynum1:2,xnum1:2]
  fmat[1:ynum,(xnum+2):(2*xnum+1)] = qmat[ynum1:2,2:xnum1]
  fmat[(ynum+2):(ynum*2+1),] = fmat[ynum:1,]

 #   fmat[floor(dim(fmat)[1]/2),] <-   fmat[floor(dim(fmat)[1]/2),] / 2
 #   fmat[,floor(dim(fmat)[1]/2)] <-   fmat[,floor(dim(fmat)[1]/2)] / 2
  
  fmat
}


##
## calculate a "more correct" dispersal rate into rectangular landscape 
##

#' calc dispersal into a cell
#'
#' calculates the amount of dispersal into a cell based on the popsize and distance of other cells
#'
#' returns a matrix with the centermost value corresponding to the source and the surrounding cells represent
#' the probability of dispersal into those cells.  Essentially this produces the dispersal kernel for a single
#' source population.  The migration matrix could be created by moving the center of this functions output to 
#' each population to determine dispersal into the others.
#' 
#' @export
integratedMigMat <- function(landx=15,landy=15,xnum=4,ynum=4,xsz=100,ysz=100,sshp,ssc,mix,nmean,nvar)
{
    if (FALSE)
    {
        landx=10;
        landy=10;
        xnum=3;
        ynum=3;
        xsz=102;
        ysz=102;
        sshp=1;
        ssc=10;
        mix=0.001;
        nmean=70;
        nvar=(nmean)*10;
        dom=seq(0,2*xsz,length=100)
        plot(dnorm(dom,nmean,sqrt(nvar))*(mix)+(1-mix)*dweibull(dom,shape=sshp,scale=ssc)~dom,type="l")
    }
    
    dm = integratedDispMat(xnum=xnum,ynum=ynum,xsz=xsz,ysz=ysz,sshp=sshp,ssc=ssc,mix=mix,nmean=nmean,nvar=nvar)
    dm=dm/sum(dm)   
    dmcol <- dim(dm)[2];     dmrow <- dim(dm)[1]

    cent=c(ceiling(dim(dm)[1]/2),ceiling(dim(dm)[2]/2))
    prop.tail.out.cent=1-dm[cent[1],cent[2]]
    print("prop.tail.out.cent")
    print(prop.tail.out.cent)
    lmat <- matrix(0,ncol=landx,nrow=landy)
   
    pops <- cbind(data.frame(pop=1:(landx*landy)),expand.grid(x=1:landx,y=1:landy))

    if (FALSE)
        {
            rmat=makeRmat(pops,dm,cent)             #rversion fine for small landscapes
        } else {
            rmat= makeRmatC(as.matrix(pops),dm,(cent-1)) #c++ version better (code in helpers.cpp)
        }
#    rmat=rmat/sum(rmat)
    rmat   
}

makeRmat <- function(pops,dm,cent=c(0,0))
{
    rmat <- matrix(0,ncol=dim(pops)[1],nrow=dim(pops)[1]) #from cols to rows transition/migration matrix
    for (i in 1:dim(pops)[1])
    {
        for (j in 1:dim(pops)[1])
        {
            offx=pops[i,"x"]-pops[j,"x"]
            offy=pops[i,"y"]-pops[j,"y"]
            
            {
                if (((offx+cent[2])>0)&((offx+cent[2])<dim(dm)[2]))     #make sure you are not running off end of matrix
                    if (((offy+cent[1])>0)&((offy+cent[1])<dim(dm)[1]))
                        rmat[j,i] <- dm[offx+cent[2],offy+cent[1]]
            }
        }
    }
    rmat
}

makeRmatBetter <- function(pops,dm,cent=c(0,0))
{

    dimnames(dm) <- dimnames(dm) <- list(1:nrow(dm),1:ncol(dm))
    dmd=data.frame(cellx=colnames(dm)[col(dm)], celly=rownames(dm)[row(dm)], dens=c(dm))

    m=as.matrix(dist(pops[,2:3],diag=T,upper=T))

    
    d=data.frame(pop1=colnames(m)[col(m)], pop2=rownames(m)[row(m)], dist=c(m))
    d$pop1 <- as.character(d$pop1)
    d$pop2 <- as.character(d$pop2)
#    d=d[d$dist<((nrow(dm)-1)/2),]
    d$offy=pops[d$pop1,"y"]-pops[d$pop2,"y"]+cent[1]
    d$offx=pops[d$pop1,"x"]-pops[d$pop2,"x"]+cent[2]
    d <- d[(d$offx>0)&(d$offx<=ncol(dm)),]
    d <- d[(d$offy>0)&(d$offy<=nrow(dm)),]

    d <- merge(d,dmd,all.x=T,by.x=c("offx","offy"),by.y=c("cellx","celly"))
    rdf=merge(expand.grid(pop1=1:max(as.numeric(d$pop1)),pop2=1:max(as.numeric(d$pop2))),d,all.x=T)
    rdf=rdf[order(as.numeric(rdf$pop1),as.numeric(rdf$pop2)),]

    rmat=matrix(rdf$dens,nrow=length(unique(rdf$pop1)),ncol=length(unique(rdf$pop2)))
    rmat[is.na(rmat)]=0
    rmat
}

makeRmatVec <- function()  #needs parameters to work, mainly here for archival code purposes

    {
                    print("this far")
#            dmdf <- cbind(expand.grid(from=1:dmcol,to=1:dmrow),dm=c(dm))
            comps <- expand.grid(to=pops$pop,from=pops$pop)
            print("this far1")
            inter <- merge(comps,pops,by.x="to",by.y="pop")
            print("this far2")
            names(inter)[3:4] <- c("to.x","to.y")
            comps <- merge(inter,pops,by.x="from",by.y="pop")
            print("this far3")
            names(comps)[5:6] <- c("from.x","from.y")
            comps <- comps[order(comps$from,comps$to),]

            comps$ox <- comps$from.x-comps$to.x
            comps$oy <- comps$from.y-comps$to.y
            
print("this far4")
            
            comps$ox[(comps$ox+cent[2]<0)|(comps$ox+cent[2]>dmcol)] <- NA
            comps$oy[(comps$oy+cent[1]<0)|(comps$oy+cent[1]>dmrow)] <- NA

            rv  <- dm[cbind(comps$ox+cent[2],comps$oy+cent[1])]
            rv[is.na(rv)] <- 0

            rmat <- matrix(rv,ncol=landx*landy)

        }
