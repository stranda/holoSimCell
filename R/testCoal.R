###
### function to return true or false if there is a (easy to identify) reason for
### non-coalescence of a coalhist object
###
### This function results from a bug JDR identified in which there is incomplete fusion back in time for
### some populations, particularly those that are early in the simulation
### In this situation, multiple coalescences from different sources into the same sink occur at the
### same time step
### 
testCoal <- function(ch)
{

    ok <- TRUE

    ### test for the multiple coalesce problem
    cho <- ch[order(ch$src,ch$time,ch$snk),]
    coalProps <- sapply(unique(cho$src),function(src)
    {
        tmp <- cho[cho$src==src,]
        if (min(tmp$time)>0)
        {
            mt <- which(tmp$time==min(tmp$time))
            if (length(mt)>1)
            {
                prp <- sum(tmp[mt,"prop"])
                if ((min(tmp[mt,"prop"])==0)&(max(tmp[mt,"prop"])==1)&(prp==1))
                    1 else 0
            } else tmp$prop[mt]
        }
        else 1
    })
    if (min(coalProps)<1) ok <- FALSE

    ### next test when ready
    
    ok
}
