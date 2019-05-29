#
# segmented regression
#
segmentGLM <- function(X,Y)
    {
        segment <- function(p,X,Y)
            {
                X1 <- X[X<=p]
                Y1 <- Y[X<=p]
                X2 <- X[X>p]
                Y2 <- Y[X>p]
                fita <- glm(Y1~X1)
                fitb <- glm(Y2~X2)
                logLik(fita)+logLik(fitb)
            }
        
        fit1 <- glm(Y~X)
        fit1_ll <- logLik(fit1)

        suppressWarnings(opt <- optimize(f=segment,X=X,Y=Y,lower=min(X,na.rm=T),upper=max(X,na.rm=T),maximum=T))
        c(breakpoint=opt$maximum,diffLL=opt$objective-fit1_ll)
    }

fake.data <- function()
    {
        X <- 1:100
        Y <- ifelse(X>20,10+20*1.4+X*0.4+rnorm(80,sd=3),
                      10+X*1.4+rnorm(20,sd=3)
            )
        data.frame(X=X,Y=Y)
    }
