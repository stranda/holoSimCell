require(holoSimCell)
require(shiny)

# Define UI for application that draws a histogram

# Define server logic required to draw a histogram
server <- function(input, output,session) {

    habgrid <- reactive({
        if (input$usehab=="Pollen")
        {
            ashpollen
        } else {
            ashenm
        }})

    observeEvent(input$usehab,{
        hg <- habgrid()
        details <- hg$details
        updateNumericInput(session,"ref1",max=details[1,"ncells"])
        updateNumericInput(session,"ref2",max=details[1,"ncells"])
        updateNumericInput(session,"xdim",value=details[1,"x.dim"])
        updateNumericInput(session,"ydim",value=details[1,"y.dim"])
        updateNumericInput(session,"gens",value=dim(hg)[1])
    })
    
    output$dispkern <- renderPlot({
        thresh = 1e-9
        dom=0:100000
        dens <- distancePDF(dom,ssh=input$ssh,ssc=(input$ssc/100)*input$sz,
                            lmn=(input$nmn/100)*input$sz,
                            lsd=sqrt(input$nmn),mix=input$mix)
        nonzero.dens <- dens[1:max(which(dens>thresh))]
        nonzero.dom <- dom[1:max(which(dens>thresh))]
        prop.dom = nonzero.dom
        plot(nonzero.dens ~ prop.dom,
             type="l", xlab="Dispersal distance in same units as cell width",
             ylab="Density")
        text(0.9*max(nonzero.dom),0.9*max(nonzero.dens),paste("Only densities >",thresh,"shown"))
    })
    
    output$histplot <- renderPlot({
        input$doplot
        isolate({
            if ((input$ref1>0)&(input$ref2>0))
            {
                refs=c(input$ref1,input$ref2)
                refsz=c(input$ref1size,input$ref2size)
            } else {
                if (input$ref1>0)
                {
                    refs=c(input$ref1)
                    refsz=c(input$ref1size)
                } else {
                    refs=c(input$ref2)
                    refsz=c(input$ref2size)
                }
                
            }

            if (input$usehab!="None") hs=habgrid() else hs=NULL
            
            ph <- getpophist.cells(h=input$xdim*input$ydim,          ##demography, num habitats
                                         xdim=input$xdim,        ##num cols
                                         ydim=input$ydim,        ##num rows
                                         maxtime=input$gens,  ##num time clicks to simulate (yrs, decades cents?) (call em decades here)
                                         lambda=input$lambda, ##intrisic rate of growth (discrete-time...right? can be modified with deltLambda)
                                        #deltLambda: see function definition
                                         K=input$K,        ##carry capacity (can be modified with deltK)
                                        #deltK: see function definition
                                         refs=refs,    ##pop ids of refugia
                                         refsz=refsz, ##sizes of refugia indicated in refs
                                         sz=input$sz,   ##number of spatial units per grid cell
                                        #dispersal traits
                                         distance.fun=distancePDF, #takes a length and 5 parameters and creates mig matrix
                                         shortscale=(input$ssc/100)*input$sz, #scale of short distance weibull
                                         longmean=(input$nmn/100)*input$sz,      #mean of long-distance norm (var is the same as the mean)
                                         shortshape=input$ssh,    #shape of short-distance
                                         mix=input$mix,         #proportion of LDD
                                         popDispInfl=function(x){log(x+1)},
                                         samptime = 1,  #!# sample every X generations
                                         CVn = NULL,       #!# coefficient of variation X% 
                                         pois.var = FALSE, #!# Poisson distribution for demographic stochasticity
                                         extFUN = NULL,
                                         hab_suit = hs
                                       )
#            print(ph$coalhist)
            print("about to plot")
                if (sum(!is.na(ph$pophist$source))>1)
                    plothist(ph)
                else
                {
                    plot(1~1,xlab="",ylab="",type="n")
                    text(1,1,"No populations colonized!" )
                }
            print("plotting done")
            
        })
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

