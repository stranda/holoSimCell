
library(shiny)
#library(Rcpp)
#sourceCpp("helpers.cpp") #this is the c++ code, sourceCpp compiles it and makes functions available
#source("pophist-cells.R") #implements getpophist.cells
#source("integrate-mig-mat.R") # does a better job of specifying colonization probs, but still needs work
#source("plothist.R") #same ol' as always


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Forward-time simulation visualizer"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
 
          actionButton("doplot","Run simulation and plot"),
         sliderInput("ssc",
                     "Scale of short-dist (percent of cell width)",
                     min = 0,
                     max = 200,
                     value = 4),
         sliderInput("ssh",
                     "Shape of short-dist",
                     min = 0.01,
                     max = 10,
                     value = 1),
         sliderInput("nmn",
                     "Scale of long-dist (percent cell width)",
                     min = 0,
                     max = 400,
                     value = 35),
         sliderInput("mix",
                     "proportion of long-dist",
                     min = 0,
                     max = 0.1,
                     value = 0.001),
         sliderInput("sz",
                     "Number of measurement units across each grid cell (km?)",
                     min = 1,
                     max = 2000,
                     value = 150),
          numericInput("lambda","lambda",value=1.01,min=0.5,max=1.5),

         sliderInput("K",
                     "carry capacity",
                     min=10,max=50000,value=500),
         sliderInput("gens",
                     "Number of time clicks:",
                     min = 1,
                     max = 100000,
                     value = 500),

         numericInput("xdim","number of X grids",value=15,min=1,max=100),
         
         numericInput("ydim","number of Y grids",value=15,min=1,max=100),
         
         numericInput("ref1","Refugium pop num",value=10),

         numericInput("ref1size","Refugium pop size",value=100),
         
         numericInput("ref2","Refugium pop num",value=0),
         
         numericInput("ref2size","Refugium pop size",value=0)
         
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
          h3("Current dispersal kernel"),
          plotOutput("dispkern"),
          h3("Map of colonization history. Hotter colors are more recent"),
         plotOutput("histplot",height="900px")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   

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
            
                pops <- getpophist.cells(h=input$xdim*input$ydim,          ##demography, num habitats
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
                                         hab_suit = NULL
                                         )
                if (sum(!is.na(pops$pophist$source))>1)
                    plothist(pops$pophist)
                else
                {
                    plot(1~1,xlab="",ylab="",type="n")
                    text(1,1,"No populations colonized!" )
                }
            
        })
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

