
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
          numericInput("lambda","lambda",value=1.03,min=0.5,max=1.5),

         sliderInput("K",
                     "carry capacity",
                     min=10,max=50000,value=500),
         sliderInput("ssc",
                     "Scale of short-dist",
                     min = 1,
                     max = 100,
                     value = 10),
         sliderInput("ssh",
                     "Shape of short-dist",
                     min = 0.01,
                     max = 10,
                     value = 0.6),
         sliderInput("nmn",
                     "Scale of long-dist",
                     min = 1,
                     max = 10000,
                     value = 150),
         sliderInput("mix",
                     "proportion of long-dist",
                     min = 0,
                     max = 0.2,
                     value = 0.00),
                 sliderInput("gens",
                     "Number of time clicks:",
                     min = 1,
                     max = 100000,
                     value = 500),

         numericInput("xdim","number of X grids",value=15,min=1,max=100),
         numericInput("ydim","number of Y grids",value=15,min=1,max=100),

         sliderInput("sz",
                     "Number of measurement units across each grid cell (km?)",
                     min = 1,
                     max = 1000,
                     value = 150),

         numericInput("ref","Refugium pop num",value=10),

         numericInput("refsize","Refugium pop size",value=100)
         
         
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
          h2("hotter colors are more recent events"),
         plotOutput("histplot",height="900px")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
    output$histplot <- renderPlot({
        input$doplot
        isolate({
            pops <- getpophist.cells(h=input$xdim*input$ydim,          ##demography, num habitats
                                     xdim=input$xdim,        ##num cols
                                     ydim=input$ydim,        ##num rows
                                     maxtime=input$gens,  ##num time clicks to simulate (yrs, decades cents?) (call em decades here)
                                     lambda=input$lambda, ##intrisic rate of growth (discrete-time...right? can be modified with deltLambda)
                                        #deltLambda: see function definition
                                     K=input$K,        ##carry capacity (can be modified with deltK)
                                        #deltK: see function definition
                                     refs=input$ref,    ##pop ids of refugia
                                     refsz=c(input$refsize), ##sizes of refugia indicated in refs
                                     sz=input$sz,   ##number of spatial units per grid cell
                                        #dispersal traits
                                     distance.fun=distancePDF, #takes a length and 5 parameters and creates mig matrix
                                     shortscale=input$ssc, #scale of short distance weibull
                                     longmean=input$nmn,      #mean of long-distance norm (var is the same as the mean)
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

