require(shiny)
require(holoSimCell)
##
ui <- fluidPage(

    
   # Application title
    titlePanel("Forward-time simulation visualizer"),
    tabsetPanel(
        tabPanel("Setup and run simulation",fluid=T,
                 sidebarLayout(
                     sidebarPanel(
                         
                         actionButton("doplot","Run simulation and plot"),
                         sliderInput("ssc",
                                     "Scale of short-dist (% cell width)",
                                     min = 0,
                                     max = 50,
                                     value = 12),
                         sliderInput("ssh",
                                     "Shape of short-dist",
                                     min = 0.01,
                                     max = 10,
                                     value = 1),
                         sliderInput("nmn",
                                     "Scale of long-dist (% cell width)",
                                     min = 0,
                                     max = 400,
                                     value = 60),
                         sliderInput("mix",
                                     "proportion of long-dist",
                                     min = 0,
                                     max = 0.2,
                                     value = 0.02),
                         numericInput("lambda","lambda",value=1.05,min=0.5,max=2.0),
                         sliderInput("K",
                                     "carry capacity",
                                     min=10,max=5000,value=1000),
                         radioButtons("refs","Refuge?",c("PA","GA","TX","ALL"))
                     ),
                     
                     mainPanel(
                         h3("Map of colonization history. Red oldest, white most recent"),
                         plotOutput("histplot",height="900px")
                     )
                 )
                 ), #end of first tabPanel
        tabPanel("Examine history in detail",fluid=T,
                  sidebarLayout(
                     sidebarPanel(
                  
                         sliderInput("timeslice",
                                     "Time point in simulation to examine",
                                     min = 1,
                                     max = 701,
                                     value = 1),
                         sliderInput("window",
                                     "number of timeclicks to examine",
                                     min = 0,
                                     max = 701,
                                     value = 0),
                         ),
                     
                     mainPanel(
                         h3("Raster colors represent hab. suitability.  Arrows are colonization events"),
                         plotOutput("histsliceplot",height="900px")
                     )
                  )
                 )
    )
)

