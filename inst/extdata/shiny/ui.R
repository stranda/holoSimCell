require(shiny)
require(holoSimCell)
##
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
                     max = 20,
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
                     value = 20),
         sliderInput("mix",
                     "proportion of long-dist",
                     min = 0,
                     max = 0.2,
                     value = 0.001),
         sliderInput("sz",
                     "Number of measurement units across cells",
                     min = 1,
                     max = 2000,
                     value = 150),
         numericInput("lambda","lambda",value=1.01,min=0.5,max=2.0),

         sliderInput("K",
                     "carry capacity",
                     min=10,max=5000,value=500),
         sliderInput("gens",
                     "Number of time clicks:",
                     min = 1,
                     max = 10000,
                     value = 500),

         checkboxInput("sc2mxTime","Scale colors to full simulation period",FALSE),

         radioButtons("usehab","Which habitat suitability",c("None","ENM","Pollen","Workshop ENM","Workshop Pollen")),
         
         numericInput("xdim","number of X grids",value=15,min=1,max=100),
         
         numericInput("ydim","number of Y grids",value=15,min=1,max=100),
         
         numericInput("ref1","Refugium pop num",value=5),

         numericInput("ref1size","Refugium pop size",value=100),
         
         numericInput("ref2","Refugium pop num",value=0),
         
         numericInput("ref2size","Refugium pop size",value=0)
         
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
          h3("Current dispersal kernel"),
          plotOutput("dispkern"),
          h3("Map of colonization history. Red oldest, white most recent"),
         plotOutput("histplot",height="900px")
      )
   )
)

