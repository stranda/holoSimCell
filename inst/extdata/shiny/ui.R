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
                     "Scale of short-dist (% cell width)",
                     min = 0,
                     max = 50,
                     value = 5),
         sliderInput("ssh",
                     "Shape of short-dist",
                     min = 0.01,
                     max = 10,
                     value = 1),
         sliderInput("nmn",
                     "Scale of long-dist (% cell width)",
                     min = 0,
                     max = 400,
                     value = 35),
         sliderInput("mix",
                     "proportion of long-dist",
                     min = 0,
                     max = 0.2,
                     value = 0.2),
         numericInput("lambda","lambda",value=1.05,min=0.5,max=2.0),
         sliderInput("K",
                     "carry capacity",
                     min=10,max=5000,value=1000),
         
         radioButtons("refs","Refuge?",c("PA","GA","TX","ALL")),
         
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
#          h3("Current dispersal kernel"),
#          plotOutput("dispkern"),
          h3("Map of colonization history. Red oldest, white most recent"),
         plotOutput("histplot",height="900px")
      )
   )
)

