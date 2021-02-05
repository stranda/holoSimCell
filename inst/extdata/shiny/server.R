require(holoSimCell)
require(shiny)

# Define UI for application that draws a histogram

# Define server logic required to draw a histogram
server <- function(input, output,session) {
    land <- reactive({
        if (FALSE)
            ashSetupLandscape(equalsuit=T,xlim=50,ylim=50,timesteps=701) else ashSetupLandscape()
    })
    getparms <- reactive({
        data.frame(shortscale=input$ssc, shortshape=input$ssh, longmean=input$nmn,
                   lambda=input$lambda,
                   mix=input$mix,
                   refs=input$refs,
                   Ne=input$K, ref_Ne=100*input$K)
    })

    ph <- reactive({
print("bout to call  landscape")
landscape <- land()
print("got past the landscape reactive")
        parms <- getparms()

print("made it past landscape and getparms")
        
        if(parms$refs == "PA") {
            refpops <- c(1000, 1001, 999, 1051, 949)
                                        #refpops <- 1000
        } else if(parms$refs == "TX") {
            refpops <- c(526, 527, 525, 577, 475)
                                        #refpops <- 526
        } else if(parms$refs == "GA") {
            refpops <- c(536, 537, 535, 587, 485)
                                        #refpops <- 536
        } else if(parms$refs == "ALL") {
            refpops <- c(536, 537, 535, 587, 485,
                         526, 527, 525, 577, 475,
                         1000, 1001, 999, 1051, 949)
                                        #refpops <- c(526,536,1000)
        }

print(parms)
print(refpops)

avgCellsz <- mean(c(res(landscape$sumrast)))  ### if the cells are not square, then use the average of the cell width and length
p <- getpophist2.cells(h = landscape$details$ncells,
                       xdim = landscape$details$x.dim,
                       ydim = landscape$details$y.dim,
                       hab_suit=landscape,
                       refs=refpops, 
                       refsz=parms$ref_Ne,
                       lambda=parms$lambda,
                       mix=parms$mix,
                       shortscale=(parms$shortscale/100)*avgCellsz,
                       shortshape=parms$shortshape,
                      longmean=(parms$longmean/100)*avgCellsz,
                       ysz=res(landscape$sumrast)[2], 
                       xsz=res(landscape$sumrast)[1], 
                      K = parms$Ne)
save(file="p.rda",p,parms,landscape)
p
    })

    
    output$histplot <- renderPlot({
        input$doplot
        isolate({
            ph1 <- ph()
            print(names(ph1))
            plothist(ph1,maxtime=nrow(land()$hab_suit))
            })
        })
   
}

# Run the application 
#shinyApp(ui = ui, server = server)

