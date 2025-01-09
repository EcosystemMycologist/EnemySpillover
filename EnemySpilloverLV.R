
library(shiny)
library(deSolve)

ui <- fluidPage(
   # Application title
   titlePanel("Lotka_Voltera model of Enemy Spillover"),

   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput(inputId = "KN",label = "Native species carrying capacity (K)",
                     min = 0,max = 1000,value = 900, step=10),  
         sliderInput(inputId = "rN",label = "Native species intrinic growth rate (r)",
                     min = 0,max = 0.1,value = 0.04, step=0.01),
         sliderInput(inputId = "KI",label = "Invasive species carrying capacity (K)",
                     min = 0,max = 1000,value = 450, step=10),  
         sliderInput(inputId = "rI",label = "Invasive species intrinic growth rate (r)",
                     min = 0,max = 0.1,value = 0.08, step=0.01),
         sliderInput(inputId = "C",label = "Native - invasive niche overlap (competition)",
                     min = 0,max = 1,value = 0, step=0.01),
         sliderInput(inputId =  "tMax",label = "How long to run simulation",
                     min = 0,max = 2000,value = 1600, step=25),
         sliderInput(inputId =  "tInvasion",label = "Time of invasion (must be < total time)",
                     min = 0,max = 2000,value = 800, step=25),
         sliderInput(inputId =  "aN",label = "Attack rate of enemies on natives",
                     min = 0,max = 0.01,value = 0.0005, step=0.0001),
         sliderInput(inputId =  "eN",label = "Efficiency of enemies on natives",
                     min = 0,max = 1,value = 0.1, step=0.005),
         sliderInput(inputId =  "aI",label = "Attack rate of enemies on invasives",
                     min = 0,max = 0.01,value = 0.0005, step=0.0001),
         sliderInput(inputId =  "eI",label = "Efficiency of enemies on invasives",
                     min = 0,max = 1,value = 0.1, step=0.005),
         sliderInput(inputId =  "s",label = "Death rate of enemies",
                     min = 0,max = 0.1,value = 0.02, step=0.005),
         sliderInput(inputId =  "P",label = "initial enemies (set to 0 to exclude)",
                     min = 0,max = 100,value = 1, step=1)
      ),
      mainPanel(
         plotOutput("lineplot")
      )
   )
)
 

shinyServer <- function(input, output) {
    
    lotka.volterra <- function(t, y, params){
      # create local variables to refer to the populations
      N <- y[1]
      I <- y[2]
      P <- y[3]
      # create local variables to refer to the different parameters
      rN <- params$rN
      KN <- params$KN
      rI <- params$rI
      KI <- params$KI
      C <- params$C
      aN <- params$aN
      eN <- params$eN
      aI <- params$aI
      eI <- params$eI
      s <- params$s
      # rates of change of the two populations
      dN.dt <- rN * N * (1 - (N + I*C) / KN) - aN * N * P
      dI.dt <- rI * I * (1 - (I + N*C) / KI) - aI * I * P
      dP.dt <- eN * aN * N * P + eI * aI * I * P - s * P 
      return(list(c(dN.dt, dI.dt, dP.dt)))
    }
    
   output$lineplot <- renderPlot({
     # determine the time scale to iterate the model experiments
     # time at which to introduce the invasive plant
     tMax <- input$tMax
     times <- seq(0, input$tInvasion, by=0.1)
     
     # set up the baseline parameters for the model.
     params<-list(rN = input$rN, 
                  KN = input$KN, 
                  rI = input$rI,
                  KI = input$KI,
                  C = input$C,
                  eN = input$eN,
                  aN = input$aN,
                  eI = input$eI,
                  aI = input$aI,
                  s = input$s)
     
     # run model pre-invasion
     y0 <- c(N=input$KN, I = 0, P = input$P)
     preInvasion <- ode(y0, times, lotka.volterra, params)
     # run moel post-invasion (starting from end point of previous)
     yI <- c(N=preInvasion[nrow(preInvasion),2], I = 1, P=preInvasion[nrow(preInvasion),4])
     timesI <- seq(input$tInvasion, tMax, by=0.1)
     postInvasion <- ode(yI, timesI, lotka.volterra, params)
     #combine two runs
     expt1 <- rbind(preInvasion, postInvasion)
     
     # plot result
     colors <- c("#78B389", "#2c2e9c", "#E64A00")
     matplot(
       expt1[,1],
       (expt1[, 2:ncol(expt1)]),
       type="l",
       xlab="Time",
       xlim=c(0,max(expt1[,1])),
       ylab="Population size",
       lwd=2, col=colors, lty=1, bty="n",
       ylim=c(0,max(expt1[,2:ncol(expt1)])+100)
     )
     abline(v=input$tInvasion, lty=2, col="grey")
     # add a legend to the figure
     legend("topright", c("Native plant","Invasive plant","Enemies"), col=colors, lty=1, lwd=2,bty="n")
     
      })
}
shinyApp(ui = ui, server = shinyServer)