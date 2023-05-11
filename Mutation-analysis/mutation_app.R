library(shiny)
library(ggplot2)
library(dplyr)
library(shinydashboard)




mut.df <- data.frame("Position" = 1:1000,
)


ui <- fluidPage(
  actionButton(inputId = "upd", label= "Update"),
  selectInput(inputId = "column", label = "Select the data for analysis", 
              choices= colnames(mut.df), selected = "AA.pos"),
  hr(),
  column(6,
  checkboxGroupInput(inputId = "eff", label = "Effect",
               choices = c("Diagnosed" = "+",
                           "Unknown" = "?",
                           "Asymptomatic" =" -")),
  checkboxGroupInput(inputId = "exposure", label = "Exposed or Buried", 
               choices = c("Exposed", "Buried"), selected = "Either")),
  column(6, checkboxGroupInput(inputId = "eff2", label = "Effect2",
                     choices = c("Diagnosed" = "+",
                                 "Unknown" = "?",
                                 "Asymptomatic" =" -")),
  checkboxGroupInput(inputId = "exposure2", label = "Exposed or Buried 2", 
                     choices = c("Exposed", "Buried"), selected = "Either"),
  radioButtons(inputId = "binning", label = "Select x axis representation",
              choices = c("Position", "Region", "Domain"), selected= "Position"),
  selectInput(inputId = "col", label = "How do you want it colored?", 
              choices = c("Domain", "Region"), selected="Region")),
  column(6, plotOutput("pop1")),
  column(6, plotOutput("pop2")),
  verbatimTextOutput("f"),
  verbatimTextOutput("t"),
  verbatimTextOutput("jf")
)

server <- function(input, output) {
  observeEvent(input$upd,{
    rv <- if(is.null(input$eff)) mut.df
      else if(length(input$eff) == 1) {
        filter(mut.df, Effect == input$eff)
      }
      else if(length(input$eff) == 2) {
        filter(mut.df, Effect == input$eff[1] | Effect == input$eff[2])
      }
      else filter(mut.df, Effect == input$eff[1] | Effect == input$eff[2] | Effect == input$eff[3])
    
    pop1 <- if(is.null(input$exposure)) rv
    else if(length(input$exposure) == 1) filter(rv, Water == input$exposure)
    else filter(rv, Water == input$exposure[1] | Water == input$exposure[2])
    
    rv2 <- if(is.null(input$eff2)) mut.df
    else if(length(input$eff2) == 1) {
      filter(mut.df, Effect == input$eff2)
    }
    else if(length(input$eff2) == 2) {
      filter(mut.df, Effect == input$eff2[1] | Effect == input$eff2[2])
    }
    else filter(mut.df, Effect == input$eff2[1] | Effect == input$eff2[2] | Effect == input$eff2[3])
    
    pop2 <- if(is.null(input$exposure2)) rv2 
    else if(length(input$exposure2) == 1) filter(rv2, Water == input$exposure2)
    else filter(rv2, Water == input$exposure2[1] | Water == input$exposure2[2])
    
    output$f <- renderPrint({
      var.test(pop1[[input$column]], pop2[[input$column]], na.rm=TRUE)
    })
    t.logical <- ifelse(var.test(pop1[[input$column]], pop2[[input$column]], na.rm=TRUE)$p.value < 0.05,
                        TRUE, FALSE)
    output$t <- renderPrint({
      t.test(pop1[input$column], pop2[input$column], var.equal = t.logical)
    })
    output$pop1 <- renderPlot({
      ggplot(pop1, aes(x= pop1$AA.pos , y= pop1[input$column])) +
        geom_point()+
        labs(x = "Position", y = paste(input$column), title = paste(input$column, "vs. Position")) +
        theme_bw()})
    output$pop2 <- renderPlot({
      ggplot(pop2, aes(x=pop2$AA.pos, y= pop2[input$column])) +
        geom_point() +
        labs(x = "Position", y = paste(input$column), title = paste(input$column, "vs. Position")) +
        theme_bw()})
    output$jf <- renderPrint({
      input$col
     })
  })
}

shinyApp(ui = ui, server = server)
