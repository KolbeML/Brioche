library(shiny)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(leaflet)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "Basic dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Widgets", tabName = "widgets", icon = icon("th"))
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              fluidRow(
                box(plotOutput("plot1", height = 250)),
                box(plotOutput("plot2", height = 250)),
                
                box(
                  title = "Controls",
                  sliderInput("myslider", "Number of observations:", 1, 50, 15)
                )
              )
      ),
      
      # Second tab content
      tabItem(tabName = "widgets",
              h2("Widgets tab content")
      )
    )
  )
)

server <- function(input, output) {
  histdata <- rnorm(500)
  
  output$plot1 <- renderPlot({
    dat <- mtcars[1:input$myslider,]
    ggplot(dat, aes(mpg)) + geom_histogram(fill="cadetblue")
  })
  output$plot2 <- renderPlot({
    dat <- mtcars[1:input$myslider,]
    ggplot(dat, aes(mpg, wt)) + geom_point() + stat_smooth()
  })
}

shinyApp(ui, server)