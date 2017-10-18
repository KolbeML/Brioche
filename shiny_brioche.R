#import all packages necessary, apply applicable libraries
list.of.packages <- c("shiny", "shinythemes", "shinyjs", "ggplot2", "leaflet", "shinydashboard")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages,function(x){library(x,character.only=TRUE)})

#my variables - testy
myDF <- ""


ui <- dashboardPage(
  dashboardHeader(title = "Brioche"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("CSV Upload Test", tabName = "csvup", icon = icon("th")),
      menuItem("CSV Contents", tabName = "csvcontent", icon = icon("th")),
      menuItem("Contact Us", tabName = "contact", icon = icon("th"))
    )#end sidebarmenu
  ),#end dashboard sidebar
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
                ) #end box
              ) #end fluidrow
      ), #end dashboard tabitem
      
      # csv upload content
      tabItem(tabName = "csvup",
              fluidRow(
                fileInput("myFile", "Choose CSV File",
                        accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv") #end accept
                ) #end fileinput
              )#end fluidrow
      ),#end tabitem 
      
      tabItem(tabName = "csvcontent",
              fluidRow(
                tableOutput("contents")
              ) #end fluidrow
      ), #end dashboard tabitem
      
      tabItem(tabName = "contact",
              h2("Contact us or view Brioche's source code on ",
                a("GitHub", target="_blank", href="http://www.GitHub.com/KolbeML/Brioche/")
                  )
              ) #end tabitem contact
    )
  )
)



server <- function(input, output) {
  histdata <- rnorm(500)
  
  output$plot1 <- renderPlot({
    dat <- mtcars[1:input$myslider,]
    ggplot(dat, aes(mpg)) + geom_histogram(fill="cadetblue")
  })
  output$contents <- renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header = input$header)
  })
  output$plot2 <- renderPlot({
    dat <- mtcars[1:input$myslider,]
    ggplot(dat, aes(mpg, wt)) + geom_point() + stat_smooth()
  })              
}

shinyApp(ui, server)