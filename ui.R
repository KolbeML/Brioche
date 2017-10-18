#ui.R
shinyUI(pageWithSidebar(
  headerPanel("Great Circle Map demo"),
  
  sidebarPanel(
    fileInput('datafile', 'Choose CSV file',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    uiOutput("fromCol"),
    uiOutput("toCol"),
    uiOutput("amountflag"),
    conditionalPanel(
      condition="input.amountflag==true",
      uiOutput("amountCol")
    ),
    conditionalPanel(
      condition="input.amountflag==true",
      uiOutput("lineSelector")
    ),
    actionButton("getgeo", "Get geodata"),
    downloadButton("downloadPlot", "Save map file")
    
  ),
  mainPanel(
    tableOutput("filetable"),
    tableOutput("geotable"),
    plotOutput("geoplot")
  )
))