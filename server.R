#server.R

shinyServer(function(input, output) {
  filedata <- reactive({
    infile <- input$datafile
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    read.csv(infile$datapath)
  })
  
  geodata <- reactive({
    if (input$getgeo == 0) return(NULL)
    df=filedata()
    if (is.null(df)) return(NULL)
    
    isolate({
      dummy=filedata()
      fr=input$from
      to=input$to
      locs=data.frame(place=unique(c(as.vector(dummy[[fr]]),as.vector(dummy[[to]]))),stringsAsFactors=F)      
      cbind(locs, t(sapply(locs$place,geocode, USE.NAMES=F))) 
    })
    
  })
  
  output$toCol <- renderUI({
    df <-filedata()
    if (is.null(df)) return(NULL)
    
    items=names(df)
    names(items)=items
    selectInput("to", "To:",items)
    
  })
  
  output$amountflag <- renderUI({
    df <-filedata()
    if (is.null(df)) return(NULL)
    
    checkboxInput("amountflag", "Use values?", FALSE)
  })
  
  output$fromCol <- renderUI({
    df <-filedata()
    if (is.null(df)) return(NULL)
    
    items=names(df)
    names(items)=items
    selectInput("from", "From:",items)
    
  })
  
  output$amountCol <- renderUI({
    df <-filedata()
    if (is.null(df)) return(NULL)
    #Let's only show numeric columns
    nums <- sapply(df, is.numeric)
    items=names(nums[nums])
    names(items)=items
    selectInput("amount", "Amount:",items)
  })
  
  output$lineSelector <- renderUI({
    radioButtons("lineselector", "Line type:",
                 c("Uniform" = "uniform",
                   "Thickness proportional" = "thickprop",
                   "Colour proportional" = "colprop"))
  })
  
  
  output$filetable <- renderTable({
    filedata()
  })
  
  
  geodata2 <- reactive({
    if (input$getgeo == 0) return(NULL)
    df=filedata()
    if (input$amountflag != 0) {
      maxval=max(df[input$amount],na.rm=T)
      minval=min(df[input$amount],na.rm=T)
      df$b8g43bds=10*df[input$amount]/maxval
    }
    gf=geodata()
    df=merge(df,gf,by.x=input$from,by.y='place')
    merge(df,gf,by.x=input$to,by.y='place')
  })
  
  output$geotable <- renderTable({
    if (input$getgeo == 0) return(NULL)
    geodata2()
  })
  
  plotter = reactive({
    if (input$getgeo == 0) return(map("world"))
    #Method pinched from: http://flowingdata.com/2011/05/11/how-to-map-connections-with-great-circles/
    map("world")
    df=geodata2()
    
    pal <- colorRampPalette(c("blue", "red"))
    colors <- pal(100)
    
    for (j in 1:nrow(df)){
      inter <- gcIntermediate(c(df[j,]$lon.x[[1]], df[j,]$lat.x[[1]]), c(df[j,]$lon.y[[1]], df[j,]$lat.y[[1]]), n=100, addStartEnd=TRUE)
      if (input$amountflag == 0) lines(inter, col="red", lwd=0.8)
      else {
        if (input$lineselector == 'colprop') {
          colindex <- round( (df[j,]$b8g43bds[[1]]/10) * length(colors) )
          lines(inter, col=colors[colindex], lwd=0.8)
        } else if (input$lineselector == 'thickprop') {
          lines(inter, col="red", lwd=df[j,]$b8g43bds[[1]])
        } else lines(inter, col="red", lwd=0.8)
      } 
    } 
  })
  
  output$geoplot<- renderPlot({
    if (input$getgeo == 0) return(map("world"))
    plotter()
  })
  
  #Not sure what's going on here; works fine with code inline but not if called as function?
  output$downloadPlot <- downloadHandler(
    filename = function() {paste('data-', Sys.Date(), '.png', sep='')},
    content = function(file) { 
      png(file)
      if (input$getgeo == 0) return(map("world"))
      #Method pinched from: http://flowingdata.com/2011/05/11/how-to-map-connections-with-great-circles/
      map("world")
      df=geodata2()
      
      pal <- colorRampPalette(c("blue", "red"))
      colors <- pal(100)
      
      for (j in 1:nrow(df)){
        inter <- gcIntermediate(c(df[j,]$lon.x[[1]], df[j,]$lat.x[[1]]), c(df[j,]$lon.y[[1]], df[j,]$lat.y[[1]]), n=100, addStartEnd=TRUE)
        if (input$amountflag == 0) lines(inter, col="red", lwd=0.8)
        else {
          if (input$lineselector == 'colprop') {
            colindex <- round( (df[j,]$b8g43bds[[1]]/10) * length(colors) )
            lines(inter, col=colors[colindex], lwd=0.8)
          } else if (input$lineselector == 'thickprop') {
            lines(inter, col="red", lwd=df[j,]$b8g43bds[[1]])
          } else lines(inter, col="red", lwd=0.8)
        }
      }
      dev.off()
    },
    contentType = 'image/png'
  )
  
})