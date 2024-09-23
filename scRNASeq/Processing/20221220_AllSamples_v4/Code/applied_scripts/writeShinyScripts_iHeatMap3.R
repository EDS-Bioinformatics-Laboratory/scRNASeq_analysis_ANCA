# Functions to write the 'server.R' and 'ui.R' files for the Shiny application for Microarray and RNASeq analyses

# Microarray
#   default objects:
#     - tt_all.rds: containing the topTable information for all genes and comparisons
#     - eset.rds: ESet containing the (normalized) expression data and experiment information
#     - BroadSets.rds: containing MSigDB genesets
#     - tabCameraFile: prefix of the CAMERA table files holding the topTable of the genesets (fx. "tab.camera.v5.0_")
#     - The existence of 'aQM_before_neqc/index.html' and 'aQM_neqc/index.html' is implied
#     - The existence of a directory 'venn' holding the results of the makeClickableVenn function is implied

# RNASeq
#   default objects:
#     - tt_all.rds: containing the topTable information for all genes and comparisons
#     - v.rds : Voom object with added exp. information added to the v$targets slot
#     - BroadSets.rds: containing MSigDB genesets
#     - tabCameraFile: prefix of the CAMERA table files holding the topTable of the genesets (fx. "tab.camera.v5.0_")
#     - The existence of a directory 'venn' holding the results of the makeClickableVenn function is implied


# TO DO (??):
#   + Now, the scripts needs RDS objects, make it suitable for .RData ??
#   + Check exsistence of objects before reading them in?
#   + Check whether all information is there?
#   + Now written for human & voom (RNASeq), but what if mouse & edgeR? (fx 'mgi_symbol' vs. 'hgnc_symbol', D$counts vs. v$E)
#   - Build in all kinds of information into the UI? Like "Human", "Voom", platform etc.? -> it becomes big and intractable?
#   + Build in possibiliy to provide an 'ExpInfo.txt' file? (if it is given, create a tab for it, otherwise not)

#### General runApp script ####
write_runShinyApp <- function(fileOut="runShinyApp.R", libLoc="RLibs"){
  # Open the file to write to
  fileConn <- file(fileOut)
  
  # Gather the txt
  txt <- paste0("
    # Script to start Shiny app
    # - Detects whether we are on the CDW
    # - Creates a local directory for the R libraries to be installed
    # - Makes sure to use the Chrome browser if detected
    
    # https://www.r-bloggers.com/deploying-desktop-apps-with-r/
    dir.create(\"",libLoc,"\")
    .libPaths(\"",libLoc,"\")
    
    # Check availability of BiocManager package
    if (!requireNamespace(\"BiocManager\") && tools:::.BioC_version_associated_with_R_version() >= 3.8 ){
      install.packages(\"BiocManager\")
    }
    
    #
    if (!require(\"shiny\")){ 
      if ( version$major >= 4 || version$minor > 5.1 ){
        BiocManager::install(\"shiny\")
      } else {
        source(\"http://bioconductor.org/biocLite.R\")
        biocLite(\"shiny\", suppressUpdates=TRUE) 
      }
    }
    
    # Check where we are, i.e. local or CDW
    checkLocation <- system(\"WHERE /F /R C:\\\\ DirectFlex.exe\", intern=TRUE)
    if ( is.null(attributes(checkLocation)) ){
      # We are on the CDW...
      myBrowser <- checkLocation
    } else {
      #Set default browser to Chrome
      myBrowser <- system(\"WHERE /F /R C:\\\\ chrome.exe\", intern=TRUE)
    }
    
    options(browser = myBrowser)
    
    # You can probably also use IE, but maybe not all figure/pages will work - CHECK!
    if ( !is.null(attributes(myBrowser)) ){
      # Set default browser to IE
      myBrowser <- system(\"WHERE /F /R C:\\\\ iexplore.exe\", intern=TRUE)
      # # Search returns multiple hits, take the one in \"Program File\" and \"Internet Explorer\" directory
      myBrowser <- myBrowser[grep(\"Program Files\\\\\",myBrowser)]
      options(browser = myBrowser[1])
    }
    
    message('library paths:\n', paste('... ', .libPaths(), sep='', collapse='\n'))
    
    launch.browser = function(appUrl, browser.path=myBrowser) {
      message('Browser path: ', browser.path)
      #  shell(sprintf('\"%s %s\"', browser.path, appUrl), intern=TRUE)
      system(\"cmd.exe\", input=paste(browser.path, appUrl, sep=\" \"), intern=TRUE)
    }
    
    shiny::runApp('./', launch.browser=launch.browser)
    
    # Perhaps remove temporary directory that holds the packages?
    # unlink(\"",libLoc,"\")
    
  ")
  
  # Write the txt to the file
  writeLines(txt,fileConn)
  
  # Close the connection
  close(fileConn)
  
}

#### Microarray ####
writeUI_Microarray <- function(fileOut="ui.R",tt="tt_all.rds", eset="eset.rds",...){
  # Open the file to write to
  fileConn <- file(fileOut)
  
  # Gather the txt
  txt <- paste0("
library(shiny)
library(ggvis)
library(Biobase)
library(markdown)

# check existence of needed files - R 3.2.0
#if(!dir.exists(\"venn\")){
#  cat(\"No Venn information found, directory \'venn\' does not exist!\\n\")
#}

# check existence of needed files - pre R 3.2.0
if(\"venn\" %in% dir()==FALSE) {
  cat(\"No Venn information found, directory \'venn\' does not exist!\\n\")
}

# read in dataset
if(file.exists(\"",tt,"\")){
  tt <- readRDS(\"",tt,"\")
} else {
  cat(\"",tt," could not be found\\n\")
}

if(file.exists(\"",eset,"\")){
  eset <- readRDS(\"",eset,"\")
} else {
  cat(\"",eset," could not be found\\n\")
}

# define title
title <- \"Microarray Analysis\"
                
# Define UI for microarray application
shinyUI(fluidPage(
# Theme (https://github.com/ClintWeathers/stuff_about_r/blob/master/datatable_links_plus_theme.md)
#  theme = shinytheme(\"cosmo\"), 
                
# Application title
headerPanel(title),
                
mainPanel(
  tabsetPanel(type=\"tabs\",
                
    # Tab for QC
    # See: http://stackoverflow.com/questions/22177974/r-shiny-using-iframe-for-local-files
    tabPanel(\"Quality Control\",
      tabsetPanel(
      # link naar html pagina's van zip-map QC
        ", 
        if (dir.exists("aQM_before_neqc")){
          if (dir.exists("aQM_neqc")){"
            tabPanel(\"before neqc\",
            htmlOutput('beforeQC')
            ),        
            "} else {"
            tabPanel(\"before neqc\",
            htmlOutput('beforeQC')
            )"
          }
        },
        if (dir.exists("aQM_neqc")){"
          tabPanel(\"after neqc\",
          htmlOutput('afterQC')
          )"
        },"
      )
    ),
                
    # Tab for Differential Expression Plots
    tabPanel(\"Differential Expression Plots\",
      fluidRow(
        column(6,   
          selectInput(\"x\", \"Comparison:\",gsub(\")\",\"\",gsub(\"P.Value \\\\(\",\"\",grep(\"P.Value \\\\(\",names(tt),value=TRUE))), width=\"600px\")
        )
      ),
      fluidRow(
        column(6,
          selectInput(\"cDEP\", \"Color by factor:\",names(pData(eset)))
        ),
        column(6,
          sliderInput(\"log2FC\", \"Log2FC cutoff:\", 0.5, 10, 1, step = NULL, round = FALSE, ticks = TRUE, animate = FALSE, width = NULL, sep = \",\", pre = NULL, post = NULL, timeFormat = NULL, timezone = NULL, dragRange = TRUE)
        )
      ),
      fluidRow(
        column(6,
          selectInput(\"fDEP\", \"Group by factor:\",names(pData(eset)))
        ),
        column(6,
          sliderInput(\"pVal\", \"p-value cutoff:\", 0.00001, 0.1, 0.05, step = NULL, round = FALSE, ticks = TRUE, animate = FALSE, width = NULL, sep = \",\", pre = NULL, post = NULL, timeFormat = NULL, timezone = NULL, dragRange = TRUE)
        )                           
      ),
      hr(),
                                                                                                        
      fluidRow(
        column(6,
          textInput(\"g\", \"Gene (Symbol or EntrezID, separate by ','):\", value = \"1\")
        )
      ),
                                                                                                        
      hr(),
                                                                                                        
      fluidRow(
        verbatimTextOutput(\"volcano_sel\"),                           
        column(6,div(style = \"height:500px;width:500px;\",
          ggvisOutput(\"volcano_plot\"))
        ),
        column(6,div(style = \"height:500px;width:500px;\",
          plotOutput(\"beeswarm_plot\",width=400, height=400) )       
        )
        # needed for proper page layout
        #    cellArgs = list(style=\"width: 400px;\")    
       )
    ),
                                                                                                        
    # Tab for Differenatial Expression Table
    tabPanel(\"Differential Expression Table\",
      fluidRow(
        column(6,   
          selectInput(\"tabX\", \"Comparison:\",gsub(\")\",\"\",gsub(\"P.Value \\\\(\",\"\",grep(\"P.Value \\\\(\",names(tt),value=TRUE))), width='600px'),
          downloadButton('downloadData', 'Download')            
        )
      ),
                                                           
      hr(),
                                                                                                                                                                                                   
      fluidRow(
        tags$head(tags$style(HTML( \".has-feedback .form-control { padding-right: 0px;}\" ))),
        column(12,
          DT::dataTableOutput(\"DE_table\")           
        )
      )
    ),
                                                                                                                                                                                                   
    # Tab for Venn diagrams
    tabPanel(\"Venn diagrams\",
      fluidRow(
        column(12,   
          selectInput(\"vennX\", \"Comparison:\",gsub(\".html\",\"\",gsub(\"venn_\",\"\",list.files(path=\"venn\", pattern=\"*all*\"))), width='600px')
        )
      ),
                                                                                                                      
      hr(),
                                                                                                                                                                       
      fluidRow(
        column(12,
          htmlOutput('venn')      
        )
      )                
    ),
                                                                                                                                                                                                   
    # Tab for Gene Set Enrichment Analysis
      tabPanel(\"GeneSet Enrichment\",
        fluidRow(
          column(6,   
            selectInput(\"tabGSEA\", \"Comparison:\",gsub(\")\",\"\",gsub(\"P.Value \\\\(\",\"\",grep(\"P.Value \\\\(\",names(tt),value=TRUE))), width='600px')
          )
        ),
        fluidRow(
          column(6,
            selectInput(\"cGSEA\", \"Color by factor:\",names(pData(eset)))
          )
        ),
        fluidRow(
          column(6,
            selectInput(\"fGSEA\", \"Group by factor:\",names(pData(eset)))
          )
        ),
                                                                                                                                                                                                                                                                                                 
        hr(),
                                                                                                                                                                                                                                                                                                 
        fluidRow(
          tags$head(tags$style(HTML( \".has-feedback .form-control { padding-right: 0px;}\" ))),
          column(12,
            DT::dataTableOutput(\"genSetTab\")
          ),
          hr(),
          column(12,
            DT::dataTableOutput(\"genSetTab_genes\")
          )
        ),
        absolutePanel(
          top = 60, right = 20, width = 300, draggable = TRUE, style = \"opacity: 0.92\",
          wellPanel(
            HTML(markdownToHTML(fragment.only=TRUE, text=c(\"Beeswarm plot of the selected Gene (`draggable`)\"))),
            plotOutput(\"beeswarm_plot2\", height=\"200px\")
          )
        ),
        absolutePanel(
          top = 20, right = 60, width = 300, draggable = TRUE, style = \"opacity: 0.92\",
          wellPanel(
            HTML(markdownToHTML(fragment.only=TRUE, text=c(\"Volcano plot of the selected GeneSet (`draggable`)\"))),
            plotOutput(\"genSetVolcanoPlot\", height=\"200px\")
          )          
        )
      ),

      # Generate Beeswarm and Probe expression barplots for selected Genes
      tabPanel(\"Selected Genes\",
        fluidRow(
          column(6,   
            textInput(\"s\", \"Gene (Symbol or EntrezID, separate by ','):\", value=\"\")
          )
        ),
        fluidRow(
          column(6,
            selectInput(\"cSel\", \"Color by factor:\",names(pData(eset)))
          )
        ),
        fluidRow(
          column(6,
            selectInput(\"fSel\", \"Group by factor:\",names(pData(eset)))
          )
        ),
        fluidRow(
          column(6,   
            textInput(\"dir\", \"Files will be saved into <dir>.zip:\",value=\"\"),
            downloadButton('downloadPictures', 'Download')                          
          )
        )
      )   
    )
  )
                                                                                                                                                                                                                                                                                                 
))
  ")
  
  # Write the txt to the file
  writeLines(txt,fileConn)
  
  # Close the connection
  close(fileConn)
  
}

writeServer_Microarray <- function(fileOut="server.R",tt="tt_all.rds", eset="eset.neqc.rds", BroadSets="BroadSets.rds", tabCameraFile="tab.camera.v4.0_",...){
  # Open the file to write to
  fileConn <- file(fileOut)
  
  # Gather the txt
  txt <- paste0("
# Check for packages
# http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
if (!require(\"pacman\")){ install.packages(\"pacman\") }
pacman::p_load(devtools,curl,shiny,shinythemes,DT,ggvis,beeswarm,plyr,markdown,RColorBrewer,randomcoloR)
if (!require(\"d3heatmap\")){ pacman::p_install_gh(\"rstudio/d3heatmap@colorlabels\") }
                
# Biobase and limma are packages from Bioconductor
if(require(\"limma\") == FALSE) {
  source(\"http://bioconductor.org/biocLite.R\")
  biocLite(\"limma\")
}
if (require(\"Biobase\")  == FALSE) {
  source(\"http://bioconductor.org/biocLite.R\")
  biocLite(\"Biobase\")
}

# Gleaned from
# - http://www.matteodefelice.name/research/2014/12/29/exploring-time-series-data-with-ggplot2-and-shiny/
# - http://web.stanford.edu/~cengel/cgi-bin/anthrospace/building-my-first-shiny-application-with-ggplot
# 
# The latter worked with reactiveText and reactivePlot, both of which are deprecated
#
# For adapting plot sizes etc.:
# - http://shiny.rstudio.com/articles/images.html
  
library(shiny)
library(shinythemes)
library(DT)
library(Biobase)
library(ggvis)
library(beeswarm)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output, session) {
  
  # Get directory
  pwd=getwd()
  
  # Set paths for QC and Venn HTMLs
  addResourcePath(\"pwd\", pwd)
  
  # Initialise
  select <- reactiveValues(id = NULL, symbol = NULL, entrez = NULL)
  
  ### Data import:
  DF <- as.data.frame(do.call(readRDS,list(\"",tt,"\")))
  eset = readRDS(\"",eset,"\")
  BroadSets = readRDS(\"",BroadSets,"\")
  
  # Specify the treatment groups (this is not generic...)
  #    treatment = factor(with(pData(eset),paste(CellType,Treatment,TimePoint,sep=\"_\")))
  # If you have added the treatment factor to pData(eset), you can now extract this info. Hence, more generic
  #    treatment = factor(pData(eset)$Treatment)
  
  # QC
  # We already have the arrayQualityMetrics output with the nicely formatted and linked HTMLs
  # See: http://stackoverflow.com/questions/22177974/r-shiny-using-iframe-for-local-files
  #      http://sikwati.com/blog/make-iframe-fit-container-height/
  output$beforeQC <- renderUI({
    tags$iframe(
      seamless=\"seamless\", align=\"center\", frameborder=\"0\", width=\"100%\", style=\"height:100em\", scrolling=\"yes\",
      src=\"pwd/aQM_before_neqc/index.html\"
    )
  })
  
  output$afterQC <- renderUI({
    tags$iframe(
      seamless=\"seamless\", align=\"center\", frameborder=\"0\", width=\"100%\", style=\"height:100em\", scrolling=\"yes\",
      src=\"pwd/aQM_neqc/index.html\"
    )
  })
  
  # Make the Volcano plots linked with Beeswarm plots for the genes for the different comparisons
  colorDEP <- reactive({
    cFactor <- input$cDEP
    cFactor
  })

  factorDEP <- reactive({
    Factor <- input$fDEP
    Factor
  })

  DF1 <- reactive({
    # check for the input variable
    yCol = paste0(\"P.Value (\",input$x,\")\")
    xCol = input$x
    
    dat = data.frame(x = DF[[input$x]], y = -log10(as.numeric(DF[[yCol]])), symbol = as.character(DF[[\"symbol_reannotated\"]]), entrez = as.character(DF[[\"entrezid_reannotated\"]]))
    dat$id = 1:nrow(dat)
    
    dat
  })
  
  click <- function(data, ...) {
    select$id <- data$id
    select$symbol <- data$symbol
    select$entrez <- data$entrez
    print(data)
    
    treatment <- factor(pData(eset)[,factorDEP()])

    myLevels <- levels(factor(pData(eset)[,colorDEP()]))
    myColorsLevels <- rainbow(length(myLevels))
    myColors <- myColorsLevels[factor(pData(eset)[,colorDEP()])]

    # Beeswarm plot
    output$beeswarm_plot <- renderImage({
      
      # Temporary placeholder
      outfile = tempfile(fileext=\".png\")
      
      # Generate the beeswarm plot
      newTitle = paste0(fData(eset)$symbol_reannotated[select$id],\" (\",DF$entrezid_reannotated[select$id],\")\",\" expression\" )
      png(outfile,width=500, height=430)
      par(mar=c(9.1, 4.1, 2.1, 12.1), xpd=TRUE)
      par(mgp=c(2,1,0))
      beeswarm(exprs(eset)[select$id,]~treatment,data=eset,method=\"swarm\",labels=levels(treatment), vertical=TRUE, pch=15, pwcol=myColors, main=newTitle, xlab=factorDEP(), ylab=\"log2(expression)\")
      legend(\"topright\",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = factorDEP())
      dev.off()
      
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 500,
           height = 430,
           alt = \"This is alternate text\")
    }, deleteFile = TRUE)
  }
  
  output$volcano_sel <- renderText({
    paste0(dim(DF1.logFC.sel())[1],\" (out of \", dim(DF1())[1],\") genes fall within the selected cutoffs (log2FC >= abs(\",input$log2FC,\"), p-value <= \",input$pVal,\")\n\")
  })

  output$ggvis_plot <- renderUI({
    ggvisOutput(\"volcano_plot\")
  })
  
  # http://stackoverflow.com/questions/25011544/how-is-data-passed-from-reactive-shiny-expression-to-ggvis-plot  
  
  #  output[['value1']] <- renderText({ names(DF) })
  
  all_values <- function(x) {
    if(is.null(x)) return(NULL)
    row <- DF1()[DF1()$id == x$id,]
    row <- row[,c(3,4,1,2)]
    names(row) <- c(\"Gene Symbol\",\"EntrezID\",\"Log2FC\",\"-log10(pvalue)\")
    paste0(names(row), \": \", format(row), collapse = \"<br />\")
  }
  
  # Selected Genes
  DF1.sel <- reactive({ 
    if (input$g == 0) return()
    selList = unlist(strsplit(input$g,\",\"))
    data <- unique(DF1()[DF1()$symbol %in% selList | DF1()$symbol_reannotated %in% selList,])
    data <- rbind(data,unique(DF1()[DF1()$entrez %in% selList | DF1()$entrez_reannotated %in% selList,]))
    data <- subset(data,data$symbol!=\"NA\")
    
    print(data)
  })
  
  # Selected Genes
  DF1.logFC.sel <- reactive({ 
    if (input$g == 0) return()
    data <- DF1()[abs(DF1()$x) >= input$log2FC & DF1()$y >= -log10(input$pVal),]
    
    print(data)
  })
  
  # Volcano plot
  reactive({
    DF1 <- DF1()
    DF1.logFC.sel <- DF1.logFC.sel()
    DF1.sel <- DF1.sel()
    DF1 %>%
    ggvis(~x, ~y) %>%
    handle_click(on_click = click) %>%
    layer_points(key := ~id, size := 30 , size.hover := 500, data = DF1) %>%
    layer_points(key := ~id, size := 30 , stroke := \"deepskyblue\", fill := \"deepskyblue\", size.hover := 500, data = DF1.logFC.sel) %>%
    layer_points(key := ~id, size := 30 , stroke := \"red\", fill := \"red\", size.hover := 500, data = DF1.sel) %>%
    add_axis(\"x\", title=\"log2FoldChange\") %>%
    add_axis(\"y\", title=\"-log10(pvalue)\") %>%
    add_tooltip(all_values, \"hover\") %>%
    layer_paths(stroke:='red', data = data.frame(x=range(DF1()$x),y=-log10(input$pVal)) ) %>%
    layer_paths(stroke:='red', data = data.frame(x=-input$log2FC,y=range(DF1()$y)) ) %>%
    layer_paths(stroke:='red', data = data.frame(x=input$log2FC,y=range(DF1()$y)) ) %>%
    set_options(width = 500, height = 400, padding = padding(10, 100, 50, 50)) }) %>%
    bind_shiny(\"volcano_plot\")
  
  default_options()
  
  # Data Table Output
  DF_tab <- reactive({
    # get column index adj.P.Val (last column before the individual comparisons start)
    indx = grep(\"^adj.P.Val$\", colnames(DF))
    indxLogFC = grep(paste0(\"^\",input$tabX,\"$\"),colnames(DF)) # Not really needed!!
    indxStart = grep(paste0(\"t \\\\(\",input$tabX,\")\"),colnames(DF))
    indxEnd = grep(paste0(\"adj.P.Val \\\\(\",input$tabX,\")\"),colnames(DF))
    vec <- c(1:indx,indxStart:indxEnd)
    
    # get URLs for the EntrezIDs
    indxEntrez = grep(\"^entrez\",colnames(DF))
    urls <- c(paste0(\"http://www.ncbi.nlm.nih.gov/gene/\",as.character(DF[,indxEntrez])))
    refs <- paste0(\"<a href='\",  urls, \"'>\",DF[,indxEntrez],\"</a>\")
    DF[,indxEntrez] <- refs
    
    dat = DF[,vec]
    # sort on P.value
    indx = grep(paste0(\"P.Value \\\\(\",input$tabX,\")\"), colnames(dat))
    dat = dat[order(dat[,indx]),]
    dat
  })
  
  indxEntrez = grep(\"^entrez\",colnames(DF))
  output$DE_table <- DT::renderDataTable({ datatable(DF_tab(),rownames = FALSE, escape=indxEntrez-1, options = list(iDisplayLength = 25, search=list(regex = TRUE)),filter=\'top\') })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$tabX, '.csv', sep='') },
    content = function(file) {
      dat = DF_tab()
      indxEntrez = grep(\"^entrez\",colnames(dat))
      noUrls <- gsub(\"<.*'>\",\"\",as.character(dat[,indxEntrez]))
      noUrls <- gsub(\"</a>\",\"\",noUrls)
      dat[,indxEntrez] <- noUrls
      
      write.csv(dat, file)
    }
  )
  
  # Table for Gene Pathway Enrichment per comparison
  colorGSEA <- reactive({
    cFactor <- input$cGSEA
    cFactor
  })
  factorGSEA <- reactive({
    Factor <- input$fGSEA
    Factor
  })
  
  DF2 <- reactive({
    # check for the input variable
    yCol = paste0(\"P.Value (\",input$tabGSEA,\")\")
    xCol = input$tabGSEA
    
    dat = data.frame(x = DF[[xCol]], y = -log10(as.numeric(DF[[yCol]])), symbol = as.character(DF[[\"symbol_reannotated\"]]), entrez = as.character(DF[[\"entrezid_reannotated\"]]))
    dat$id = 1:nrow(dat)
    
    dat
  })
  
  click2 <- function(data, ...) {
    select$id <- data$id
    select$symbol <- data$symbol
    select$entrez <- data$entrez
    print(data)
    
    treatment <- factor(pData(eset)[,factorGSEA()])

    myLevels <- levels(factor(pData(eset)[,colorGSEA()]))
    myColorsLevels <- rainbow(length(myLevels))
    myColors <- myColorsLevels[factor(pData(eset)[,colorGSEA()])]

    # Beeswarm plot
    output$beeswarm_plot2 <- renderImage({
      
      # Temporary placeholder
      outfile = tempfile(fileext=\".png\")
      
      # Generate the beeswarm plot
      newTitle = paste0(fData(eset)$symbol_reannotated[select$id],\" (\",DF$entrezid_reannotated[select$id],\")\",\" expression\" )
      png(outfile,width=500, height=430, pointsize=12)
      par(mar=c(9.1, 4.1, 2.1, 12.1), xpd=TRUE)
      par(mgp=c(2,1,0))
      beeswarm(exprs(eset)[select$id,]~treatment,data=eset,method=\"swarm\",labels=levels(treatment), cex=c(1.5), vertical=TRUE, pch=15,  pwcol=myColors, main=\"\", xlab=\"\", ylab=\"\", yaxt=\"n\")
      legend(\"topright\",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = factorGSEA())
      mtext(as.expression(substitute(bold(newTitle))), side=3, cex = 2.0)
      mtext(factorGSEA(), side=1, line=2.2, cex=2)
      mtext(\"log2(expression)\", side=2, line=2.6, cex=2)
      axis(2,cex.axis=1.5)
      dev.off()
      
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 250,
           height = 215,
           alt = \"This is alternate text\")
    }, deleteFile = TRUE)
  }
  
  output$ggvis_plot <- renderUI({
    ggvisOutput(\"genSetVolcanoPlot\")
  })
  
  all_values2 <- function(x) {
    if(is.null(x)) return(NULL)
    row <- DF2()[DF2()$id == x$id,]
    row <- row[,c(3,4,1,2)]
    names(row) <- c(\"Gene Symbol\",\"EntrezID\",\"Log2FC\",\"-log10(pvalue)\")
    paste0(names(row), \": \", format(row), collapse = \"<br />\")
  }
  
  # Need data frames for ablines
  h_abline <- reactive({ data.frame(x=range(DF2()$x),y=-log10(0.05)) })
  v_abline1 <- reactive({ data.frame(x=-log2(2),y=range(DF2()$y)) })
  v_abline2 <- reactive({ data.frame(x=log2(2),y=range(DF2()$y)) })
  
  # Selected Genes
  #  DF2.sel <- reactive({ 
  #    if(input$g != \"\") {
  #      selList = unlist(strsplit(input$g,\",\"))
  #      data <- unique(DF2()[DF2()$symbol %in% selList,])
  #      data <- rbind(data,unique(DF2()[DF2()$entrez %in% selList,]))
  #      data <- subset(data,data$symbol!=\"NA\")
  #    
  #      print(data)
  #    }
  #  })
  
  # Data Table Output
  genSetTab <- reactive({
    xCol = input$tabGSEA
    
    tab.camera = read.table(file=paste0(\"",tabCameraFile,"\",xCol,\".txt\"),header=TRUE, sep=\"\t\")
      
    # get URLs for the GeneSets@MSigDB and open in new tab
    urls <- c(paste0(\"http://www.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=\",as.character(gsub(\"^C[1-8]_|^H_\",\"\",as.character(tab.camera[,1])))))
    refs <- paste0(\"<a href='\",  urls, \"' target=_blank>\",tab.camera[,1],\"</a>\")
    
    # make the NGenes column clickable to show 
    # - table of genes in the geneset, 
    # - color these genes in the volcano-plot and 
    # - when clicked in the volcano plot display the beeswarm plot
    #    linkGeneSet <- c(paste(sets[[as.character(tab.camera[,1])]],sep=\":\"))
    
    if (nrow(tab.camera)>0){
      #    rownames(tab.camera) <- as.character(tab.camera[,1])
      tab.camera[,1] <- refs
      #    tab.camera[,2] <- linkNGenes
      tab.camera
    } else {
      shiny::showNotification(\"No genesets for this category\", type = \"error\")
      NULL
    }

  })
  
  # Again get the indices from the Broadset coupled with the eset?? And read those in?
  # It would mean we provide the Broadsets.rds as well (which might make sense for provenance...)
  output$genSetTab <- DT::renderDataTable(genSetTab(), rownames=FALSE, escape=-c(1), selection = 'single', options = list(iDisplayLength = 10, search=list(regex = TRUE)), filter=\'top\')
  
  click3 <- function(data, ...) {
    # Beeswarm plot
    output$beeswarm_plot2 <- renderImage({
      # Temporary placeholder
      outfile = tempfile(fileext=\".png\")
      
      png(outfile,width=500, height=430, pointsize=12)
      dev.off()
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 250,
           height = 215,
           alt = \"This is alternate text\")
    }, deleteFile = TRUE )
  }
  
  DF2.sel <- reactive({
    s = input$genSetTab_row_last_clicked
    if (length(s)) {
      s <- strsplit(strsplit(as.character(genSetTab()[s,1]),\">\")[[1]][2],\"<\")[[1]][1]
      cat('These rows were selected:\\n\t')
      cat(s, sep = '\\n')
    }
    # This will couple the Entrez IDs to the BroadSets defined earlier
    indx = unlist(BroadSets[s])
    data <- unique(DF2()[DF2()$entrez %in% indx,])
    data <- subset(data,data$entrez!=\"NA\")
    click3()
    
    data
  })
  
  DF2 %>%
    ggvis(~x, ~y) %>%
    handle_click(on_click = click2) %>%
    layer_points(key := ~id, size := 30 , size.hover := 500, data = DF2) %>%
    layer_points(key := ~id, size := 30 , stroke := \"red\", fill := \"red\", size.hover := 500, data = DF2.sel) %>%
    add_axis(\"x\", title=\"log2FoldChange\") %>%
    add_axis(\"y\", title=\"-log10(pvalue)\") %>%
    add_tooltip(all_values2, \"hover\") %>%
    layer_paths(stroke:='red', data=h_abline) %>%
    layer_paths(stroke:='red', data=v_abline1) %>%
    layer_paths(stroke:='red', data=v_abline2) %>%
    set_options(width = 800, height = 400, padding = padding(10, 50, 50, 50)) %>%
    bind_shiny(\"genSetVolcanoPlot\")
  
  default_options()
  
  DF2a.sel <- reactive({
    s = input$genSetTab_row_last_clicked
    if (length(s)) {
      s <- strsplit(strsplit(as.character(genSetTab()[s,1]),\">\")[[1]][2],\"<\")[[1]][1]
      cat('These rows were selected:\\n\t')
      cat(s, sep = '\\n')
    }
    # This will couple the Entrez IDs to the BroadSets defined earlier
    indx = unlist(BroadSets[s])
    data <- unique(DF2()[DF2()$entrez %in% indx,])
    data <- subset(data,data$entrez!=\"NA\")
    data <- data[,c(3,4,1,2)]
    rownames(data) <- data$id
    data$id <- NULL
    
    # Generate URL for EntrezID
    if (length(data[,2])) {
      urls <- c(paste0(\"http://www.ncbi.nlm.nih.gov/gene/\",as.character(data[,2])))
      refs <- paste0(\"<a href='\",  urls, \"'>\",data[,2],\"</a>\")
      data[,\"entrez\"] <- refs
    }
    
    names(data)[3] <- \"LogFC\"
    names(data)[4] <- \"-log10(pVal)\"
    data
  })
  
  
  output$genSetTab_genes <- DT::renderDataTable(DF2a.sel(), rownames=FALSE, escape=FALSE, options = list(iDisplayLength = 10, search=list(regex = TRUE)), filter=\'top\')
  # Replace initial VolcanoPlot with a MA-plot?
  
  # Venn diagrams
  # We already have the Venn diagrams with the nicely formatted and linked HTMLs
  # See: http://stackoverflow.com/questions/22177974/r-shiny-using-iframe-for-local-files
  DF_Venn <- reactive({
    # check for the input variable
    comp = paste0(\"venn_\",input$vennX,\".html\")
    
    comp
  })
  
  getPage<-function(x) {
    return(
      tags$iframe(
        src = paste0(\"pwd/venn/\",x), seamless=\"seamless\", width=\"100%\", style=\"height:50em\",  frameborder=\"0\", scrolling=\"yes\"
      )
    )
  }
  
  output$venn <- renderUI({
    #    tags$iframe(
    #      seamless=\"seamless\", align=\"center\", frameborder=\"0\", style=\"width:100%\", height=\"100%\", scrolling=\"yes\",
    #      src=paste0(\"pwd/venn/\",DF_Venn())
    #    )
    x=paste0(\"venn_\",input$vennX,\".html\")
    getPage(x)
  })
  
  # Generate Beeswarm, Probe Expression Barplots for selected genes and save to zip file
  # See: http://stackoverflow.com/questions/28228892/download-multiple-csv-files-in-a-zipped-folder-in-shiny
  colorSel <- reactive({
    cFactor <- input$cSel
    cFactor
  })
  factorSel <- reactive({
    Factor <- input$fSel
    Factor
  })

  DF4.sel <- reactive({ 
    if(input$s != \"\") {
      selList = unlist(strsplit(input$s,\",\"))
      print(selList)
      indxEntrez = grep(\"^entrez\",colnames(DF))
      indxSYMBOL = grep(\"^SYMBOL|^Symbol|^symbol\",colnames(DF))
      data <- unique(DF[DF[,indxSYMBOL] %in% selList,])
      data <- rbind(data,unique(DF[DF[,indxEntrez] %in% selList,]))
      #      data <- subset(data,data[,indxSYMBOL]!=NA)
      
    }
  })
  
  # Get a filename for the directory that will hold the pictures (zipped dir)
  zipDir <- reactive({
    if (input$dir != \"\"){
      tmp = paste0(input$dir,\".zip\")
    }
    tmp
  })
  
  output$downloadPictures <- downloadHandler(
    filename = zipDir,
    content = function(fname) {
      tmpdir <- tempdir()
      setwd(tempdir())
      
      treatment <- factor(pData(eset)[,factorSel()])

      # We can order the labels based on the \"treatment\" factor: just assume that this has been constructed based on the
      # importance of classification
      #indx.df = as.data.frame(strsplit(as.character(pData(eset)$Treatment),\"_\"))
      #indx.lab=order(pData(eset)$CellType,pData(eset)$Treatment,pData(eset)$TimePoint)
      # or just by the last chosen Factor
      indx.df = as.data.frame(strsplit(as.character(treatment),\"_\"))
      nrFactors = length(indx.df[[1]])
      indx.df.t = as.data.frame(t(indx.df))
      indx.df.t$id = 1:nrow(indx.df.t)
      for (i in 1:nrFactors) {
        indx.lab = indx.df.t[order(indx.df.t[,i]),]
      }
      indx.lab = indx.lab$id

      myLevels <- levels(factor(pData(eset)[,colorSel()]))
      myColorsLevels <- rainbow(length(myLevels))
      myColors <- myColorsLevels[factor(pData(eset)[,colorSel()])]

      indxSYMBOL = grep(\"^SYMBOL|^Symbol|^symbol\",colnames(DF4.sel()))
      fs = c()
      for (gene in DF4.sel()[,indxSYMBOL[1]]){
        indx = featureNames(eset)[fData(eset)[,indxSYMBOL[1]] %in% gene]
        if (length(indx) >= 1){
          for (i in 1:length(indx)){
            newTitle = paste0(gene,\" (\",indx[i],\")\",\" expression\" )
            png(paste0(\"filteredProbes\",gene,\"_\",indx[i],\".png\") )
            par(mar=c(9.1, 4.1, 4.1, 17.1), xpd=TRUE)
            beeswarm(exprs(eset)[indx[i],]~treatment,data=eset,method=\"swarm\",labels=levels(treatment), vertical=TRUE, pch=15,  pwcol=myColors, main=newTitle, xlab=factorSel(), ylab=\"log2(expression)\")
            legend(\"topright\", inset=c(-1.4,0),legend = myLevels, pch = 15, col = myColorsLevels, title = \"Sample\")
            dev.off()
            jpeg(paste0(\"filteredProbes\",gene,\"_\",indx[i],\"_barplot.jpg\") )
            par(mar=c(15.1, 4.1, 4.1, 4.1), xpd=TRUE)
            barplot(exprs(eset)[indx[i],indx.lab],names.arg=with(pData(eset),rownames(pData(eset))), main=newTitle, las=2, cex.names=0.8, ylab=\"log2(expression)\")
            dev.off()
            
            fs = c(fs,paste0(\"filteredProbes\",gene,\"_\",indx[i],\".png\"), paste0(\"filteredProbes\",gene,\"_\",indx[i],\"_barplot.jpg\"))
          }
        }
      }
      
      toCSV.dat <- data.frame(GeneName=rownames(exprs(eset)[grep(paste(DF4.sel()[,indxSYMBOL[1]],collapse=\"|\"), fData(eset)[,indxSYMBOL[1]]),]),
                              exprs(eset)[grep(paste(DF4.sel()[,indxSYMBOL[1]],collapse=\"|\"), fData(eset)[,indxSYMBOL[1]]),])
      
      write.table(toCSV.dat,file=\"selectedGenes.txt\", row.names=FALSE, quote=FALSE, sep=\"\t\")
      
      fs = c(fs,\"selectedGenes.txt\")

     if (length(fs) > 0){
        zip(zipfile=fname, files=fs)
        
        if(file.exists(paste0(fname, \".zip\"))) {
          file.rename(paste0(fname, \".zip\"), fname)
        }
      }
      
    },
    contentType = \"application/zip\"
  )
  
    
})
  ")
  
  # Write the txt to the file
  writeLines(txt,fileConn)

  # Close the connection
  close(fileConn)

}

#### RNASeq ####
writeUI_RNASeq <- function(fileOut="ui.R",tt="tt_all.rds", v="voom.rds", object="voom", 
                           venn = "venn", venn_gs = "venn_gs", expinfo=NULL, ...){
  # Open the file to write to
  fileConn <- file(fileOut)
  
  # Read the specs for edgeR & voom objects
  if ( object == "edgeR"){
    PVAL = "PValue"
    TARGETS = "names(v$targets)"
  } else {
    if ( object == "voom"){
      PVAL = "P.Value"
#      TARGETS = "names(v$targets)[-c(1:3)]"
      TARGETS = "factorNames"
    }
  }
  
  # Is there any experimental info present?
  if ( is.null(expinfo) ){
    expInfo =""
  } else {
    expInfo = paste0("# Tab for Experimental Info
         tabPanel(\"Exp. Info\",
                  fluidRow(
                    column(12,
                      pre(includeText(\"",expinfo,"\"))     
                    )
                  )        
        ),") 
  }
    
  # Gather the txt
  txt <- paste0("
# Check for packages and suppress updating packages
# http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
if (!requireNamespace(\"BiocManager\") && tools:::.BioC_version_associated_with_R_version() >= 3.8 ){
    install.packages(\"BiocManager\")
}

if (!require(\"pacman\")){ 
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"pacman\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"pacman\", suppressUpdates=TRUE)
  }
}

# install packages using \"pacman\"
pacman::p_load(devtools,curl,shiny,shinythemes,DT,ggvis,ggplot2,beeswarm,plyr,markdown,RColorBrewer,
               randomcoloR,colourpicker,cwhmisc,zip,gplots,
               update = FALSE
       )

# Biobase, limma and edgeR are packages from Bioconductor
if (require(\"Biobase\") == FALSE){ 
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"Biobase\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"Biobase\", suppressUpdates=TRUE)
  }
}
if(require(\"edgeR\") == FALSE) {
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"edgeR\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"edgeR\", suppressUpdates=TRUE)
  }
}
if (require(\"limma\")==  FALSE) {
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"limma\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"limma\", suppressUpdates=TRUE)
  }
}
if (require(\"ComplexHeatmap\") == FALSE) {
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"ComplexHeatmap\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"ComplexHeatmap\", suppressUpdates=TRUE)
  }
}

library(shiny)
library(shinythemes)
library(DT)
library(Biobase)
library(ggvis)
library(ggplot2)
library(beeswarm)
library(plyr)
library(limma)
library(edgeR)
library(RColorBrewer)
library(randomcoloR)
library(colourpicker)
library(zip)

# Extra function for the slider
# https://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale

# logifySlider javascript function
JS.logify <-
\"
// function to logify a sliderInput
function logifySlider (sliderId, sci = false) {
  if (sci) {
    // scientific style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
    })
  } else {
    // regular number style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return (Math.pow(10, num)); }
    })
  }
}\"
                
# call logifySlider for each relevant sliderInput
JS.onload <-
\"
// execute upon document loading
$(document).ready(function() {
  // wait a few ms to allow other scripts to execute
  setTimeout(function() {
    // include call for each slider
    logifySlider('log_slider', sci = false)
    logifySlider('pVal', sci = true)
  }, 5)})
\"

# check existence of needed files - R 3.2.0
#if(!dir.exists(\"",venn,"\")){
#  cat(\"No Venn information found, directory \'venn\' does not exist!\\n\")
#}
                
# check existence of needed files - pre R 3.2.0
if(dir.exists(\"", venn,"\")  == FALSE) {
  cat(\"No Venn information found, directory \'",venn, "\' does not exist!\\n\")
}

# See whether there is information about the experiment
if(file.exists(\"",expinfo,"\")){
  expInfo <- 1
} else {
  cat(\"",expinfo," could not be found\\n\")
}
                
# read in tt
if(file.exists(\"",tt,"\")){
  tt <- readRDS(\"",tt,"\")
} else {
  cat(\"",tt," could not be found\\n\")
}

# read in v
if(file.exists(\"",v,"\")){
  v <- readRDS(\"",v,"\")
  indx <- grep(c(\"group|lib.size|norm.factors\"), colnames(v$targets))
  if (length(indx) > 0){
    nrLevels <- apply(subset(v$targets, select=-indx), 2, function(x) length(levels(factor(x))))
    factorNames <- names(subset(v$targets, select=-indx))
  } else {
    nrLevels <- apply(v$targets, 2, function(x) length(levels(factor(x))))
    factorNames <- names(v$targets)
  }
  nrLevels <- nrLevels[which(nrLevels/nrow(v$targets) <= 0.6)]

  legendWidth <- sum(rev(sort(apply(v$targets, 2, function(x) max(nchar(x)))))[c(1:2)])
  legendHeight <- nrow(v$targets)

  # Get the contrasts and design (needed for the Venn diagram function)
  my.contrasts <- v$contrasts
  design <- v$design
  # Get the list of indices of the contrasts for the Venn diagrams - not always present!!
  # if not present, try to get the information from the original venn directory?
  if (!is.null(v$vennList)){
    vennList <<- v$vennList
  } else {
    presentVenns <- c(dir(\"", venn,"\",pattern=\"*.png\"))
    presentVenns <- gsub(\"^venn_\",\"\",gsub(\"\\\\.png\",\"\",presentVenns))
    # Extract the names, count the \"_vs_\"s and split on 2nd, 4th etc.
    indx <- as.data.frame(strsplit2(presentVenns,\"_vs_|_and_\"))
    # Make combinations
    indx1 <- apply(indx,1,function(x) paste(x[seq(1,length(x),1)]))
    if ( dim(indx)[2] >= 2) {
      indx2 <- apply(indx,1,function(x) paste(x[seq(1,length(x),2)],x[seq(2,length(x),2)], sep=\"_vs_\"))
      indxAll <- rbind(indx1,indx2)
      if ( dim(indx)[2] >= 3) {
         indx3 <- apply(indx,1,function(x) paste(x[seq(1,length(x),2)],x[seq(2,length(x),2)], x[seq(3,length(x),3)],sep=\"_vs_\"))
         indxAll <- rbind(indxAll,indx3)
         if ( dim(indx)[2] >= 4) {
            indx4 <- apply(indx,1,function(x) paste(x[seq(1,length(x),2)],x[seq(2,length(x),2)], x[seq(3,length(x),3)], x[seq(4,length(x),4)],sep=\"_vs_\"))
            indxAll <- rbind(indxAll,indx4)
         }
      }
    } else {
      indxAll <- indx1
    }
    # Now match with colnames of my.contrasts to get the indices
    indx <- as.data.frame(apply(as.data.frame(indxAll),2,function(x) match(x,colnames(my.contrasts))))
    vennList <<- list()
    for (i in 1:dim(indx)[2]){
      indices <- c(indx[,i])
      indices <- indices[!is.na(indices)]
      vennList[[i]] <<- indices
    }
    rm(indices,indx,presentVenns)
  }
} else {
  cat(\"",v," could not be found\\n\")
}

# define title
title <- \"RNASeq Analysis\"

# Define UI for RNASeq application
shinyUI(fluidPage(
  # Theme (https://github.com/ClintWeathers/stuff_about_r/blob/master/datatable_links_plus_theme.md)
  theme = shinytheme(\"cosmo\"), 
  
  # Set the color of the Warning/Error messages
  tags$head(
    tags$style(HTML(\"
      .shiny-output-error-validation { color: blue; }
    \"))
  ),

  # Application title
  headerPanel(title),
  
  mainPanel(
    tabsetPanel(type=\"tabs\",
                # Tab for Introduction of app
                tabPanel(\"Introduction\",
                    fluidRow(column(
                      h2(\"Welcome to the RNASeq Analysis Shiny App!\"),
                      br(),
                      tags$p(HTML(
                       \"The RNASeq Analysis Shiny App takes the R objects resulting from an analysis performed on your RNASeq data and 
                        allows interactive visualization and inspection of the obtained results.\nThe results are organized using different tabs:\"
                      )),
                      tags$p(HTML(
                        '<ol start=\\\'1\\\'> 
                         <li><font color=\"blue\">\"Exp. Info\"</font> shows the experimental information regarding your data, the analysis steps performed, some observations regarding the data and possibly some conclusions.</li>
                         <li><font color=\"blue\">\"FastQC of final BAM files\"</font> <i><b>(if present)</b></i> displays the results of the FastQC analysis of the BAM files obtained after alignment. Each colored tile is clickable and will show you the results of the analysis for that sample for that particular statistic. By clicking on the links in the rows starting with \"Collated ..\", you can obtain the information on that statistic for all samples. This way you can easily find possible outliers on sample level.</li>
                         <li><font color=\"blue\">\"After normalization & filtering\"</font> shows you a multi-dimensional scaling (MDS) plot in which you can color samples based on 2 different factors to identify groups of samples. You can also click on samples to get their name. A static heatmap (\"Heatmap\" tab) is shown and a somewhat more dynamic one (\"iHeatmap\" tab) in which you can specify which factors you want to show..</li>
                         <li><font color=\"blue\">\"Differential Expression Plots\"</font> gives you the results of the differentially expressed genes (DEG) analysis. Via a drop-down menu you can select which comparison you want to investigate. The corresponding \"volcano plot\" is then loaded and you can click on the dots to obtain a plot of the expression values of that particular gene in the various samples. The \"Color by Factor\" and \"Group by Factor\" menus allow you to change the appearance of this \"beeswarm plot\" accordingly. You can also display your pet gene by selecting it by name using different types of identiers. Finally you can save the plots in different formats and resolutions for use in presentations.\n<b><font color=\"red\">The generated figures are by no means publication ready and should therefore definitely not be used as such !!</font></b></li>
                         <li><font color=\"blue\">\"Differential Expression Table\"</font> shows you the data underlying all plots. Again, you can specify the comparison you are interested in. You can sort on column information by clicking the small arrows displayed next to the column header. You can filter by clicking in the boxes just below the column headers. This will give you the filter options. You can also click on gene identifiers to obtain more information at the NCBI or Ensembl webpages.</li>
                         <li><font color=\"blue\">\"Venn diagrams - Genes\"</font> allows you to see the intersections and differences of genes significantly expressed in different comparisons. The initial plot has been pre-calculated, but you can change this to your liking using different colors and/or cut-offs. The numbers in the plot are clickable and will show the corresponding genes. Genes have to show the same sign of expression up/down regulation in order to be counted as overlapping. For a statistically sound analysis it is advised to use the Benjamin-Hochberg multiple testing correction and to not use a cut-off on the log2 fold change  (that is, stick to the default value of 0).</li>
                         <li><font color=\"blue\">\"Geneset Enrichment\"</font> shows the result of a CAMERA <a href=\"https://academic.oup.com/nar/article/40/17/e133/2411151\"> (Wu & Smyth, 2012)</a> analysis using genesets obtained from the <a href=\"http://software.broadinstitute.org/gsea/msigdb/index.jsp\">MSigDB</a>. You can again select the comparison of interest. By clicking on the name of the geneset, you will be led to the corresponding page at the MSigDB site. Here you can obtain much more information on the geneset. Just clicking the row of the geneset will highlight the genes in that particular geneset in the volcanoplot (in red). Information on these genes will be displayed below the geneset table. You can then again click on these genes to get the corresponding \"beeswarm plot\". Coloring and grouping can be performed as well. Selection of specific \"Geneset Categories\" is also an option. One can also provide the EntrezGene ID of gene(s) and get only the genesets these occur in.</li> 
                         <li><font color=\"blue\">\"Venn diagrams - Gene Sets\"</font> allows you to see the intersections and differences of genesets significantly enriched in different comparisons. The initial plot has been pre-calculated, but you can change this to your liking using different colours and/or cut-offs. The numbers in the plot are clickable and will show the corresponding gene sets. Gene sets have to show the same sign of up/down regulation in order to be counted as overlapping.</li>
                         <li><font color=\"blue\">\"Selected Genes\"</font> finally will allow you to obtain some information on a gene or set of genes you are interested in. It currently consists of a barplot and beeswarm plot of the expression.</li>
                         </ol>'
                      )),
                      br(),
                      tags$h4(HTML('<font color=\"red\">This Shiny app is always under construction...!!</font>')),
                      br(),
                      tags$p(HTML(
                        'For questions, bug reporting, feature requests, remarks, comments, installation problems etc. please do not hesitate to
                                  <a href=\"mailto:a.jongejan@amc.uva.nl\">contact me!</a>'
                      )),
                      width = 12)
                    )
                ),

                ",expInfo,"

                # Tab for Fastqc files
                tabPanel(\"FastQC of final BAM files\",
                    fluidRow(
                      column(12,
                        htmlOutput('fastqc')      
                      )
                    )                
                ),
                # Tab for MDS and heatmap
                tabPanel(\"After normalization\\n& filtering\",
                         tabsetPanel(
                           tabPanel(\"MDS plot\",
                                    fluidRow(
                                      column(6,
                                        selectInput(\"f1MDSplot\", \"Color by factor 1:\",",TARGETS,")
                                      )
                                    ),
                                    fluidRow(
                                      column(6,
                                        selectInput(\"f2MDSplot\", \"Color by factor 2:\",",TARGETS,")
                                      )
                                    ),
                                    hr(),
                                    fluidRow(
                                      column(width = 7,
                                        plotOutput('mdsPlot', click=\"mdsPlot_click\", width=550, height=550)
                                      ),    
                                      column(width = 5, offset=0.3,
                                        verbatimTextOutput(\"mdsPlot_clickinfo\")
                                      )
                                    ),
                                    absolutePanel(
                                      bottom = 120, right = 20, width = paste0(legendWidth*9,\"px\"), height = paste0(legendHeight*18,\"px\"), draggable = TRUE,
                                      wellPanel(
                                        HTML(markdownToHTML(fragment.only=TRUE, text=c(  \"Legend (`draggable`)\" ))),
                                        plotOutput(\"legendMDSplot\", height=paste0(legendHeight*18,\"px\"))
                                      ),
                                      style = \"opacity: 0.92\"
                                    )
                           ),
                           tabPanel(\"Heatmap\",
                                    plotOutput('aheatmapPlot')
                           ),
                           tabPanel(\"iHeatmap\",
                                    fluidRow(
                                      checkboxGroupInput(\"fD3Heatmap\", \"Color by factor:\",names(nrLevels), selected = names(nrLevels)[1],inline = TRUE),
                                      checkboxInput(\"hD3Heatmap\",label=\"Hide factors\", value=FALSE, width=NULL)
                                    ),
                                    hr(),
                                    #  uiOutput(\"ui_heatmap\")
                                    plotOutput(\"heatmap\",height = 800),
                                    h3(\"Save pictures:\"),
                                    fluidRow(
                                      column(
                                        width = 12,
                                        column(width = 2, selectInput(\"fres_heatmap\", \"Res\", choices=c(\"100\",\"200\",\"300\"), selected = \"100\")),
                                        column(width = 2, selectInput(\"fformat_heatmap\", \"File Type\", choices=c(\"png\",\"tiff\",\"jpeg\",\"pdf\"), 
                                                                      selected = \"png\", multiple = FALSE, selectize = TRUE))
                                      ),
                                      tags$head(tags$style(\".butt{background-color:#337ab7;} .butt{color: #fff;}\")),
                                      column( 3, offset=3,
                                              downloadButton('ch_download', 'Download Plots', class = \"butt\")
                                      )
                                    ),
                                    hr()
                              )
                         )
                ),          
                
                # Tab for Differential Expression Plots
                tabPanel(\"Differential Expression Plots\",
                         fluidRow(
                           column(6,   
                              selectInput(\"x\", \"Comparison:\",gsub(\")\",\"\",gsub(\"",PVAL," \\\\(\",\"\",grep(\"",PVAL," \\\\(\",names(tt),value=TRUE))), width=\"600px\")
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput(\"cDEP\", \"Color by factor:\",",TARGETS,")
                           ),
                           column(6,
                             sliderInput(\"log2FC\", \"Log2FC cutoff:\", 0.5, 10, 1, step = NULL, round = FALSE, ticks = TRUE, animate = FALSE, width = NULL, sep = \",\", pre = NULL, post = NULL, timeFormat = NULL, timezone = NULL, dragRange = TRUE)
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput(\"fDEP\", \"Group by factor:\",",TARGETS,")
                           ),
                           # https://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale
                           tags$head(tags$script(HTML(JS.logify))),
                           tags$head(tags$script(HTML(JS.onload))),
                           tags$head(tags$script(HTML('Shiny.addCustomMessageHandler(\"jsCode\", function(message) { eval(message.value); });'))),
                           column(6,
#                             sliderInput(\"pVal\", \"p-value cutoff:\", 0.0000000001, 0.1, 0.05, step = NULL, round = FALSE, ticks = TRUE, animate = FALSE, width = NULL, sep = ",", pre = NULL, post = NULL, timeFormat = NULL, timezone = NULL, dragRange = TRUE)

                           sliderInput(\"pVal\", \"pVal cutoff (log10):\", min = -10, max = 0, value = -0.5, step = 0.5)
                           )
                         ),
                         
                         hr(),
                         
                         fluidRow(
                           column(6,
                                  textInput(\"g\", \"Gene (SYMBOL, EntrezID or Ensembl Gene ID, separate by ','):\", value = \"1\")
                           )
                         ),
                         
                         hr(),
                         
                         fluidRow(
                           verbatimTextOutput(\"volcano_sel\"),                           
                           column(7,div(style = \"height:500px;width:550px;\",
                                        ggvisOutput(\"volcano_plot\"))
                           ),
                           column(5,div(style = \"height:500px;width:500px;\",
                                        plotOutput(\"beeswarm_plot\",width=400, height=400) )       
                           )
                           # needed for proper page layout
                           #    cellArgs = list(style=\"width: 400px;\")    
                         ),
                         h3(\"Save pictures:\"),
                         fluidRow(
                           column(
                             width = 12,
                             column(width = 2, numericInput(\"fheight_volcano\", \"Height (cm)\", min=2, max=15, step=1, value = 10)),
                             column(width = 2, numericInput(\"fwidth_volcano\", \"Width (cm)\", min=2, max=15, step=1, value = 10)),
                             column(width = 2, selectInput(\"fres_volcano\", \"Res\", choices=c(\"100\",\"200\",\"300\"), selected = \"100\")),
                             column(width = 2, selectInput(\"fformat_volcano\", \"File Type\", choices=c(\"png\",\"tiff\",\"jpeg\",\"pdf\"), 
                                                    selected = \"png\", multiple = FALSE, selectize = TRUE))
                           ),
                           tags$head(tags$style(\".butt{background-color:#337ab7;} .butt{color: #fff;}\")),
                           column( 3, offset=3,
                             downloadButton('bn_download', 'Download Plots', class = \"butt\")
                           )
                           # needed for proper page layout
                           #    cellArgs = list(style=\"width: 400px;\")    
                         ),
                         hr()
  
                ),
                
                # Tab for Differential Expression Table
                tabPanel(\"Differential Expression Table\",
                         fluidRow(
                           column(6,   
                                  selectInput(\"tabX\", \"Comparison:\",gsub(\")\",\"\",gsub(\"",PVAL," \\\\(\",\"\",grep(\"",PVAL," \\\\(\",names(tt),value=TRUE))), width='600px'),
                                  downloadButton('downloadData', 'Download')            
                           )
                         ),
                         
                         hr(),
                         
                         fluidRow(
                           # see https://stackoverflow.com/questions/34863466/data-table-column-filters-not-displaying-properly-in-shiny-using-dt-package
                           # The server=FALSE does not work ...
                           tags$head(tags$style(HTML( \".has-feedback .form-control { padding-right: 0px;}\" ))),
                           column(12,
                                  DT::dataTableOutput(\"DE_table\")
                                  
                           )
                         )
                ),
                
                # Tab for Venn diagrams - Genes
                tabPanel(\"Venn diagrams - Genes\",
                         fluidRow(
                           column(2, numericInput(\"pVal_venn\", \"Cut off p-value:\",0.05)),
                           column(2, numericInput(\"lfc_venn\", \"log Fold-change:\",value=0.0, min=0.0)),
                           column(3, selectInput(\"adjustM_venn\", \"Multiple testing correction\",c(\"none\",\"Benjamini-Hochberg\")))
                         ),
                         fluidRow(
                           column(
                             width = 7,
                             # https://stackoverflow.com/questions/44784083/shiny-label-position-textinput
                             tags$form(
                               class=\"form-horizontal\",
                               tags$div(
                                 class=\"form-group\",
                                 tags$label(class = \"col-sm-4 control-label\", `for` = \"colUp\", br(), \"Color for genes:    Up\"),
                                 column(width = 2, colourInput(\"colUp\", \"\", \"red\", palette=\"limited\", showColour=\"background\")),
                                 tags$label(class = \"col-sm-2 control-label\", `for` = \"colDown\", br(), \"Down\"),
                                 column(width = 2, colourInput(\"colDown\", \"\", \"blue\", palette=\"limited\", showColour=\"background\"))
                               )
                             )
                           )
                         ),
                         fluidRow(actionButton(\"calculate\", \"Calculate/Update\",icon(\"sync\"), 
                           style=\"color: #fff; background-color: #337ab7; border-color: #2e6da4\")
                         ),

                         hr(),

                         fluidRow(
                           column(12,   
                             selectInput(\"vennX\", \"Comparison:\",gsub(\"venn_|.html\",\"\",
                                         list.files(path=\"", venn, "\", pattern=\"all\")), width='600px')
                           )
                         ),
  
                         hr(),
  
                         fluidRow(
                           column(12,htmlOutput('venn'))
                         )                
                ),

                # Tab for Gene Set Enrichment Analysis
                tabPanel(\"GeneSet Enrichment\",
                         fluidRow(
                           column(6,   
                                  selectInput(\"tabGSEA\", \"Comparison:\",gsub(\"\\\\)\",\"\",gsub(\"",PVAL," \\\\(\",\"\",grep(\"",PVAL," \\\\(\",names(tt),value=TRUE)))),
                                  # Update selection options via server.r, as this one reads in the BroadSets...
                                  selectInput(\"tabGSEA_Category\", \"Geneset Category:\",\"All\"),
                                  # Select a particular gene or multiple genes
                                  textInput(\"geneSetByGene\", \"Gene (EntrezID, separate by \',\'):\", value=\"\"),
                                  selectInput(\"cGSEA\", \"Color by factor:\",",TARGETS,"),
                                  selectInput(\"fGSEA\", \"Group by factor:\",",TARGETS,")
                           ),
                           column(width=6, ggvisOutput(\"genSetVolcanoPlot\"))
                         ),
                         
                         hr(),
                         
                         fluidRow(
                           tags$head(tags$style(HTML( \".has-feedback .form-control { padding-right: 0px;}\" ))),
                           column(12,
                                  DT::dataTableOutput(\"genSetTab\")
                           ),
                           hr(),
                           column(12,
                                  DT::dataTableOutput(\"genSetTab_genes\")
                           )
                         ),
                         h3(\"Save pictures:\"),
                         fluidRow(
                           column(
                             width = 12,
                             column(width = 2, numericInput(\"fheight_genset\", \"Height (cm)\", min=2, max=15, step=1, value = 10)),
                             column(width = 2, numericInput(\"fwidth_genset\", \"Width (cm)\", min=2, max=15, step=1, value = 10)),
                             column(width = 2, selectInput(\"fres_genset\", \"Res\", choices=c(\"100\",\"200\",\"300\"), selected = \"100\")),
                             column(width = 2, selectInput(\"fformat_genset\", \"File Type\", choices=c(\"png\",\"tiff\",\"jpeg\",\"pdf\"), 
                                                    selected = \"png\", multiple = FALSE, selectize = TRUE))
                           ),
                           tags$head(tags$style(\".butt{background-color:#337ab7;} .butt{color: #fff;}\")),
                           column( 3, offset=3,
                             downloadButton('gs_download', 'Download Plots', class = \"butt\")
                           )
                           # needed for proper page layout
                           #    cellArgs = list(style=\"width: 400px;\")    
                         ),
                         hr(),
  
                         absolutePanel(
                           top = 60, right = 20, width = 500, draggable = TRUE,
                           wellPanel(
                             HTML(markdownToHTML(fragment.only=TRUE, text=c(  \"Beeswarm plot of the selected Gene (`draggable`)\" ))),
                             plotOutput(\"beeswarm_plot2\", height=\"400px\")
                           ),
                           style = \"opacity: 0.92\"
                         ) #,
                         # absolutePanel(
                         #   top = 20, right = 60, width = 300, draggable = TRUE,
                         #   wellPanel(
                         #     HTML(markdownToHTML(fragment.only=TRUE, text=c(\"Volcano plot of the selected GeneSet (`draggable`)\" ))),
                         #     plotOutput(\"genSetVolcanoPlot\", height=\"200px\")
                         #   ),
                         #   style = \"opacity: 0.92\"
                         # )
                ),
                
                # Tab for Venn diagrams of Gene Sets
                tabPanel(\"Venn diagrams - Gene Sets\",
                         fluidRow(
                           column(2, numericInput(\"pVal_venn_gs\", \"Cut off p-value:\",0.05, min=0.0, max=1.0)),
                           column(4, selectInput(\"fdr_venn_gs\", \"Multiple testing correction\",c(\"none\",\"Benjamini-Hochberg\")))
                         ),
                         fluidRow(
                           column(width = 7,
                           # https://stackoverflow.com/questions/44784083/shiny-label-position-textinput
                             tags$form(
                               class=\"form-horizontal\",
                               tags$div(
                                 class=\"form-group\",
                                 tags$label(class = \"col-sm-4 control-label\", `for` = \"colUp_gs\", br(), \"Color for genes:    Up\"),
                                 column(width = 2, colourInput(\"colUp_gs\", \"\", \"red\", palette=\"limited\", showColour=\"background\")),
                                 tags$label(class = \"col-sm-2 control-label\", `for` = \"colDown_gs\", br(), \"Down\"),
                                 column(width = 2, colourInput(\"colDown_gs\", \"\", \"blue\", palette=\"limited\", showColour=\"background\"))
                               )
                             )
                           )
                         ),
                         fluidRow(actionButton(\"calculate_gs\", \"Calculate/Update\",icon(\"sync\"), 
                           style=\"color: #fff; background-color: #337ab7; border-color: #2e6da4\")
                         ),

                         hr(),

                         fluidRow(
                           column(12,   
                             selectInput(\"vennX_gs\", \"Comparison:\",gsub(\"venn_|.html\",\"\",
                                                    list.files(path=\"", venn_gs, "\", pattern=\"all\")), width='600px')
                             )
                         ),

                         hr(),

                         fluidRow(
                           column(12,htmlOutput('venn_gs'))
                         )                
                ),

                # Generate Beeswarm and Probe expression barplots for selected Genes
                tabPanel(\"Selected Genes\",
                         fluidRow(
                           column(6,   
                                  textInput(\"s\", \"Gene (SYMBOL, EntrezID or Ensembl Gene ID, separate by ','):\", value=\"\")
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput(\"cSel\", \"Color by factor:\",",TARGETS,")
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput(\"fSel\", \"Group by factor:\",",TARGETS,")
                           )
                         ),
                         fluidRow(
                           column(6,   
                                  textInput(\"dir\", \"Files will be saved into <dir>.zip:\",value=\"\"),
                                  downloadButton('downloadPictures', 'Download')                          
                           )
                         )
                )   
    )
  )
  
))
  ")
  
  # Write the txt to the file
  writeLines(txt,fileConn)
  
  # Close the connection
  close(fileConn)

}

writeServer_RNASeq <- function(fileOut="server.R",tt="tt_all.rds", v="voom.rds", 
                               fit="fit.rds", BroadSets="BroadSets.rds", tabCameraFile="tab.camera.v5.1_",
                               object="voom", speciesSymbol="hgnc_symbol",fastqc="fastqc_mergedBAMs",
                               venn = "venn", venn_gs = "venn_gs", scriptDir = "applied_scripts",
                               aheatmap = "aheatmap.png", ...){
  # Open the file to write to
  fileConn <- file(fileOut)
  
  # Read the specs for edgeR & voom objects
  if ( object == "edgeR"){
    PVAL = "PValue"
    COUNTS = "cpm(v, log=TRUE)"
    COUNTS_E = "cpm(v, log=TRUE)"
    COUNTS_Y = "log2 transformed counts"
    DT_INDICES =         "indx = grep(\"^LogCPM$\", colnames(DF))
    indxStart = grep(paste0(\"logFC \\\\(\",input$tabX,\")\"),colnames(DF))
    if (\"^CPM\" %in% colnames(DF)){
      indxEnd = grep(paste0(\"^CPM\"),colnames(DF), value=FALSE)
      indxEnd = indxEnd[1] - 1
    } else {
      indxEnd = dim(DF)[2]
    }
    vec <- c(1:indx,indxStart:(indxEnd))"
#    XCOL = "paste0(\"logFC (\",input$x,\")\")"
#    XCOL_GSEA = "paste0(\"logFC (\",input$tabGSEA,\")\")"
  } else {
    if ( object == "voom"){
      PVAL = "P.Value"
      COUNTS = "v"
      COUNTS_E = "v$E"
      COUNTS_Y = "norm. log2 expressions"
      DT_INDICES =     "if (\"AveExpr\" %in% colnames(DF)){ 
    indx = grep(\"AveExpr\", colnames(DF))
    } else {
      indx = grep(\"^adj.P.Val$\", colnames(DF))
    }
    indxLogFC = grep(paste0(\"logFC \\\\(\",input$tabX,\")\"),colnames(DF))
    indxStart = grep(paste0(\"t \\\\(\",input$tabX,\")\"),colnames(DF))
    indxEnd = grep(paste0(\"adj.P.Val \\\\(\",input$tabX,\")\"),colnames(DF))
    vec <- c(1:indx,indxLogFC,indxStart:indxEnd)"
#      XCOL = "paste0(\"logFC (\",input$x,\")\")"
#      XCOL_GSEA = "paste0(\"logFC (\",input$tabGSEA,\")\")"
# AJ - 0701016     XCOL = "input$x" 
# AJ - 0701016     XCOL_GSEA = "input$tabGSEA"

    }
  }

  # Set the species
  if ( speciesSymbol != ""){
    speciesSYMBOL = speciesSymbol
  }
  if ( speciesSymbol != ""){
    if (speciesSymbol == "mgi_symbol"){
       ENSEMBL_SPECIES = "Mus_musculus"
       ENSEMBL = "ENSMUSG000"
    } else if (speciesSymbol == "hgnc_symbol"){
      ENSEMBL_SPECIES = "Homo_sapiens"
      ENSEMBL = "ENSG000"
    } else if (speciesSymbol == "rgd_symbol"){
      ENSEMBL_SPECIES = "Rattus_norvegicus"
      ENSEMBL = "ENSRNOG000"
    } else {
      ENSEMBL_SPECIES = "Caenorhabditis_elegans"
      ENSEMBL = "ENS" 
    }
  } else {
    ENSEMBL_SPECIES = "Homo_sapiens"
    ENSEMBL = "ENSG000"
  }
  
  # Gather the txt
  txt <- paste0("  
# Check for packages
# http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
if (!requireNamespace(\"BiocManager\") && tools:::.BioC_version_associated_with_R_version() >= 3.8 ){
    install.packages(\"BiocManager\")
}

if (!require(\"pacman\")){ 
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"pacman\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"pacman\", suppressUpdates=TRUE)
  }
}

# install packages using \"pacman\"
pacman::p_load(devtools,curl,shiny,shinythemes,DT,ggvis,ggplot2,beeswarm,plyr,markdown,RColorBrewer,
               randomcoloR,colourpicker,cwhmisc,zip,gplots,
               update = FALSE
       )

# Biobase, limma and edgeR are packages from Bioconductor
if (require(\"Biobase\") == FALSE){ 
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"Biobase\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"Biobase\", suppressUpdates=TRUE)
  }
}
if(require(\"edgeR\") == FALSE) {
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"edgeR\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"edgeR\", suppressUpdates=TRUE)
  }
}
if (require(\"limma\")==  FALSE) {
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"limma\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"limma\", suppressUpdates=TRUE)
  }
}
if (require(\"ComplexHeatmap\") == FALSE) {
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install(\"ComplexHeatmap\")
  } else {
    source(\"http://bioconductor.org/biocLite.R\")
    biocLite(\"ComplexHeatmap\", suppressUpdates=TRUE)
  }
}


# Gleaned from
# - http://www.matteodefelice.name/research/2014/12/29/exploring-time-series-data-with-ggplot2-and-shiny/
# - http://web.stanford.edu/~cengel/cgi-bin/anthrospace/building-my-first-shiny-application-with-ggplot
# 
# The latter worked with reactiveText and reactivePlot, both of which are deprecated
#
# For adapting plot sizes etc.:
# - http://shiny.rstudio.com/articles/images.html

library(shiny)
library(shinythemes)
library(DT)
library(Biobase)
library(ggvis)
library(ggplot2)
library(beeswarm)
library(plyr)
library(limma)
library(edgeR)
library(RColorBrewer)
library(randomcoloR)
library(colourpicker)
library(ComplexHeatmap)
library(zip)

# Extra function for the slider
# https://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale

# logifySlider javascript function
JS.logify <-
\"
// function to logify a sliderInput
function logifySlider (sliderId, sci = false) {
  if (sci) {
    // scientific style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
    })
  } else {
    // regular number style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return (Math.pow(10, num)); }
    })
  }
}\"
                
# call logifySlider for each relevant sliderInput
JS.onload <-
\"
// execute upon document loading
$(document).ready(function() {
  // wait a few ms to allow other scripts to execute
  setTimeout(function() {
    // include call for each slider
    logifySlider('pVal', sci = true)
  }, 5)})
\"

# Define directory that hold functions
scriptDir <- \"", scriptDir, "\"

# Define server for RNASeq
shinyServer(function(input, output, session) {
  
  # Get directory
  pwd=getwd()
  
  # Set paths for QC and Venn HTMLs
  # 20201218 - AJ Changed from \'pwd\' to \'fastqc\'
  if (dir.exists(\"", fastqc, "\" )) {
    addResourcePath(\"pwd\",\"", fastqc,"\")
  } else {
    addResourcePath(\"pwd\",\".\")
  }
  
  # Initialise list holding values with their defaults, both DE genes and GeneSets
  outputZip <- reactiveValues(filename=\"test\", geneSymbol=\"test\", geneEntrezID=1, id=1, setGeneIds=1, geneSets=\"test\")

  #### Data import ####
  DF <- as.data.frame(do.call(readRDS,list(\"",tt,"\")),stringsAsFactors=FALSE)
  v = readRDS(\"",v,"\")
  BroadSets = readRDS(\"",BroadSets,"\")

  # Get the MSigDB categories; obtain the version from the CAMERA table
  # BroadSets variable can also include mouse etc, so more difficult to extract version number??
  # Can I not just get/copy the MSigDBCategory object form the right directory?
  msigdbVersion = gsub(\"_\",\"\", gsub(\"tab.camera.\", \"\",basename(\"",tabCameraFile,"\")))
  if (!file.exists(paste0(dirname(\"",BroadSets,"\"), \"/MSigDB_Categories.\", msigdbVersion,\".rds\"))){
    # Try to get it from the DropBox location (only works for Aldo and Perry....)
    if (exists(\"dropbox\")){
      file.copy(paste0(dropbox,\"/Support/MSigDB/\",msigdbVersion, \"/MSigDB_Categories.\", msigdbVersion,\".rds\"),\".\")
    } else {
      cat(\"Do not have a file with the MSigDB_Categories available\\n\") 
    }
  }

  if (file.exists(paste0(dirname(\"",BroadSets,"\"), \"/MSigDB_Categories.\", msigdbVersion,\".rds\"))){
    MSigDB_Categories = readRDS(paste0(dirname(\"",BroadSets,"\"), \"/MSigDB_Categories.\", msigdbVersion,\".rds\"))
    # Can also set the label and select items
    updateSelectInput(session, \"tabGSEA_Category\", choices=c(\"All\", names(MSigDB_Categories)), selected=\"All\")
  } else {
    # Can also set the label and select items
    updateSelectInput(session, \"tabGSEA_Category\", choices=c(\"Not available\"), selected=\"Not available\")
  }

  if (file.exists(\"",fit,"\")){
    fit = readRDS(\"",fit, "\")
  } else {
    # This is just a standard fitting ....
    fit = lmFit(v, v$design)
  }
  
  if (file.exists(\"",aheatmap,"\")){
    aheatmapFig = \"",aheatmap,"\"
  } else {
    aheatmapFig = \"aheatmap.png\"
  }

  source(paste(scriptDir,\"vennReports.r\", sep=\"/\"), local=TRUE)
  source(paste(scriptDir,\"vennReportsGeneSets.r\", sep=\"/\"), local=TRUE)

  # Initialize list holding the values with their defaults
  values <- reactiveValues(pVal = 0.05, lfc=0.0, adjustM=\"none\", vennDir=\"",venn,"\", vennDir_gs=\"",venn_gs,"\")

  # Initialise
  select <- reactiveValues(id = NULL, symbol = NULL, entrez = NULL, ensembl=NULL, mydf = v$targets)
  
  #### FastQC ####
  output$fastqc <- renderUI({
    tags$iframe(
    seamless=\"seamless\", align=\"center\", frameborder=\"0\", style=\"height:50em; width:100%\", height=\"100%\", scrolling=\"yes\",
    src=\"pwd//fastqc.html\")
  })

  #### MDSplot ####
  pMDS <- reactive(plotMDS(v,top=500, cex=1.0))
    
  mdsColor1 <- reactive({
   Factor1 <- input$f1MDSplot
   Factor1
  })

  mdsColor2 <- reactive({
   Factor2 <- input$f2MDSplot
   Factor2
  })

  output$mdsPlot <- renderPlot({
    Factor1 <- mdsColor1()
    Factor2 <- mdsColor2()

    # Make sure v$targets is a 'normal' data.frame!!
    v$targets <- as.data.frame(v$targets)
    if (Factor1 != Factor2){
      v$targets[,\"Factor\"] <- factor(paste(v$targets[,Factor1], v$targets[,Factor2], sep=\" : \"))
      values$mydf$Factor <- factor(paste(v$targets[,Factor1], v$targets[,Factor2], sep=\" : \"))
    } else {
      v$targets[,\"Factor\"] <- factor(v$targets[,Factor1])
      values$mydf$Factor <- factor(v$targets[,Factor1])
    }
    jColors <- with(v$targets,data.frame(Factor=levels(v$targets[,\"Factor\"]),colors = rainbow(length(levels(v$targets[,\"Factor\"])))))
    colnames(jColors) <- c(\"Factor\",\"colors\")
    v$targets$colors <- NULL
    colorsNew <- plyr::join(v$targets,jColors)
    colLegend <<- rainbow(length(levels(v$targets[,\"Factor\"])))
    colors <- colLegend[v$targets[,\"Factor\"]]

    # make test plot to get the x and y ranges and add/subtract 0.5 for space
    pmds.test <- pMDS()
    x.range <- range(pmds.test$x) + c(-0.5,0.5)
    y.range <- range(pmds.test$y) + c(-0.5,0.5)

    pmds <- pMDS()
    plot(pmds,xlab=\"Leading logFC dim 1\", ylab=\"Leading logFC dim 2\", 
         xlim=x.range,ylim=y.range,pch=16,col=colors, cex=1.0,
         main=\"Top 500 genes\")
    # add legend
    #legend(\"topright\", col=colLegend, 
    #       legend=levels(v$targets[,\"Factor\"]), pch = 16, cex = 1.0, title=Factor)
    
  })
  
  output$mdsPlot_clickinfo <- renderText({
    if (is.null(input$mdsPlot_click$x)) return()
    Factor1 <- mdsColor1()
    Factor2 <- mdsColor2()
    pmds <- pMDS()
    a <- data.frame(xcor=pmds$x,ycor=pmds$y)
    #a$sample <- rownames(a)
    a$sample <- paste0(v$targets[,Factor1],\":\",v$targets[,Factor2])
    colnames(a)<-c(\"xcor\",\"ycor\",\"sample\")
    hit <- a$sample[a$xcor >= (input$mdsPlot_click$x -0.1) & a$xcor <= (input$mdsPlot_click$x + 0.1) & 
                      a$ycor >= (input$mdsPlot_click$y -0.1) & a$ycor <= (input$mdsPlot_click$y + 0.1) ]
    paste0(\"sample:\t\",hit,\"\\n\")
  })
  
  output$legendMDSplot <- renderPlot({
    par(mar = c(0, 0, 1, 1), oma = c(0, 0, 0, 0))
    plot(0, 0, type = \"n\", axes = FALSE, xlab = \"\", ylab = \"\")
    legend(\"top\", col=colLegend, 
           legend=levels(values$mydf$Factor), pch = 16, cex = 1.0, title=\"Condition\")
  })
  
  #### aheatmap ####
  output$aheatmapPlot <- renderImage({
    return(list(
      src = aheatmapFig,
      filetype = \"image/png\",
      width = 600,
      height = 900,
      alt = \"This is a heatmap\"
    ))
  }, deleteFile = FALSE)
  
  #### ComplexHeatmap ####
  if (length(v$distanceMatrix) > 0){
  output$ui_heatmap <- renderUI({
    d3heatmapOutput(\"heatmap\", width = \"100%\", height = \"800px\")
  })

  D3Color <- reactive({
    D3Factor <- input$fD3Heatmap
  })

  # get levels in targets file for d3heatmap factor selection
  # if the ratio of factors/samples > 0.6 then just use the rainbow colors
  # Skip the \"group\", \"lib.size\" and \"norm.factors\" if present

  # Remove columns with NA
  indx <- c(which(is.na(colnames(v$targets))), grep(\"^NA.\", colnames(v$targets)))
  if (length(indx) > 0){
    v$targets <- v$targets[,-indx]
  }

  indx <- grep(c(\"group|lib.size|norm.factors\"), names(v$targets))
  if (length(indx) > 0){
    nrLevels <- apply(subset(v$targets, select=-indx), 2, function(x) length(levels(factor(x))))
  } else {
    nrLevels <- apply(v$targets, 2, function(x) length(levels(factor(x))))
  }
  nrLevels <- nrLevels[which(nrLevels / nrow(v$targets) <= 0.6)]
  n <- sum(nrLevels)
  palette <- distinctColorPalette(n)
  paletteCol <- list()
  n <- 1
  for (i in names(nrLevels)) {
    paletteCol[[i]] <- palette[n:(n + nrLevels[i] - 1)]
    n <- n + nrLevels[i]
  }
  names(paletteCol) <- names(nrLevels)

  # ComplexHeatmap
  complexheatmap <- function()({
    # color for the heatmap
    hmcol = colorRampPalette(brewer.pal(9, \"GnBu\"))(100)
    # get the annotatin rows based on the selected factors
    if (is.null(input$fD3Heatmap) || input$hD3Heatmap) {
      # Create the Heatmap
      ht <- Heatmap(v$distanceMatrix,
                    name           = \"euclidean distance\", 
                    heatmap_legend_param = list(legend_direction = \"horizontal\", legend_width = unit(10, \"cm\")),
                    col            = rev(hmcol),               
                    column_title   = \"Sample-to-sample distances\"
           )

      # Draw the blank heatmap
      # 09042019 - AJ: Original code, but gives error; Error in [.unit: index out of bounds ('unit' subsetting)
      #draw(ht, heatmap_legend_side = \"bottom\", padding = unit(4,\"mm\"))
      draw(ht, heatmap_legend_side = \"bottom\")
    } else {
      annotation = data.frame(v$targets[, D3Color()[1]])
      names(annotation) = D3Color()[1]
      for (i in D3Color()) {
        if (!i %in% colnames(annotation)) {
          annotation[, paste(i)] = factor(v$targets[, i])
        }
      }
  
      # Really make sure all columns are factors!
      # Sometimes reversion to class 'AsIs' was observed and colors would not be drawn
      annotation[] <- lapply(annotation, factor)
  
      # TODO: get colors, should be most distinctive...
      # see http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
        annColors = list()
        for (i in 1:ncol(annotation)) {
          if (names(annotation)[i] %in% names(nrLevels)) {
            annColors[[i]] = paletteCol[[names(annotation)[i]]]
            names(annColors[[i]]) = levels(factor(annotation[,i]))
          } else {
            annColors[[i]] = col2hex(rainbow(length(levels(factor(annotation[, i])))))
            setNames(annColors[[i]], levels(factor(annotation[,i])))
          }
        }
        names(annColors) = colnames(annotation)
  
        annLegendParams = c()
        for (i in colnames(annotation)){
          tmpList = list(list(title_gp=gpar(fontsize=14, fontface=\"bold\"), labels_gp=gpar(fontsize=12)))
          names(tmpList) = i
          annLegendParams = c(annLegendParams,tmpList)
        }
        
        # Get the HeatmapAnnotation object 
        ha <- HeatmapAnnotation(df = annotation, 
                                col = annColors,
                                which = \"column\",
                                annotation_height = unit(rep(0.5, ncol(annotation)), \"cm\"),
                                show_annotation_name = rep(TRUE, ncol(annotation)),
                                annotation_name_gp = gpar(fontface=\"bold\"),
                                annotation_legend_param = annLegendParams
        )
        
        # Create the Heatmap
        ht <- Heatmap(v$distanceMatrix,
                      top_annotation = ha,
                      name           = \"euclidean distance\", 
                      heatmap_legend_param = list(legend_direction = \"horizontal\", legend_width = unit(10, \"cm\")),
                      col            = rev(hmcol),               
                      column_title   = \"Sample-to-sample distances\"
        )
        
        # Draw the heatmap with the desired location for the annotation
        # 09042019 - AJ: Original code, but gives error; Error in [.unit: index out of bounds ('unit' subsetting)
        #draw(ht, heatmap_legend_side = \"bottom\", annotation_legend_side = \"right\", padding = unit(4,\"mm\"))
        draw(ht, heatmap_legend_side = \"bottom\", annotation_legend_side = \"right\")        

        # Set output file name
        outputZip$filename <- paste0(\"heatmap_\",paste0(colnames(annotation),collapse=\"_\"),\".zip\")
     }
  })

  output$heatmap <- renderPlot({ 
    complexheatmap()
  })
    
  ch_downloadname <- reactive({
    tmp = outputZip$filename
    tmp
  })
    
  # download handler
  output$ch_download <- downloadHandler(
    filename = ch_downloadname,
    content = function(fname) {
      tmpdir <- tempdir()
      setwd(tempdir())
      
      fs=c()
      
      # Generate the complex heatmap
      outfile = \"complexHeatmap\"
        
#     fheight <- input$fheight_heatmap
#     fwidth <- input$fwidth_heatmap
      fheight <- 20
      fwidth <- 20
      fres <- as.numeric(input$fres_heatmap)
        
      if(input$fformat_heatmap==\"pdf\") fheight <- round(fheight*0.3937,2)
      if(input$fformat_heatmap==\"pdf\") fwidth <- round(fwidth*0.3937,2)
        
      if(input$fformat_heatmap==\"png\") png(paste0(outfile,\".png\"), height=fheight, width=fwidth, res=fres, units=\"cm\")
      if(input$fformat_heatmap==\"tiff\") tiff(paste0(outfile,\".tiff\"), height=fheight, width=fwidth, res=fres, units=\"cm\",compression=\"lzw\")
      if(input$fformat_heatmap==\"jpeg\") jpeg(paste0(outfile,\".jpeg\"), height=fheight, width=fwidth, res=fres, units=\"cm\",quality=100)
      if(input$fformat_heatmap==\"pdf\") pdf(paste0(outfile,\".pdf\"), height=fheight, width=fwidth)
      complexheatmap()
      dev.off()
        
      fs = paste0(outfile,\".\",input$fformat_heatmap)
        
      zip(zipfile=fname, files=fs)
        
      if(file.exists(paste0(fname, \".zip\"))) {
        file.rename(paste0(fname, \".zip\"), fname)
      }
        
    },
    contentType = \"application/zip\"
    )
  }

  #### Tab: Data Expression Plots ####
  # Make the Volcano plots linked with Beeswarm plots for the genes for the different comparisons
  colorDEP <- reactive({
    cFactor <- input$cDEP
    cFactor
  })

  factorDEP <- reactive({
    Factor <- input$fDEP
    Factor
  })

  DF1 <- reactive({
    # check for the input variable
    yCol = paste0(\"",PVAL," (\",input$x,\")\")
    xCol = paste0(\"logFC (\",input$x,\")\")
    
    dat = data.frame(x = DF[[xCol]], y = -log10(as.numeric(DF[[yCol]])), 
               symbol = as.character(DF[[\"",speciesSYMBOL,"\"]]), 
               entrez = as.character(DF[[\"entrezgene\"]]), ensembl = as.character(DF[[\"ensembl_gene_id\"]]),
               stringsAsFactors = FALSE)
    dat$id = 1:nrow(dat)
    
    dat
  })
  
  click <- function(data, ...) {
    select$id <- data$id
    select$symbol <- data$symbol
    select$entrez <- data$entrez
    select$ensembl <- data$ensembl
    print(data)
    
    treatment <- factor(v$targets[,factorDEP()])
    # Use predefined levels if present
    if (!is.null(levels(v$targets[,factorDEP()]))){
      levels(treatment) <- levels(v$targets[,factorDEP()])
    }

    myLevels <- levels(factor(v$targets[,colorDEP()]))
    myColorsLevels <- rainbow(length(myLevels))
    myColors <- myColorsLevels[factor(v$targets[,colorDEP()])]

    # Beeswarm plot
    output$beeswarm_plot <- renderImage({
      
      # Temporary placeholder
      outfile = tempfile(fileext=\".png\")
      
      # get the clicked gene
      outputZip$geneSymbol <- DF$",speciesSYMBOL,"[select$id]
      outputZip$geneEntrezID <- DF$entrezgene[select$id]
      outputZip$id <- select$id
      outputZip$filename <- paste0(outputZip$geneSymbol,\".zip\")

      # Generate the beeswarm plot
      newTitle = paste0(DF$",speciesSYMBOL,"[select$id],\" (\",DF$entrezgene[select$id],\")\",\" expression\" )
      png(outfile,width=500, height=430)
      par(mar=c(9.1, 4.1, 2.1, 12.1), xpd=TRUE)
      par(mgp=c(2,1,0))
      par(las=2)
      beeswarm(",COUNTS_E,"[select$id,]~treatment,data=",COUNTS_E,",method=\"swarm\",labels=levels(treatment), vertical=TRUE, pch=15, pwcol=myColors, main=newTitle, xlab=factorDEP(), ylab=\"",COUNTS_Y,"\")
      legend(\"topright\",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = colorDEP())
      dev.off()
      
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 500,
           height = 430,
           alt = \"This is alternate text\")
    }, deleteFile = TRUE)
  }
  
  # http://stackoverflow.com/questions/25011544/how-is-data-passed-from-reactive-shiny-expression-to-ggvis-plot  
  
  #  output[['value1']] <- renderText({ names(DF) })
  
  all_values <- function(x) {
    if(is.null(x)) return(NULL)
    row <- DF1()[DF1()$id == x$id,]
    row <- row[,c(3,4,5,1,2)]
    names(row) <- c(\"SYMBOL\",\"EntrezID\",\"EnsemblID\",\"Log2FC\",\"-log10(pvalue)\")
    paste0(names(row), \": \", format(row), collapse = \"<br />\")
  }
  
  # Selected Genes
  DF1.sel <- reactive({ 
    if (input$g == 0) return()
    selList = unlist(strsplit(input$g,\",\"))
    data <- unique(DF1()[DF1()$symbol %in% selList,])
    data <- rbind(data,unique(DF1()[DF1()$entrez %in% selList,]))
    data <- rbind(data,unique(DF1()[DF1()$ensembl %in% selList,]))
    data <- subset(data,data$symbol!=\"NA\")
    
    print(data)
  })
  
  # Selected Genes
  DF1.logFC.sel <- reactive({ 
    if (input$g == 0) return()
    data <- DF1()[abs(DF1()$x) >= input$log2FC & DF1()$y >= -log10(10^input$pVal),]
    
    print(data)
  })
  
  # Volcano plot
  reactive({
    DF1 <- DF1()
    DF1.logFC.sel <- DF1.logFC.sel()
    DF1.sel <- DF1.sel()
    DF1 %>%
    ggvis(~x, ~y) %>%
    handle_click(on_click = click) %>%
    layer_points(key := ~id, size := 30 , size.hover := 500, data = DF1) %>%
    layer_points(key := ~id, size := 30 , stroke := \"deepskyblue\", fill := \"deepskyblue\", size.hover := 500, data = DF1.logFC.sel) %>%
    layer_points(key := ~id, size := 30 , stroke := \"red\", fill := \"red\", size.hover := 500, data = DF1.sel) %>%
    add_axis(\"x\", title=\"log2FoldChange\") %>%
    add_axis(\"y\", title=\"-log10(pvalue)\") %>%
    add_tooltip(all_values, \"hover\") %>%
    layer_paths(stroke:='red', data = data.frame(x=range(DF1()$x),y=-log10(10^input$pVal)) ) %>%
    layer_paths(stroke:='red', data = data.frame(x=-input$log2FC,y=range(DF1()$y)) ) %>%
    layer_paths(stroke:='red', data = data.frame(x=input$log2FC,y=range(DF1()$y)) ) %>%
    set_options(width = 500, height = 400, padding = padding(10, 100, 50, 50)) }) %>%
    bind_shiny(\"volcano_plot\")
  
  default_options()
  
  output$volcano_sel <- renderText({
    paste0(dim(DF1.logFC.sel())[1],\" (out of \", dim(DF1())[1],\") genes fall within the selected cutoffs (log2FC >= abs(\",input$log2FC,\"), p-value <= 1e\",input$pVal,\")\n\")
  })

  output$ggvis_plot <- renderUI({
    ggvisOutput(\"volcano_plot\")
  })

  # download function
  # https://stackoverflow.com/questions/30879084/r-shiny-download-different-image-formats
  #
  fn_downloadname <- reactive({
    print(outputZip$geneSymbol)
    tmp = outputZip$filename
    tmp
  })

  # download handler
  output$bn_download <- downloadHandler(
   filename = fn_downloadname,
   content = function(fname) {
     tmpdir <- tempdir()
     setwd(tempdir())

     fs=c()

     DF1 <- DF1()

     treatment <- factor(v$targets[,factorDEP()])
     # Use predefined levels if present
     if (!is.null(levels(v$targets[,factorDEP()]))){
       levels(treatment) <- levels(v$targets[,factorDEP()])
     }

     myLevels <- levels(factor(v$targets[,colorDEP()]))
     myColorsLevels <- rainbow(length(myLevels))
     myColors <- myColorsLevels[factor(v$targets[,colorDEP()])]

     # Generate the beeswarm plot
     newTitle = paste0(outputZip$geneSymbol,\" (\",outputZip$geneEntrezID,\")\",\" expression\" )
     outfile = paste0(\"beeswarm_\",outputZip$geneSymbol)

     fheight <- input$fheight_volcano
     fwidth <- input$fwidth_volcano
     fres <- as.numeric(input$fres_volcano)

     if(input$fformat_volcano==\"pdf\") fheight <- round(fheight*0.3937,2)
     if(input$fformat_volcano==\"pdf\") fwidth <- round(fwidth*0.3937,2)

     if(input$fformat_volcano==\"png\") png(paste0(outfile,\".png\"), height=fheight, width=fwidth, res=fres, units=\"cm\")
     if(input$fformat_volcano==\"tiff\") tiff(paste0(outfile,\".tiff\"), height=fheight, width=fwidth, res=fres, units=\"cm\",compression=\"lzw\")
     if(input$fformat_volcano==\"jpeg\") jpeg(paste0(outfile,\".jpeg\"), height=fheight, width=fwidth, res=fres, units=\"cm\",quality=100)
     if(input$fformat_volcano==\"pdf\") pdf(paste0(outfile,\".pdf\"), height=fheight, width=fwidth)
       par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
       par(mgp=c(2,1,0))
       par(las=2)
       beeswarm(",COUNTS_E,"[outputZip$id,]~treatment,data=",COUNTS_E,",method=\"swarm\",
                labels=levels(treatment), vertical=TRUE, pch=15, pwcol=myColors, main=newTitle, 
                xlab=factorDEP(), ylab=\"",COUNTS_Y,"\")
       legend(\"topright\",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = colorDEP())
     dev.off()

     # Generate the Volcano plot
     outfile = paste0(\"volcanoPlot_\",outputZip$geneSymbol)

     geneToHighlight <- outputZip$id
     DF1$highlight <- ifelse(abs(DF1$x) >= input$log2FC & DF1$y >= -log10(10^input$pVal),\"highlight\",\"normal\")
     DF1$size <- 30
     DF1$highlight[geneToHighlight] <- \"hit\"
     DF1$size[geneToHighlight] <- 500
     geneToHighlight.df <- DF1[geneToHighlight, ]
     mycolours <- c(\"highlight\"=\"deepskyblue\", \"normal\"=\"black\", \"hit\"=\"red\")

     p <- ggplot(DF1, aes(x=x, y=y)) + geom_point(aes(colour = highlight, size = size))
     p <- p + scale_color_manual(\"Status\", values = mycolours)
     p <- p + labs(x=\"log2FoldChange\", y=\"-log10(pvalue)\", title=input$x)
     p <- p + geom_vline(aes(xintercept=-input$log2FC), color=\"red\")
     p <- p + geom_vline(aes(xintercept=input$log2FC), color=\"red\") 
     p <- p + geom_hline(aes(yintercept=-log10(10^input$pVal)), color=\"red\")
     p <- p + geom_label(data = geneToHighlight.df, aes(x = x, y = y * 0.90, label = symbol), alpha=0.5, size=3)
     p <- p + theme(legend.position = \"none\", plot.title = element_text(hjust = 0.5), plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), \"points\")) + theme()

     if(input$fformat_volcano==\"png\") png(paste0(outfile,\".png\"), height=fheight, width=fwidth, res=fres, units=\"cm\")
     if(input$fformat_volcano==\"tiff\") tiff(paste0(outfile,\".tiff\"), height=fheight, width=fwidth, res=fres, units=\"cm\",compression=\"lzw\")
     if(input$fformat_volcano==\"jpeg\") jpeg(paste0(outfile,\".jpeg\"), height=fheight, width=fwidth, res=fres, units=\"cm\",quality=100)
       plot(p)
     dev.off()

     if(input$fformat_volcano==\"pdf\") ggsave( paste0(outfile,\".pdf\"), plot = p, device=input$fformat_volcano, 
        height=fheight, width=fwidth, dpi=res)

     fs = c(paste0(\"beeswarm_\",outputZip$geneSymbol,\".\",input$fformat_volcano),paste0(\"volcanoPlot_\",outputZip$geneSymbol,\".\",input$fformat_volcano))

     zip(zipfile=fname, files=fs)

     if(file.exists(paste0(fname, \".zip\"))) {
       file.rename(paste0(fname, \".zip\"), fname)
     }

    },
    contentType = \"application/zip\"
  )
  
  #### Tab: Data Table Output ####
  DF_tab <- reactive({
    # get column index adj.P.Val (last column before the individual comparisons start)
    ",DT_INDICES,"

    # get URLs for the EntrezIDs
    indxEntrez = grep(\"^entrez\",colnames(DF))
    urls <- c(paste0(\"http://www.ncbi.nlm.nih.gov/gene/\",as.character(DF[,indxEntrez])))
    refs <- paste0(\"<a href='\",  urls, \"'>\",DF[,indxEntrez],\"</a>\")
    DF[,indxEntrez] <- refs
    
    # get URLs for the Ensembl Gene IDs
    indxEnsembl = grep(\"^ensembl\",colnames(DF))
    urls <- c(paste0(\"http://www.ensembl.org/",ENSEMBL_SPECIES,"/Gene/Summary?db=core;g=\",as.character(DF[,indxEnsembl])))
    refs <- paste0(\"<a href='\",  urls, \"'>\",DF[,indxEnsembl],\"</a>\")
    DF[,indxEnsembl] <- refs

    dat = DF[,vec]
    dat
  })
  
  indxEntrez = grep(\"^entrez\",colnames(DF))
  output$DE_table <- DT::renderDataTable({ datatable(DF_tab(),rownames = TRUE, escape=indxEntrez, options = list(iDisplayLength = 25, search=list(regex = TRUE)), filter=\'top\') })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$tabX, '.csv', sep='') },
    content = function(file) {
      dat = DF_tab()
      indxEntrez = grep(\"^entrez\",colnames(dat))
      noUrls <- gsub(\"<.*'>\",\"\",as.character(dat[,indxEntrez]))
      noUrls <- gsub(\"</a>\",\"\",noUrls)
      dat[,indxEntrez] <- noUrls
      
      indxEnsembl = grep(\"^ensembl\",colnames(dat))
      noUrls <- gsub(\"<.*'>\",\"\",as.character(dat[,indxEnsembl]))
      noUrls <- gsub(\"</a>\",\"\",noUrls)
      dat[,indxEnsembl] <- noUrls

      write.csv(dat, file)
    }
  )
  
  #### Tab: Gene Set Enrichment - Table and plots ####
  # Table for Gene Pathway Enrichment per comparison
  colorGSEA <- reactive({
    cFactor <- input$cGSEA
    cFactor
  })
  
  factorGSEA <- reactive({
    Factor <- input$fGSEA
    Factor
  })
  
  DF2 <- reactive({
    # check for the input variable
    yCol = paste0(\"",PVAL," (\",input$tabGSEA,\")\")
    xCol = paste0(\"logFC (\",input$tabGSEA,\")\")
    
    dat = data.frame(x = DF[[xCol]], y = -log10(as.numeric(DF[[yCol]])), 
                     symbol = as.character(DF[[\"",speciesSYMBOL,"\"]]), 
                     entrez = as.character(DF[[\"entrezgene\"]]), ensembl = as.character(DF[[\"ensembl_gene_id\"]]))
    if (\"",speciesSYMBOL,"\" == \"mgi_symbol\"){
      dat$hgnc_entrez = as.character(DF[[\"entrezgene\"]])
    }

    dat$id = 1:nrow(dat)
    
    dat
  })
  
  click2 <- function(data, ...) {
    select$id <- data$id
    select$symbol <- data$symbol
    select$entrez <- data$entrez
    select$ensembl <- data$ensembl
    print(data)
    
    treatment <- factor(v$targets[,factorGSEA()])
    # Use predefined levels if present
    if (!is.null(levels(v$targets[,factorGSEA()]))){
      levels(treatment) <- levels(v$targets[,factorGSEA()])
    }

    myLevels <- levels(factor(v$targets[,colorGSEA()]))
    myColorsLevels <- rainbow(length(myLevels))
    myColors <- myColorsLevels[factor(v$targets[,colorGSEA()])]

    # Beeswarm plot
    output$beeswarm_plot2 <- renderImage({
      
      # Temporary placeholder
      outfile = tempfile(fileext=\".png\")
      
      # get the clicked gene
      outputZip$geneSymbol <- DF$",speciesSYMBOL,"[select$id]
      outputZip$geneEntrezID <- DF$entrezgene[select$id]
      outputZip$id <- select$id

      # Generate the beeswarm plot
      newTitle = paste0(DF$",speciesSYMBOL,"[select$id],\" (\",DF$entrez[select$id],\")\",\" expression\" )
      png(outfile,width=500, height=430, pointsize=12)
      par(mar=c(9.1, 4.1, 2.1, 12.1), xpd=TRUE)
      par(mgp=c(2,1,0))
      beeswarm(",COUNTS_E,"[select$id,]~treatment,data=",COUNTS_E,",method=\"swarm\",labels=levels(treatment), las=2, cex=c(1.5), vertical=TRUE, pch=15,  pwcol=myColors, main=\"\", xlab=\"\", ylab=\"\", yaxt=\"n\")
      legend(\"topright\",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = colorGSEA())
      mtext(as.expression(substitute(bold(newTitle))), side=3, cex = 2.0)
      mtext(factorGSEA(), side=1, line=2.2, cex=2)
      mtext(\"",COUNTS_Y,"\", side=2, line=2.6, cex=2)
      axis(2,cex.axis=1.5)
      dev.off()
      
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 480,
           height = 400,
           alt = \"This is alternate text\")
    }, deleteFile = TRUE)
  }
  
  output$ggvis_plot <- renderUI({
    ggvisOutput(\"genSetVolcanoPlot\")
  })
  
  all_values2 <- function(x) {
    if(is.null(x)) return(NULL)
    row <- DF2()[DF2()$id == x$id,]
    row <- row[,c(3,4,5,1,2)]
    names(row) <- c(\"SYMBOL\",\"EntrezID\",\"EnsemblID\",\"Log2FC\",\"-log10(pvalue)\")
    paste0(names(row), \": \", format(row), collapse = \"<br />\")
  }
  
  # Need data frames for ablines
  h_abline <- reactive({ data.frame(x=range(DF2()$x),y=-log10(0.05)) })
  v_abline1 <- reactive({ data.frame(x=-log2(2),y=range(DF2()$y)) })
  v_abline2 <- reactive({ data.frame(x=log2(2),y=range(DF2()$y)) })
  
  # Selected Genes
  #  DF2.sel <- reactive({ 
  #    if(input$g != \"\") {
  #      selList = unlist(strsplit(input$g,\",\"))
  #      data <- unique(DF2()[DF2()$symbol %in% selList,])
  #      data <- rbind(data,unique(DF2()[DF2()$entrez %in% selList,]))
  #      data <- subset(data,data$symbol!=\"NA\")
  #    
  #      print(data)
  #    }
  #  })
  
  #### Tab: Geneset Enrichment - Data Table Output ####
  geneSetByGene.sel <- reactive({ 
    if(input$geneSetByGene != \"\") {
      selList = unlist(strsplit(input$geneSetByGene,\",\"))
      selList = input$geneSetByGene

      selList
    }
  })

  genSetTab <- reactive({
    xCol = input$tabGSEA
    
    geneSetCategory = input$tabGSEA_Category
    
    selList <- geneSetByGene.sel()
    
    if ( geneSetCategory == \"All\"){
      if ( length(selList) != 0 ){
        idx_names <- c()
        for (gene_id in selList){
          idx_names <- c(idx_names, names(which(sapply(BroadSets, function(y) gene_id %in% y))))
        }
        tab.camera = read.table(file=paste0(pwd,\"/",tabCameraFile,"\",xCol,\".txt\"),header=TRUE, sep=\"\t\")
        tab.camera = tab.camera[as.character(tab.camera[,1]) %in% idx_names, ]
      } else {
        tab.camera = read.table(file=paste0(pwd,\"/",tabCameraFile,"\",xCol,\".txt\"),header=TRUE, sep=\"\t\")
      }
    } else {
      tab.camera = read.table(file=paste0(pwd,\"/",tabCameraFile,"\",xCol,\".txt\"),header=TRUE, sep=\"\t\")
      tab.camera = tab.camera[gsub(\"^C[1-8]_|^H_\",\"\",as.character(tab.camera[,1])) %in% MSigDB_Categories[[geneSetCategory]], ]
      if ( length(selList) != 0 ){
        idx_names <- c()
        for (gene_id in selList){
          idx_names <- c(idx_names, names(which(sapply(BroadSets, function(y) gene_id %in% y))))
        }
        tab.camera = tab.camera[as.character(tab.camera[,1]) %in% idx_names, ]
        if (!nrow(tab.camera) > 0){
          shiny::showNotification(\"This gene is absent in this category\", type = \"error\")
          NULL
        }

      }      
    }

    # get URLs for the GeneSets@MSigDB and open in new tab
    urls <- c(paste0(\"http://www.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=\",as.character(gsub(\"^C[1-8]_|^H_\",\"\",as.character(tab.camera[,1])))))
    refs <- paste0(\"<a href='\",  urls, \"' target=_blank>\",tab.camera[,1],\"</a>\")
    
    # make the NGenes column clickable to show 
    # - table of genes in the geneset, 
    # - color these genes in the volcano-plot and 
    # - when clicked in the volcano plot display the beeswarm plot
    #    linkGeneSet <- c(paste(sets[[as.character(tab.camera[,1])]],sep=\":\"))
    
    if (nrow(tab.camera)>0){
      #    rownames(tab.camera) <- as.character(tab.camera[,1])
      tab.camera[,1] <- refs
      #    tab.camera[,2] <- linkNGenes
      tab.camera
    } else {
      shiny::showNotification(\"No genesets for this category\", type = \"error\")
      NULL
    }
})
  
  #### Tab: Geneset Enrichment - Beeswarm for download ####
  # Again get the indices from the Broadset coupled with the Voom/EdgeR object?? And read those in?
  # It would mean we provide the Broadsets.rds as well (which might make sense for provenance...)
  output$genSetTab <- DT::renderDataTable(genSetTab(), rownames=FALSE, escape=-c(1), selection = 'single', options = list(iDisplayLength = 10, search=list(regex = TRUE)),filter=\'top\')
  
  click3 <- function(data, ...) {
    # Beeswarm plot
    output$beeswarm_plot2 <- renderImage({
      # Temporary placeholder
      outfile = tempfile(fileext=\".png\")
      
      png(outfile,width=500, height=430, pointsize=12)
      dev.off()
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 250,
           height = 215,
           alt = \"This is alternate text\")
    }, deleteFile = TRUE )
  }
  
  DF2.sel <- reactive({
    selList <- geneSetByGene.sel()
    s = input$genSetTab_row_last_clicked
    if (length(s) & length(selList) == 0) {
      s <- strsplit(strsplit(as.character(genSetTab()[s,1]),\">\")[[1]][2],\"<\")[[1]][1]
      cat('These rows were selected:\\n\t')
      cat(s, sep = '\\n')
      outputZip$filename <- paste0(s,\".zip\")
      outputZip$geneSet <- s
    }
    # This will couple the Entrez IDs to the BroadSets defined earlier
    indx = unlist(BroadSets[s])
    data <- unique(DF2()[DF2()$entrez %in% indx,])
    data <- subset(data,data$entrez!=\"NA\")
    if (\"",speciesSYMBOL,"\" == \"mgi_symbol\"){
      data <- unique(DF2()[DF2()$entrez %in% indx,])
      data <- subset(data,data$entrez!=\"NA\")
    }
    outputZip$setGeneIDs <- data$id
    click3()
    
    data
  })
  
  #### Tab: Geneset Enrichment - Volcano plot for download ####
  DF2a.sel <- reactive({
    selList <- geneSetByGene.sel()
    s = input$genSetTab_row_last_clicked
    if (length(s) & length(selList) == 0) {
      s <- strsplit(strsplit(as.character(genSetTab()[s,1]),\">\")[[1]][2],\"<\")[[1]][1]
      cat('These rows were selected:\\n\t')
      cat(s, sep = '\\n')
    }
    # This will couple the Entrez IDs to the BroadSets defined earlier
    indx = unlist(BroadSets[s])
    data <- unique(DF2()[DF2()$entrez %in% indx,])
    data <- subset(data,data$entrez!=\"NA\")
    data <- data[,c(3,4,1,2)]
    if (\"",speciesSYMBOL,"\" == \"mgi_symbol\"){
      data <- unique(DF2()[DF2()$entrez %in% indx,])
      data <- subset(data,data$entrez!=\"NA\")
    }
    rownames(data) <- data$id
    data$id <- NULL

    DF2 %>%
      ggvis(~x, ~y) %>%
      handle_click(on_click = click2) %>%
      layer_points(key := ~id, size := 30 , size.hover := 500, data = DF2) %>%
      layer_points(key := ~id, size := 30 , stroke := \"red\", fill := \"red\", size.hover := 500, data = DF2.sel) %>%
      add_axis(\"x\", title=\"log2FoldChange\") %>%
      add_axis(\"y\", title=\"-log10(pvalue)\") %>%
      add_tooltip(all_values2, \"hover\") %>%
      layer_paths(stroke:='red', data=h_abline) %>%
      layer_paths(stroke:='red', data=v_abline1) %>%
      layer_paths(stroke:='red', data=v_abline2) %>%
      set_options(width = 400, height = 400, padding = padding(10, 50, 50, 50)) %>%
      bind_shiny(\"genSetVolcanoPlot\")
  
    default_options()
  
    # Generate URL for EntrezID
    if (length(data[,2])) {
      urls <- c(paste0(\"http://www.ncbi.nlm.nih.gov/gene/\",as.character(data[,\"entrez\"])))
      refs <- paste0(\"<a href='\",  urls, \"'>\",data[,\"entrez\"],\"</a>\")
      data[,\"entrez\"] <- refs
    }
    
    indx <- which(colnames(data) %in% c(\"x\", \"y\"))
    colnames(data)[indx] <- c(\"LogFC\", \"-log10(pVal)\")
    
    data
  })
  
  output$genSetTab_genes <- DT::renderDataTable(DF2a.sel(), rownames=FALSE, escape=FALSE, options = list(iDisplayLength = 10, search=list(regex = TRUE)), filter=\'top\')
  # Replace initial VolcanoPlot with a MA-plot?
  
  #### Tab: Geneset Enrichment - Download function ####
  # download function
  # https://stackoverflow.com/questions/30879084/r-shiny-download-different-image-formats
  #
  gn_downloadname <- reactive({
    print(outputZip$geneSet)
    tmp = outputZip$filename
    tmp
  })

  # download handler
  output$gs_download <- downloadHandler(
    filename = gn_downloadname,
    content = function(fname) {
      tmpdir <- tempdir()
      setwd(tempdir())

      fs=c()

      DF1 <- DF2()

      treatment <- factor(v$targets[,factorGSEA()])
      # Use predefined levels if present
      if (!is.null(levels(v$targets[,factorGSEA()]))){
        levels(treatment) <- levels(v$targets[,factorGSEA()])
      }

      myLevels <- levels(factor(v$targets[,colorGSEA()]))
      myColorsLevels <- rainbow(length(myLevels))
      myColors <- myColorsLevels[factor(v$targets[,colorGSEA()])]

      # Generate the beeswarm plot
      newTitle = paste0(outputZip$geneSymbol,\" (\",outputZip$geneEntrezID,\")\",\" expression\" )
      outfile = paste0(\"beeswarm_\",outputZip$geneSymbol)

      fheight <- input$fheight_genset
      fwidth <- input$fwidth_genset
      fres <- as.numeric(input$fres_genset)

      if(input$fformat_genset==\"pdf\") fheight <- round(fheight*0.3937,2)
      if(input$fformat_genset==\"pdf\") fwidth <- round(fwidth*0.3937,2)

      if(input$fformat_genset==\"png\") png(paste0(outfile,\".png\"), height=fheight, width=fwidth, res=fres, units=\"cm\")
      if(input$fformat_genset==\"tiff\") tiff(paste0(outfile,\".tiff\"), height=fheight, width=fwidth, res=fres, units=\"cm\",compression=\"lzw\")
      if(input$fformat_genset==\"jpeg\") jpeg(paste0(outfile,\".jpeg\"), height=fheight, width=fwidth, res=fres, units=\"cm\",quality=100)
      if(input$fformat_genset==\"pdf\") pdf(paste0(outfile,\".pdf\"), height=fheight, width=fwidth)
        par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
        par(mgp=c(2,1,0))
        par(las=2)
        beeswarm(",COUNTS_E,"[outputZip$id,]~treatment,data=",COUNTS_E,",method=\"swarm\",labels=levels(treatment), 
                 vertical=TRUE, pch=15, pwcol=myColors, main=newTitle, xlab=factorDEP(), ylab=\"",COUNTS_Y,"\")
                 legend(\"topright\",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = colorDEP())
      dev.off()

      # Generate the Volcano plot - ToDO: highlight the one gene!
      outfile = paste0(\"volcanoPlot_\",outputZip$geneSet)

      geneSetToHighlight <- outputZip$setGeneIDs
      geneToHighlight <- outputZip$id
      DF1$highlight <- \"normal\"
      DF1$highlight[geneSetToHighlight] <- \"highlight\"
      DF1$size <- 100
      DF1$highlight[geneToHighlight] <- \"hit\"
      DF1$size[geneToHighlight] <- 500
      geneToHighlight.df <- DF1[geneToHighlight, ]
      mycolours <- c(\"highlight\"=\"red\", \"normal\"=\"black\", \"hit\"=\"deepskyblue\")

      DF1_layer1 <- DF1[DF1$highlight == \"normal\",]
      DF1_layer2 <- DF1[DF1$highlight != \"normal\",]
      p <- ggplot(DF1, aes(x=x, y=y)) + geom_point(data=DF1_layer1, aes(x=x,y=y, colour = highlight, size = size))
      p <- p + geom_point(data=DF1_layer2, aes(x=x,y=y, colour = highlight, size = size))
      p <- p + scale_color_manual(\"Status\", values = mycolours)
      p <- p + labs(x=\"log2FoldChange\", y=\"-log10(pvalue)\", title=strwrap(outputZip$geneSet, width=60),subtitle=paste0(\"(\",input$tabGSEA,\")\"))
      p <- p + geom_vline(aes(xintercept=-log2(2)), color=\"gray\")
      p <- p + geom_vline(aes(xintercept=log2(2)), color=\"gray\") 
      p <- p + geom_hline(aes(yintercept=0.05), color=\"gray\")
      p <- p + geom_label(data = geneToHighlight.df, aes(x = x, y = y * 0.85, label = symbol), alpha=0.3, size=3)
      p <- p + theme(legend.position = \"none\", plot.title = element_text(size=14,hjust = 0.5), 
      plot.subtitle = element_text(size=12,hjust = 0.5), plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), \"points\")) + theme()

      if(input$fformat_genset==\"png\") png(paste0(outfile,\".png\"), height=fheight, width=fwidth, res=fres, units=\"cm\")
      if(input$fformat_genset==\"tiff\") tiff(paste0(outfile,\".tiff\"), height=fheight, width=fwidth, res=fres, units=\"cm\",compression=\"lzw\")
      if(input$fformat_genset==\"jpeg\") jpeg(paste0(outfile,\".jpeg\"), height=fheight, width=fwidth, res=fres, units=\"cm\",quality=100)
        plot(p)
      dev.off()

      if(input$fformat_genset==\"pdf\") ggsave( paste0(outfile,\".pdf\"), plot = p, device=input$fformat_genset, 
         height=fheight, width=fwidth, dpi=res)

      fs = c(paste0(\"beeswarm_\",outputZip$geneSymbol,\".\",input$fformat_genset),paste0(\"volcanoPlot_\",outputZip$geneSet,\".\",input$fformat_genset))

      zip(zipfile=fname, files=fs)

      if(file.exists(paste0(\"geneSet_\",fname, \".zip\"))) {
        file.rename(paste0(\"geneSet_\",fname, \".zip\"), fname)
      }

    },
    contentType = \"application/zip\"
  )

  #### Tab: Venn diagrams ####
  # We already have the Venn diagrams with the nicely formatted and linked HTMLs
  # See: http://stackoverflow.com/questions/22177974/r-shiny-using-iframe-for-local-files
  DF_Venn <- reactive({
    # check for the input variable
    comp = paste0(\"venn_\",input$vennX,\"_all.html\")
    comp
  })

  # Observe whether the \"Calculate\" button is pressed and put the new values into the list
  # and make the new Venn diagrams with these new cut-offs.
  # NOTE: off course, we could also have the usr specify the colos for \"up\" and \"down\" -> future?!
  #       https://deanattali.com/2015/06/28/introducing-shinyjs-colourinput/
  #
  myVenn <- \"",venn,"\"
  observe({
   if (input$calculate > 0) {
     pVal_venn <- isolate(input$pVal_venn)
     lfc_venn <- isolate(input$lfc_venn)
     adjustM_venn <- isolate(ifelse(input$adjustM_venn==\"Benjamini-Hochberg\",\"BH\",\"none\"))
     colUp <- isolate(input$colUp)
     colDown <- isolate(input$colDown)

     if (lfc_venn != 0){
       showNotification(\"WARNING: decideTests is used in recalculating the Venn diagrams \n\n
                          Although this function enables users to set p-value and lfc cutoffs simultaneously,\n
                          this combination criterion is not usually recommended. Unless the fold changes and p-values \n
                          are very highly correlated, the addition of a fold change cutoff can increase the \n
                          family-wise error rate or false discovery rate above the nominal level. \n
                          Users wanting to use fold change thresholding are recommended to use treat instead of eBayes \n
                          and to leave lfc at the default value when using decideTests\",
                          duration=NULL, type=\"warning\", closeButton = TRUE
       )
     }

     ### Check if calculation has not already been performed and a directory with results is already present

     # Make sure there is no directory backslash present in the venn_gs directory name...
     myVenn <- ifelse(substring(\"", venn,"\", nchar(\"",venn,"\")) == \"/\", substr(\"", venn,"\", 1, nchar(\"", venn,"\")-1),\"", venn,"\")

     values$vennDir <- paste(myVenn, gsub(\"\\\\.\",\"_\",pVal_venn), gsub(\"\\\\.\",\"_\",lfc_venn), adjustM_venn, sep=\"_\") 
     myVenn <- values$vennDir
     
     if (dir.exists(values$vennDir)){
       print(paste0(\"This calculation has already been performed and results can be found in \", values$vennDir))
     } else {
       print(paste0(\"Results can be found in \",values$vennDir ))

       my.contrasts <- v$contrasts

       withProgress(message = 'Generating Venn diagrams', value = 0.1, {
         for (i in 1:length(vennList)){
           indices = unlist(vennList[i])
           vennReports(fit=fit,contrasts=my.contrasts, indices=indices, myTitle=\"vennDiagram_constrast_1\",
                       reportsDir=values$vennDir, col=c(colUp,colDown), 
                       p.value=pVal_venn, adjust.method=adjustM_venn, lfc=lfc_venn, scriptDir=\"", scriptDir,"\"
           )
           incProgress(length(vennList)/i, message = paste(\"diagram \", i))
         } 
       })
            
       # Move the PNGs to the right place, i.e. not the values$vennDir, but one level up
       Sys.sleep(5)
       toBeMoved <- list.files(\".\", pattern=\"*.png\", full.names=TRUE)
       file.copy(toBeMoved, to=paste0(values$vennDir,\"/../\"))
       file.remove(toBeMoved)
     }
   }
  })

  # Get the name of the directory holding the new Venn diagrams
  vennDir <- renderText({values$vennDir})

  getPage<-function(x) {
    myVenn <- gsub(\"^..\\\\/\", \"\", gsub(\"^.\\\\/\", \"\",gsub(\"\\\\/$\", \"\", values$vennDir)))
    addResourcePath(myVenn, values$vennDir)
    #addResourcePath(\"vennHTML\",vennDir())
    return(
      tags$iframe(
        #src = paste0(\"vennHTML\",\"/venn_\",x,\".html\"), seamless=\"seamless\", width=\"100%\", style=\"height:50em\",  frameborder=\"0\", scrolling=\"yes\"
        src = paste0(myVenn,\"//venn_\",x,\".html\"), seamless=\"seamless\", width=\"100%\", style=\"height:50em\",  frameborder=\"0\", scrolling=\"yes\"
      )
    )
  }
  
  output$venn <- renderUI({
    #    tags$iframe(
    #      seamless=\"seamless\", align=\"center\", frameborder=\"0\", style=\"width:100%\", height=\"100%\", scrolling=\"yes\",
    #      src=paste0(\"pwd/venn/\",DF_Venn())
    #    )
    #x=paste0(\"venn_\",input$vennX,\".html\")
    getPage(input$vennX)
  })

  #### Tab: Venn diagrams - Gene Sets ####
  # We already have the Venn diagrams with the nicely formatted and linked HTMLs
  # See: http://stackoverflow.com/questions/22177974/r-shiny-using-iframe-for-local-files
  DF_Venn_gs <- reactive({
    # check for the input variable
    comp = paste0(\"venn_\",input$vennX_gs,\"_all.html\")
    comp
  })

  # Observe whether the \"Calculate\" button is pressed and put the new values into the list
  # and make the new Venn diagrams with these new cut-offs.
  # NOTE: off course, we could also have the user specify the colors for \"up\" and \"down\" -> future?!
  #       https://deanattali.com/2015/06/28/introducing-shinyjs-colourinput/
  #
  myVennGS <- \"",venn_gs,"\"

  observe({
    if (input$calculate_gs > 0) {
      pVal_venn_gs <- isolate(input$pVal_venn_gs)
      fdr_venn_gs <- isolate(ifelse(input$fdr_venn_gs==\"Benjamini-Hochberg\",\"BH\",\"none\"))
      colUp_gs <- isolate(input$colUp_gs)
      colDown_gs <- isolate(input$colDown_gs)

      if ( is.null(pVal_venn_gs)){ pVal_venn_gs = \"0\" }
      
      # Make sure there is no directory backslash present in the venn_gs directory name...
      myVennGS <- ifelse(substring(\"", venn_gs,"\", nchar(\"",venn_gs,"\")) == \"/\", substr(\"", venn_gs,"\", 1, nchar(\"", venn_gs,"\")-1),\"", venn_gs,"\")

      values$vennDir_gs <- paste(myVennGS, gsub(\"\\\\.\",\"_\",pVal_venn_gs), gsub(\"\\\\.\",\"_\",fdr_venn_gs), sep=\"_\") 
      myVennGS <- values$vennDir_gs
      
      if (dir.exists(values$vennDir_gs)){
        print(paste0(\"This calculation has already been performed and results can be found in \", values$vennDir_gs))
      } else {
        print(paste0(\"Results can be found in \",values$vennDir_gs ))
      
        if ( fdr_venn_gs == \"BH\"){ 
          fdr_venn_gs = pVal_venn_gs
          pVal_venn_gs = NULL 
        }
        my.contrasts <- v$contrasts
      
        if ( fdr_venn_gs == \"none\"){ fdr_venn_gs = NULL } else { pVal_venn_gs = NULL }

        withProgress(message = 'Generating Venn diagrams', value = 0.1, {
          for (i in 1:length(vennList)){
            indices = unlist(vennList[i])
            vennReportsGeneSets(tab=paste0(dirname(\"",v,"\"),\"/tt.genesets.rds\"), filePrefix=paste0(pwd,\"/",tabCameraFile,"\"),
              contrasts=my.contrasts, indices=indices, myTitle=\"vennDiagram_genesets_1\",
              reportsDir=values$vennDir_gs, col=c(colUp_gs,colDown_gs), 
              p.value=pVal_venn_gs, fdr=fdr_venn_gs, scriptDir=\"", scriptDir,"\"
            )
            incProgress(length(vennList)/i, message = paste(\"diagram \", i))
          } 
         })
       
         # Move the PNGs to the right place, i.e. not the values$vennDir_gs, but one level up
         Sys.sleep(5)
         toBeMoved <- list.files(\".\", pattern=\"*.png\", full.names=TRUE)
         file.copy(toBeMoved, to=paste0(values$vennDir_gs,\"/../\"))
         file.remove(toBeMoved)
      }
    }
  })

  # Get the name of the directory holding the new Venn diagrams
  vennDir_gs <- renderText({values$vennDir_gs})

  getPage_gs<-function(x) {
    myVennGS <- gsub(\"^..\\\\/\", \"\", gsub(\"^.\\\\/\", \"\", gsub(\"\\\\/$\", \"\", values$vennDir_gs)))
    addResourcePath(myVennGS, values$vennDir_gs)
    # addResourcePath(\"venn_gs_HTML\",vennDir_gs())
    return(
      tags$iframe(
        #src = paste0(\"venn_gs_HTML\", \"/venn_\",x,\".html\"), seamless=\"seamless\", width=\"100%\", style=\"height:50em\",  frameborder=\"0\", scrolling=\"yes\"
        src = paste0(myVennGS, \"//venn_\",x,\".html\"), seamless=\"seamless\", width=\"100%\", style=\"height:50em\",  frameborder=\"0\", scrolling=\"yes\"
      )
    )
  }

  output$venn_gs <- renderUI({
  #    tags$iframe(
  #      seamless=\"seamless\", align=\"center\", frameborder=\"0\", style=\"width:100%\", height=\"100%\", scrolling=\"yes\",
  #      src=paste0(\"pwd/venn/\",DF_Venn())
  #    )
  #x=paste0(\"venn_\",input$vennX,\".html\")
    getPage_gs(input$vennX_gs)
  })

  
  #### Tab: Selected Genes ####
  # Generate Beeswarm, Probe Expression Barplots for selected genes and save to zip file
  # See: http://stackoverflow.com/questions/28228892/download-multiple-csv-files-in-a-zipped-folder-in-shiny
  factorSel <- reactive({
    Factor <- input$fSel
    Factor
  })

  colorSel <- reactive({
    cFactor <- input$cSel
    cFactor
  })

  DF4.sel <- reactive({ 
    if(input$s != \"\") {
      selList = unlist(strsplit(input$s,\",\"))
      print(selList)
      indxEntrez = grep(\"^entrez\",colnames(DF))
      indxSYMBOL = grep(\"^",speciesSYMBOL,"\",colnames(DF))
      indxEnsembl = grep(\"^ensembl\",colnames(DF))
      data <- unique(DF[DF[,indxSYMBOL[1]] %in% selList,])
      data <- rbind(data,unique(DF[DF[,indxEntrez] %in% selList,]))
      data <- rbind(data,unique(DF[DF[,indxEnsembl] %in% selList,]))
      #      data <- subset(data,data[,indxSYMBOL]!=NA)
      
    }
  })
  
  # We order the labels based on the \"treatment\" factor: just assume that this has been constructed based on the
  # importance of classification
  # R doesn't like \";\" in the columnnames of the design, so these have been replaced by \"_\", but these also have been
  # used to paste different subsets to each other...assuming the first \"_\" to be the most important, replace that by \";\"
  # and split on this
  # Some designs do not contain a \"splitter\", though...
  indx.df = as.data.frame(strsplit(sub(\"_\",\";\",as.character(colnames(v$design))),\";\"))
  nrFactors = length(indx.df[[1]])
  indx.df.t = as.data.frame(t(indx.df))
  indx.df.t$id = 1:nrow(indx.df.t)
  for (i in 1:nrFactors) {
    indx.lab = indx.df.t[order(indx.df.t[,i]),]
  }
  indx.lab = indx.lab$id
  
  # Get a filename for the directory that will hold the pictures (zipped dir)
  zipDir <- reactive({
    if (input$dir != \"\"){
      tmp = paste0(input$dir,\".zip\")
    }
    tmp
  })
  
  output$downloadPictures <- downloadHandler(
    filename = zipDir,
    content = function(fname) {
      tmpdir <- tempdir()
      setwd(tempdir())
      
      treatment <- factor(v$targets[,factorSel()])
      # Use predefined levels if present
      if (!is.null(levels(v$targets[,factorSel()]))){
        levels(treatment) <- levels(v$targets[,factorSel()])
      }

      myLevels <- levels(factor(v$targets[,colorSel()]))
      myColorsLevels <- rainbow(length(myLevels))
      myColors <- myColorsLevels[factor(v$targets[,colorSel()])]

      indxSYMBOL = grep(\"^",speciesSYMBOL,"$\",colnames(DF4.sel()))
      indxEnsembl = grep(\"^ensembl_gene_id$\",colnames(DF4.sel()))

      fs = c()
      for (gene in DF4.sel()[,indxSYMBOL[1]]){
        newTitle = paste0(gene,\" expression\" )
        if(grepl(rownames(v)[1], \"",ENSEMBL,"\" )){
          geneOld = gene
          gene = rownames(",COUNTS,"[grep(paste(DF4.sel()[,indxEnsembl[1]],collapse=\"|\"), v$genes$ensembl_gene_id),])
        }
        png(paste0(\"filteredProbes\",gene,\".png\") )
        par(mar=c(9.1, 4.1, 4.1, 17.1), xpd=TRUE)
        beeswarm(",COUNTS_E,"[grep(paste0(\"^\",gene,\"$\"), v$genes$",speciesSYMBOL,"),]~treatment,data=",COUNTS_E,",method=\"swarm\",labels=levels(treatment), vertical=TRUE, pch=15, pwcol=myColors, main=newTitle, xlab=factorSel(), ylab=\"",COUNTS_Y,"\")
        legend(\"topright\", inset=c(-1.4,0),legend = myLevels, pch = 15, col = myColorsLevels, title = \"Sample\")
        dev.off()
        jpeg(paste0(\"filteredProbes\",gene,\"_barplot.jpg\") )
        par(mar=c(15.1, 4.1, 4.1, 4.1), xpd=TRUE)
#        barplot(",COUNTS,"[grep(paste0(\"^\",gene,\"$\"), rownames(v)),indx.lab],names.arg=with(",COUNTS,",rownames(v$targets)[indx.lab]), main=newTitle, las=2, cex.names=0.8, ylab=\"",COUNTS_Y,"\")
        barplot(",COUNTS_E,"[grep(paste0(\"^\",gene,\"$\"), v$genes$",speciesSYMBOL,"),],names.arg=with(",COUNTS,",rownames(v$targets)), main=newTitle, las=2, cex.names=0.8, ylab=\"",COUNTS_Y,"\")
        dev.off()
        
        fs = c(fs,paste0(\"filteredProbes\",gene,\".png\"), paste0(\"filteredProbes\",gene,\"_barplot.jpg\"))
      }

      if(grepl(\"",ENSEMBL,"\" ,rownames(v)[1])){
        toCSV.dat <- data.frame(GeneName=rownames(",COUNTS_E,"[grep(paste(DF4.sel()[,indxEnsembl[1]],collapse=\"|\"), v$genes$",speciesSYMBOL,"),]),
                              ",COUNTS_E,"[grep(paste(DF4.sel()[,indxEnsembl[1]],collapse=\"|\"), v$genes$",speciesSYMBOL,"),])     
      } else {
        toCSV.dat <- data.frame(GeneName=rownames(",COUNTS_E,"[grep(paste(DF4.sel()[,indxSYMBOL[1]],collapse=\"|\"), v$genes$",speciesSYMBOL,"),]),
                              ",COUNTS_E,"[grep(paste(DF4.sel()[,indxSYMBOL[1]],collapse=\"|\"), v$genes$",speciesSYMBOL,"),])
      }

      write.table(toCSV.dat,file=\"selectedGenes.txt\", row.names=FALSE, quote=FALSE, sep=\"\t\")
      
      fs = c(fs,\"selectedGenes.txt\")

      zip(zipfile=fname, files=fs)
      
      if(file.exists(paste0(fname, \".zip\"))) {
        file.rename(paste0(fname, \".zip\"), fname)
      }
      
    },
    contentType = \"application/zip\"
  )
   
})
  ")
  
  # Write the txt to the file
  writeLines(txt,fileConn)
  
  # Close the connection
  close(fileConn)

}