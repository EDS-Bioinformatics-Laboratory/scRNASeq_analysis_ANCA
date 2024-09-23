  
# Check for packages
# http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
if (!requireNamespace("BiocManager") && tools:::.BioC_version_associated_with_R_version() >= 3.8 ){
    install.packages("BiocManager")
}

if (!require("pacman")){ 
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install("pacman")
  } else {
    source("http://bioconductor.org/biocLite.R")
    biocLite("pacman", suppressUpdates=TRUE)
  }
}

# install packages using "pacman"
pacman::p_load(devtools,curl,shiny,shinythemes,DT,ggvis,ggplot2,beeswarm,plyr,markdown,RColorBrewer,
               randomcoloR,colourpicker,cwhmisc,zip,gplots,
               update = FALSE
       )

# Biobase, limma and edgeR are packages from Bioconductor
if (require("Biobase") == FALSE){ 
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install("Biobase")
  } else {
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biobase", suppressUpdates=TRUE)
  }
}
if(require("edgeR") == FALSE) {
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install("edgeR")
  } else {
    source("http://bioconductor.org/biocLite.R")
    biocLite("edgeR", suppressUpdates=TRUE)
  }
}
if (require("limma")==  FALSE) {
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install("limma")
  } else {
    source("http://bioconductor.org/biocLite.R")
    biocLite("limma", suppressUpdates=TRUE)
  }
}
if (require("ComplexHeatmap") == FALSE) {
  if ( version$major >= 4 || version$minor > 5.1 ){
    BiocManager::install("ComplexHeatmap")
  } else {
    source("http://bioconductor.org/biocLite.R")
    biocLite("ComplexHeatmap", suppressUpdates=TRUE)
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
"
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
}"
                
# call logifySlider for each relevant sliderInput
JS.onload <-
"
// execute upon document loading
$(document).ready(function() {
  // wait a few ms to allow other scripts to execute
  setTimeout(function() {
    // include call for each slider
    logifySlider('pVal', sci = true)
  }, 5)})
"

# Define directory that hold functions
scriptDir <- "applied_scripts/"

# Define server for RNASeq
shinyServer(function(input, output, session) {
  
  # Get directory
  pwd=getwd()
  
  # Set paths for QC and Venn HTMLs
  # 20201218 - AJ Changed from 'pwd' to 'fastqc'
  if (dir.exists("../Results/SummaryTables/QC_fastqc_final" )) {
    addResourcePath("pwd","../Results/SummaryTables/QC_fastqc_final")
  } else {
    addResourcePath("pwd",".")
  }
  
  # Initialise list holding values with their defaults, both DE genes and GeneSets
  outputZip <- reactiveValues(filename="test", geneSymbol="test", geneEntrezID=1, id=1, setGeneIds=1, geneSets="test")

  #### Data import ####
  DF <- as.data.frame(do.call(readRDS,list("14_libra_deg/tt_all.AllComparisons.rds")),stringsAsFactors=FALSE)
  v = readRDS("14_libra_deg/v.AllComparisons.rds")
  BroadSets = readRDS("14_libra_deg/BroadSets_v2023.1.Hs.rds")

  # Get the MSigDB categories; obtain the version from the CAMERA table
  # BroadSets variable can also include mouse etc, so more difficult to extract version number??
  # Can I not just get/copy the MSigDBCategory object form the right directory?
  msigdbVersion = gsub("_","", gsub("tab.camera.", "",basename("14_libra_deg/reports/tab.camera.v2023.1.Hs_")))
  if (!file.exists(paste0(dirname("14_libra_deg/BroadSets_v2023.1.Hs.rds"), "/MSigDB_Categories.", msigdbVersion,".rds"))){
    # Try to get it from the DropBox location (only works for Aldo and Perry....)
    if (exists("dropbox")){
      file.copy(paste0(dropbox,"/Support/MSigDB/",msigdbVersion, "/MSigDB_Categories.", msigdbVersion,".rds"),".")
    } else {
      cat("Do not have a file with the MSigDB_Categories available\n") 
    }
  }

  if (file.exists(paste0(dirname("14_libra_deg/BroadSets_v2023.1.Hs.rds"), "/MSigDB_Categories.", msigdbVersion,".rds"))){
    MSigDB_Categories = readRDS(paste0(dirname("14_libra_deg/BroadSets_v2023.1.Hs.rds"), "/MSigDB_Categories.", msigdbVersion,".rds"))
    # Can also set the label and select items
    updateSelectInput(session, "tabGSEA_Category", choices=c("All", names(MSigDB_Categories)), selected="All")
  } else {
    # Can also set the label and select items
    updateSelectInput(session, "tabGSEA_Category", choices=c("Not available"), selected="Not available")
  }

  if (file.exists("14_libra_deg/fit..AllComparisons.rds")){
    fit = readRDS("14_libra_deg/fit..AllComparisons.rds")
  } else {
    # This is just a standard fitting ....
    fit = lmFit(v, v$design)
  }
  
  if (file.exists("14_libra_deg/aheatmap_initial.png")){
    aheatmapFig = "14_libra_deg/aheatmap_initial.png"
  } else {
    aheatmapFig = "aheatmap.png"
  }

  source(paste(scriptDir,"vennReports.r", sep="/"), local=TRUE)
  source(paste(scriptDir,"vennReportsGeneSets.r", sep="/"), local=TRUE)

  # Initialize list holding the values with their defaults
  values <- reactiveValues(pVal = 0.05, lfc=0.0, adjustM="none", vennDir="14_libra_deg/venn/", vennDir_gs="14_libra_deg/venn_gs/")

  # Initialise
  select <- reactiveValues(id = NULL, symbol = NULL, entrez = NULL, ensembl=NULL, mydf = v$targets)
  
  #### FastQC ####
  output$fastqc <- renderUI({
    tags$iframe(
    seamless="seamless", align="center", frameborder="0", style="height:50em; width:100%", height="100%", scrolling="yes",
    src="pwd//fastqc.html")
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
      v$targets[,"Factor"] <- factor(paste(v$targets[,Factor1], v$targets[,Factor2], sep=" : "))
      values$mydf$Factor <- factor(paste(v$targets[,Factor1], v$targets[,Factor2], sep=" : "))
    } else {
      v$targets[,"Factor"] <- factor(v$targets[,Factor1])
      values$mydf$Factor <- factor(v$targets[,Factor1])
    }
    jColors <- with(v$targets,data.frame(Factor=levels(v$targets[,"Factor"]),colors = rainbow(length(levels(v$targets[,"Factor"])))))
    colnames(jColors) <- c("Factor","colors")
    v$targets$colors <- NULL
    colorsNew <- plyr::join(v$targets,jColors)
    colLegend <<- rainbow(length(levels(v$targets[,"Factor"])))
    colors <- colLegend[v$targets[,"Factor"]]

    # make test plot to get the x and y ranges and add/subtract 0.5 for space
    pmds.test <- pMDS()
    x.range <- range(pmds.test$x) + c(-0.5,0.5)
    y.range <- range(pmds.test$y) + c(-0.5,0.5)

    pmds <- pMDS()
    plot(pmds,xlab="Leading logFC dim 1", ylab="Leading logFC dim 2", 
         xlim=x.range,ylim=y.range,pch=16,col=colors, cex=1.0,
         main="Top 500 genes")
    # add legend
    #legend("topright", col=colLegend, 
    #       legend=levels(v$targets[,"Factor"]), pch = 16, cex = 1.0, title=Factor)
    
  })
  
  output$mdsPlot_clickinfo <- renderText({
    if (is.null(input$mdsPlot_click$x)) return()
    Factor1 <- mdsColor1()
    Factor2 <- mdsColor2()
    pmds <- pMDS()
    a <- data.frame(xcor=pmds$x,ycor=pmds$y)
    #a$sample <- rownames(a)
    a$sample <- paste0(v$targets[,Factor1],":",v$targets[,Factor2])
    colnames(a)<-c("xcor","ycor","sample")
    hit <- a$sample[a$xcor >= (input$mdsPlot_click$x -0.1) & a$xcor <= (input$mdsPlot_click$x + 0.1) & 
                      a$ycor >= (input$mdsPlot_click$y -0.1) & a$ycor <= (input$mdsPlot_click$y + 0.1) ]
    paste0("sample:	",hit,"\n")
  })
  
  output$legendMDSplot <- renderPlot({
    par(mar = c(0, 0, 1, 1), oma = c(0, 0, 0, 0))
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend("top", col=colLegend, 
           legend=levels(values$mydf$Factor), pch = 16, cex = 1.0, title="Condition")
  })
  
  #### aheatmap ####
  output$aheatmapPlot <- renderImage({
    return(list(
      src = aheatmapFig,
      filetype = "image/png",
      width = 600,
      height = 900,
      alt = "This is a heatmap"
    ))
  }, deleteFile = FALSE)
  
  #### ComplexHeatmap ####
  if (length(v$distanceMatrix) > 0){
  output$ui_heatmap <- renderUI({
    d3heatmapOutput("heatmap", width = "100%", height = "800px")
  })

  D3Color <- reactive({
    D3Factor <- input$fD3Heatmap
  })

  # get levels in targets file for d3heatmap factor selection
  # if the ratio of factors/samples > 0.6 then just use the rainbow colors
  # Skip the "group", "lib.size" and "norm.factors" if present

  # Remove columns with NA
  indx <- c(which(is.na(colnames(v$targets))), grep("^NA.", colnames(v$targets)))
  if (length(indx) > 0){
    v$targets <- v$targets[,-indx]
  }

  indx <- grep(c("group|lib.size|norm.factors"), names(v$targets))
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
    hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    # get the annotatin rows based on the selected factors
    if (is.null(input$fD3Heatmap) || input$hD3Heatmap) {
      # Create the Heatmap
      ht <- Heatmap(v$distanceMatrix,
                    name           = "euclidean distance", 
                    heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(10, "cm")),
                    col            = rev(hmcol),               
                    column_title   = "Sample-to-sample distances"
           )

      # Draw the blank heatmap
      # 09042019 - AJ: Original code, but gives error; Error in [.unit: index out of bounds ('unit' subsetting)
      #draw(ht, heatmap_legend_side = "bottom", padding = unit(4,"mm"))
      draw(ht, heatmap_legend_side = "bottom")
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
          tmpList = list(list(title_gp=gpar(fontsize=14, fontface="bold"), labels_gp=gpar(fontsize=12)))
          names(tmpList) = i
          annLegendParams = c(annLegendParams,tmpList)
        }
        
        # Get the HeatmapAnnotation object 
        ha <- HeatmapAnnotation(df = annotation, 
                                col = annColors,
                                which = "column",
                                annotation_height = unit(rep(0.5, ncol(annotation)), "cm"),
                                show_annotation_name = rep(TRUE, ncol(annotation)),
                                annotation_name_gp = gpar(fontface="bold"),
                                annotation_legend_param = annLegendParams
        )
        
        # Create the Heatmap
        ht <- Heatmap(v$distanceMatrix,
                      top_annotation = ha,
                      name           = "euclidean distance", 
                      heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(10, "cm")),
                      col            = rev(hmcol),               
                      column_title   = "Sample-to-sample distances"
        )
        
        # Draw the heatmap with the desired location for the annotation
        # 09042019 - AJ: Original code, but gives error; Error in [.unit: index out of bounds ('unit' subsetting)
        #draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right", padding = unit(4,"mm"))
        draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right")        

        # Set output file name
        outputZip$filename <- paste0("heatmap_",paste0(colnames(annotation),collapse="_"),".zip")
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
      outfile = "complexHeatmap"
        
#     fheight <- input$fheight_heatmap
#     fwidth <- input$fwidth_heatmap
      fheight <- 20
      fwidth <- 20
      fres <- as.numeric(input$fres_heatmap)
        
      if(input$fformat_heatmap=="pdf") fheight <- round(fheight*0.3937,2)
      if(input$fformat_heatmap=="pdf") fwidth <- round(fwidth*0.3937,2)
        
      if(input$fformat_heatmap=="png") png(paste0(outfile,".png"), height=fheight, width=fwidth, res=fres, units="cm")
      if(input$fformat_heatmap=="tiff") tiff(paste0(outfile,".tiff"), height=fheight, width=fwidth, res=fres, units="cm",compression="lzw")
      if(input$fformat_heatmap=="jpeg") jpeg(paste0(outfile,".jpeg"), height=fheight, width=fwidth, res=fres, units="cm",quality=100)
      if(input$fformat_heatmap=="pdf") pdf(paste0(outfile,".pdf"), height=fheight, width=fwidth)
      complexheatmap()
      dev.off()
        
      fs = paste0(outfile,".",input$fformat_heatmap)
        
      zip(zipfile=fname, files=fs)
        
      if(file.exists(paste0(fname, ".zip"))) {
        file.rename(paste0(fname, ".zip"), fname)
      }
        
    },
    contentType = "application/zip"
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
    yCol = paste0("P.Value (",input$x,")")
    xCol = paste0("logFC (",input$x,")")
    
    dat = data.frame(x = DF[[xCol]], y = -log10(as.numeric(DF[[yCol]])), 
               symbol = as.character(DF[["hgnc_symbol"]]), 
               entrez = as.character(DF[["entrezgene"]]), ensembl = as.character(DF[["ensembl_gene_id"]]),
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
      outfile = tempfile(fileext=".png")
      
      # get the clicked gene
      outputZip$geneSymbol <- DF$hgnc_symbol[select$id]
      outputZip$geneEntrezID <- DF$entrezgene[select$id]
      outputZip$id <- select$id
      outputZip$filename <- paste0(outputZip$geneSymbol,".zip")

      # Generate the beeswarm plot
      newTitle = paste0(DF$hgnc_symbol[select$id]," (",DF$entrezgene[select$id],")"," expression" )
      png(outfile,width=500, height=430)
      par(mar=c(9.1, 4.1, 2.1, 12.1), xpd=TRUE)
      par(mgp=c(2,1,0))
      par(las=2)
      beeswarm(v$E[select$id,]~treatment,data=v$E,method="swarm",labels=levels(treatment), vertical=TRUE, pch=15, pwcol=myColors, main=newTitle, xlab=factorDEP(), ylab="norm. log2 expressions")
      legend("topright",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = colorDEP())
      dev.off()
      
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 500,
           height = 430,
           alt = "This is alternate text")
    }, deleteFile = TRUE)
  }
  
  # http://stackoverflow.com/questions/25011544/how-is-data-passed-from-reactive-shiny-expression-to-ggvis-plot  
  
  #  output[['value1']] <- renderText({ names(DF) })
  
  all_values <- function(x) {
    if(is.null(x)) return(NULL)
    row <- DF1()[DF1()$id == x$id,]
    row <- row[,c(3,4,5,1,2)]
    names(row) <- c("SYMBOL","EntrezID","EnsemblID","Log2FC","-log10(pvalue)")
    paste0(names(row), ": ", format(row), collapse = "<br />")
  }
  
  # Selected Genes
  DF1.sel <- reactive({ 
    if (input$g == 0) return()
    selList = unlist(strsplit(input$g,","))
    data <- unique(DF1()[DF1()$symbol %in% selList,])
    data <- rbind(data,unique(DF1()[DF1()$entrez %in% selList,]))
    data <- rbind(data,unique(DF1()[DF1()$ensembl %in% selList,]))
    data <- subset(data,data$symbol!="NA")
    
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
    layer_points(key := ~id, size := 30 , stroke := "deepskyblue", fill := "deepskyblue", size.hover := 500, data = DF1.logFC.sel) %>%
    layer_points(key := ~id, size := 30 , stroke := "red", fill := "red", size.hover := 500, data = DF1.sel) %>%
    add_axis("x", title="log2FoldChange") %>%
    add_axis("y", title="-log10(pvalue)") %>%
    add_tooltip(all_values, "hover") %>%
    layer_paths(stroke:='red', data = data.frame(x=range(DF1()$x),y=-log10(10^input$pVal)) ) %>%
    layer_paths(stroke:='red', data = data.frame(x=-input$log2FC,y=range(DF1()$y)) ) %>%
    layer_paths(stroke:='red', data = data.frame(x=input$log2FC,y=range(DF1()$y)) ) %>%
    set_options(width = 500, height = 400, padding = padding(10, 100, 50, 50)) }) %>%
    bind_shiny("volcano_plot")
  
  default_options()
  
  output$volcano_sel <- renderText({
    paste0(dim(DF1.logFC.sel())[1]," (out of ", dim(DF1())[1],") genes fall within the selected cutoffs (log2FC >= abs(",input$log2FC,"), p-value <= 1e",input$pVal,")
")
  })

  output$ggvis_plot <- renderUI({
    ggvisOutput("volcano_plot")
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
     newTitle = paste0(outputZip$geneSymbol," (",outputZip$geneEntrezID,")"," expression" )
     outfile = paste0("beeswarm_",outputZip$geneSymbol)

     fheight <- input$fheight_volcano
     fwidth <- input$fwidth_volcano
     fres <- as.numeric(input$fres_volcano)

     if(input$fformat_volcano=="pdf") fheight <- round(fheight*0.3937,2)
     if(input$fformat_volcano=="pdf") fwidth <- round(fwidth*0.3937,2)

     if(input$fformat_volcano=="png") png(paste0(outfile,".png"), height=fheight, width=fwidth, res=fres, units="cm")
     if(input$fformat_volcano=="tiff") tiff(paste0(outfile,".tiff"), height=fheight, width=fwidth, res=fres, units="cm",compression="lzw")
     if(input$fformat_volcano=="jpeg") jpeg(paste0(outfile,".jpeg"), height=fheight, width=fwidth, res=fres, units="cm",quality=100)
     if(input$fformat_volcano=="pdf") pdf(paste0(outfile,".pdf"), height=fheight, width=fwidth)
       par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
       par(mgp=c(2,1,0))
       par(las=2)
       beeswarm(v$E[outputZip$id,]~treatment,data=v$E,method="swarm",
                labels=levels(treatment), vertical=TRUE, pch=15, pwcol=myColors, main=newTitle, 
                xlab=factorDEP(), ylab="norm. log2 expressions")
       legend("topright",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = colorDEP())
     dev.off()

     # Generate the Volcano plot
     outfile = paste0("volcanoPlot_",outputZip$geneSymbol)

     geneToHighlight <- outputZip$id
     DF1$highlight <- ifelse(abs(DF1$x) >= input$log2FC & DF1$y >= -log10(10^input$pVal),"highlight","normal")
     DF1$size <- 30
     DF1$highlight[geneToHighlight] <- "hit"
     DF1$size[geneToHighlight] <- 500
     geneToHighlight.df <- DF1[geneToHighlight, ]
     mycolours <- c("highlight"="deepskyblue", "normal"="black", "hit"="red")

     p <- ggplot(DF1, aes(x=x, y=y)) + geom_point(aes(colour = highlight, size = size))
     p <- p + scale_color_manual("Status", values = mycolours)
     p <- p + labs(x="log2FoldChange", y="-log10(pvalue)", title=input$x)
     p <- p + geom_vline(aes(xintercept=-input$log2FC), color="red")
     p <- p + geom_vline(aes(xintercept=input$log2FC), color="red") 
     p <- p + geom_hline(aes(yintercept=-log10(10^input$pVal)), color="red")
     p <- p + geom_label(data = geneToHighlight.df, aes(x = x, y = y * 0.90, label = symbol), alpha=0.5, size=3)
     p <- p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "points")) + theme()

     if(input$fformat_volcano=="png") png(paste0(outfile,".png"), height=fheight, width=fwidth, res=fres, units="cm")
     if(input$fformat_volcano=="tiff") tiff(paste0(outfile,".tiff"), height=fheight, width=fwidth, res=fres, units="cm",compression="lzw")
     if(input$fformat_volcano=="jpeg") jpeg(paste0(outfile,".jpeg"), height=fheight, width=fwidth, res=fres, units="cm",quality=100)
       plot(p)
     dev.off()

     if(input$fformat_volcano=="pdf") ggsave( paste0(outfile,".pdf"), plot = p, device=input$fformat_volcano, 
        height=fheight, width=fwidth, dpi=res)

     fs = c(paste0("beeswarm_",outputZip$geneSymbol,".",input$fformat_volcano),paste0("volcanoPlot_",outputZip$geneSymbol,".",input$fformat_volcano))

     zip(zipfile=fname, files=fs)

     if(file.exists(paste0(fname, ".zip"))) {
       file.rename(paste0(fname, ".zip"), fname)
     }

    },
    contentType = "application/zip"
  )
  
  #### Tab: Data Table Output ####
  DF_tab <- reactive({
    # get column index adj.P.Val (last column before the individual comparisons start)
    if ("AveExpr" %in% colnames(DF)){ 
    indx = grep("AveExpr", colnames(DF))
    } else {
      indx = grep("^adj.P.Val$", colnames(DF))
    }
    indxLogFC = grep(paste0("logFC \\(",input$tabX,")"),colnames(DF))
    indxStart = grep(paste0("t \\(",input$tabX,")"),colnames(DF))
    indxEnd = grep(paste0("adj.P.Val \\(",input$tabX,")"),colnames(DF))
    vec <- c(1:indx,indxLogFC,indxStart:indxEnd)

    # get URLs for the EntrezIDs
    indxEntrez = grep("^entrez",colnames(DF))
    urls <- c(paste0("http://www.ncbi.nlm.nih.gov/gene/",as.character(DF[,indxEntrez])))
    refs <- paste0("<a href='",  urls, "'>",DF[,indxEntrez],"</a>")
    DF[,indxEntrez] <- refs
    
    # get URLs for the Ensembl Gene IDs
    indxEnsembl = grep("^ensembl",colnames(DF))
    urls <- c(paste0("http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",as.character(DF[,indxEnsembl])))
    refs <- paste0("<a href='",  urls, "'>",DF[,indxEnsembl],"</a>")
    DF[,indxEnsembl] <- refs

    dat = DF[,vec]
    dat
  })
  
  indxEntrez = grep("^entrez",colnames(DF))
  output$DE_table <- DT::renderDataTable({ datatable(DF_tab(),rownames = TRUE, escape=indxEntrez, options = list(iDisplayLength = 25, search=list(regex = TRUE)), filter='top') })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$tabX, '.csv', sep='') },
    content = function(file) {
      dat = DF_tab()
      indxEntrez = grep("^entrez",colnames(dat))
      noUrls <- gsub("<.*'>","",as.character(dat[,indxEntrez]))
      noUrls <- gsub("</a>","",noUrls)
      dat[,indxEntrez] <- noUrls
      
      indxEnsembl = grep("^ensembl",colnames(dat))
      noUrls <- gsub("<.*'>","",as.character(dat[,indxEnsembl]))
      noUrls <- gsub("</a>","",noUrls)
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
    yCol = paste0("P.Value (",input$tabGSEA,")")
    xCol = paste0("logFC (",input$tabGSEA,")")
    
    dat = data.frame(x = DF[[xCol]], y = -log10(as.numeric(DF[[yCol]])), 
                     symbol = as.character(DF[["hgnc_symbol"]]), 
                     entrez = as.character(DF[["entrezgene"]]), ensembl = as.character(DF[["ensembl_gene_id"]]))
    if ("hgnc_symbol" == "mgi_symbol"){
      dat$hgnc_entrez = as.character(DF[["entrezgene"]])
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
      outfile = tempfile(fileext=".png")
      
      # get the clicked gene
      outputZip$geneSymbol <- DF$hgnc_symbol[select$id]
      outputZip$geneEntrezID <- DF$entrezgene[select$id]
      outputZip$id <- select$id

      # Generate the beeswarm plot
      newTitle = paste0(DF$hgnc_symbol[select$id]," (",DF$entrez[select$id],")"," expression" )
      png(outfile,width=500, height=430, pointsize=12)
      par(mar=c(9.1, 4.1, 2.1, 12.1), xpd=TRUE)
      par(mgp=c(2,1,0))
      beeswarm(v$E[select$id,]~treatment,data=v$E,method="swarm",labels=levels(treatment), las=2, cex=c(1.5), vertical=TRUE, pch=15,  pwcol=myColors, main="", xlab="", ylab="", yaxt="n")
      legend("topright",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = colorGSEA())
      mtext(as.expression(substitute(bold(newTitle))), side=3, cex = 2.0)
      mtext(factorGSEA(), side=1, line=2.2, cex=2)
      mtext("norm. log2 expressions", side=2, line=2.6, cex=2)
      axis(2,cex.axis=1.5)
      dev.off()
      
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 480,
           height = 400,
           alt = "This is alternate text")
    }, deleteFile = TRUE)
  }
  
  output$ggvis_plot <- renderUI({
    ggvisOutput("genSetVolcanoPlot")
  })
  
  all_values2 <- function(x) {
    if(is.null(x)) return(NULL)
    row <- DF2()[DF2()$id == x$id,]
    row <- row[,c(3,4,5,1,2)]
    names(row) <- c("SYMBOL","EntrezID","EnsemblID","Log2FC","-log10(pvalue)")
    paste0(names(row), ": ", format(row), collapse = "<br />")
  }
  
  # Need data frames for ablines
  h_abline <- reactive({ data.frame(x=range(DF2()$x),y=-log10(0.05)) })
  v_abline1 <- reactive({ data.frame(x=-log2(2),y=range(DF2()$y)) })
  v_abline2 <- reactive({ data.frame(x=log2(2),y=range(DF2()$y)) })
  
  # Selected Genes
  #  DF2.sel <- reactive({ 
  #    if(input$g != "") {
  #      selList = unlist(strsplit(input$g,","))
  #      data <- unique(DF2()[DF2()$symbol %in% selList,])
  #      data <- rbind(data,unique(DF2()[DF2()$entrez %in% selList,]))
  #      data <- subset(data,data$symbol!="NA")
  #    
  #      print(data)
  #    }
  #  })
  
  #### Tab: Geneset Enrichment - Data Table Output ####
  geneSetByGene.sel <- reactive({ 
    if(input$geneSetByGene != "") {
      selList = unlist(strsplit(input$geneSetByGene,","))
      selList = input$geneSetByGene

      selList
    }
  })

  genSetTab <- reactive({
    xCol = input$tabGSEA
    
    geneSetCategory = input$tabGSEA_Category
    
    selList <- geneSetByGene.sel()
    
    if ( geneSetCategory == "All"){
      if ( length(selList) != 0 ){
        idx_names <- c()
        for (gene_id in selList){
          idx_names <- c(idx_names, names(which(sapply(BroadSets, function(y) gene_id %in% y))))
        }
        tab.camera = read.table(file=paste0(pwd,"/14_libra_deg/reports/tab.camera.v2023.1.Hs_",xCol,".txt"),header=TRUE, sep="	")
        tab.camera = tab.camera[as.character(tab.camera[,1]) %in% idx_names, ]
      } else {
        tab.camera = read.table(file=paste0(pwd,"/14_libra_deg/reports/tab.camera.v2023.1.Hs_",xCol,".txt"),header=TRUE, sep="	")
      }
    } else {
      tab.camera = read.table(file=paste0(pwd,"/14_libra_deg/reports/tab.camera.v2023.1.Hs_",xCol,".txt"),header=TRUE, sep="	")
      tab.camera = tab.camera[gsub("^C[1-8]_|^H_","",as.character(tab.camera[,1])) %in% MSigDB_Categories[[geneSetCategory]], ]
      if ( length(selList) != 0 ){
        idx_names <- c()
        for (gene_id in selList){
          idx_names <- c(idx_names, names(which(sapply(BroadSets, function(y) gene_id %in% y))))
        }
        tab.camera = tab.camera[as.character(tab.camera[,1]) %in% idx_names, ]
        if (!nrow(tab.camera) > 0){
          shiny::showNotification("This gene is absent in this category", type = "error")
          NULL
        }

      }      
    }

    # get URLs for the GeneSets@MSigDB and open in new tab
    urls <- c(paste0("http://www.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=",as.character(gsub("^C[1-8]_|^H_","",as.character(tab.camera[,1])))))
    refs <- paste0("<a href='",  urls, "' target=_blank>",tab.camera[,1],"</a>")
    
    # make the NGenes column clickable to show 
    # - table of genes in the geneset, 
    # - color these genes in the volcano-plot and 
    # - when clicked in the volcano plot display the beeswarm plot
    #    linkGeneSet <- c(paste(sets[[as.character(tab.camera[,1])]],sep=":"))
    
    if (nrow(tab.camera)>0){
      #    rownames(tab.camera) <- as.character(tab.camera[,1])
      tab.camera[,1] <- refs
      #    tab.camera[,2] <- linkNGenes
      tab.camera
    } else {
      shiny::showNotification("No genesets for this category", type = "error")
      NULL
    }
})
  
  #### Tab: Geneset Enrichment - Beeswarm for download ####
  # Again get the indices from the Broadset coupled with the Voom/EdgeR object?? And read those in?
  # It would mean we provide the Broadsets.rds as well (which might make sense for provenance...)
  output$genSetTab <- DT::renderDataTable(genSetTab(), rownames=FALSE, escape=-c(1), selection = 'single', options = list(iDisplayLength = 10, search=list(regex = TRUE)),filter='top')
  
  click3 <- function(data, ...) {
    # Beeswarm plot
    output$beeswarm_plot2 <- renderImage({
      # Temporary placeholder
      outfile = tempfile(fileext=".png")
      
      png(outfile,width=500, height=430, pointsize=12)
      dev.off()
      # Return a list containing the filename
      list(src = outfile,
           contentType = 'image/png',
           width = 250,
           height = 215,
           alt = "This is alternate text")
    }, deleteFile = TRUE )
  }
  
  DF2.sel <- reactive({
    selList <- geneSetByGene.sel()
    s = input$genSetTab_row_last_clicked
    if (length(s) & length(selList) == 0) {
      s <- strsplit(strsplit(as.character(genSetTab()[s,1]),">")[[1]][2],"<")[[1]][1]
      cat('These rows were selected:\n	')
      cat(s, sep = '\n')
      outputZip$filename <- paste0(s,".zip")
      outputZip$geneSet <- s
    }
    # This will couple the Entrez IDs to the BroadSets defined earlier
    indx = unlist(BroadSets[s])
    data <- unique(DF2()[DF2()$entrez %in% indx,])
    data <- subset(data,data$entrez!="NA")
    if ("hgnc_symbol" == "mgi_symbol"){
      data <- unique(DF2()[DF2()$entrez %in% indx,])
      data <- subset(data,data$entrez!="NA")
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
      s <- strsplit(strsplit(as.character(genSetTab()[s,1]),">")[[1]][2],"<")[[1]][1]
      cat('These rows were selected:\n	')
      cat(s, sep = '\n')
    }
    # This will couple the Entrez IDs to the BroadSets defined earlier
    indx = unlist(BroadSets[s])
    data <- unique(DF2()[DF2()$entrez %in% indx,])
    data <- subset(data,data$entrez!="NA")
    data <- data[,c(3,4,1,2)]
    if ("hgnc_symbol" == "mgi_symbol"){
      data <- unique(DF2()[DF2()$entrez %in% indx,])
      data <- subset(data,data$entrez!="NA")
    }
    rownames(data) <- data$id
    data$id <- NULL

    DF2 %>%
      ggvis(~x, ~y) %>%
      handle_click(on_click = click2) %>%
      layer_points(key := ~id, size := 30 , size.hover := 500, data = DF2) %>%
      layer_points(key := ~id, size := 30 , stroke := "red", fill := "red", size.hover := 500, data = DF2.sel) %>%
      add_axis("x", title="log2FoldChange") %>%
      add_axis("y", title="-log10(pvalue)") %>%
      add_tooltip(all_values2, "hover") %>%
      layer_paths(stroke:='red', data=h_abline) %>%
      layer_paths(stroke:='red', data=v_abline1) %>%
      layer_paths(stroke:='red', data=v_abline2) %>%
      set_options(width = 400, height = 400, padding = padding(10, 50, 50, 50)) %>%
      bind_shiny("genSetVolcanoPlot")
  
    default_options()
  
    # Generate URL for EntrezID
    if (length(data[,2])) {
      urls <- c(paste0("http://www.ncbi.nlm.nih.gov/gene/",as.character(data[,"entrez"])))
      refs <- paste0("<a href='",  urls, "'>",data[,"entrez"],"</a>")
      data[,"entrez"] <- refs
    }
    
    indx <- which(colnames(data) %in% c("x", "y"))
    colnames(data)[indx] <- c("LogFC", "-log10(pVal)")
    
    data
  })
  
  output$genSetTab_genes <- DT::renderDataTable(DF2a.sel(), rownames=FALSE, escape=FALSE, options = list(iDisplayLength = 10, search=list(regex = TRUE)), filter='top')
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
      newTitle = paste0(outputZip$geneSymbol," (",outputZip$geneEntrezID,")"," expression" )
      outfile = paste0("beeswarm_",outputZip$geneSymbol)

      fheight <- input$fheight_genset
      fwidth <- input$fwidth_genset
      fres <- as.numeric(input$fres_genset)

      if(input$fformat_genset=="pdf") fheight <- round(fheight*0.3937,2)
      if(input$fformat_genset=="pdf") fwidth <- round(fwidth*0.3937,2)

      if(input$fformat_genset=="png") png(paste0(outfile,".png"), height=fheight, width=fwidth, res=fres, units="cm")
      if(input$fformat_genset=="tiff") tiff(paste0(outfile,".tiff"), height=fheight, width=fwidth, res=fres, units="cm",compression="lzw")
      if(input$fformat_genset=="jpeg") jpeg(paste0(outfile,".jpeg"), height=fheight, width=fwidth, res=fres, units="cm",quality=100)
      if(input$fformat_genset=="pdf") pdf(paste0(outfile,".pdf"), height=fheight, width=fwidth)
        par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
        par(mgp=c(2,1,0))
        par(las=2)
        beeswarm(v$E[outputZip$id,]~treatment,data=v$E,method="swarm",labels=levels(treatment), 
                 vertical=TRUE, pch=15, pwcol=myColors, main=newTitle, xlab=factorDEP(), ylab="norm. log2 expressions")
                 legend("topright",inset=c(-0.65,0),legend = myLevels, pch = 15, col = myColorsLevels, title = colorDEP())
      dev.off()

      # Generate the Volcano plot - ToDO: highlight the one gene!
      outfile = paste0("volcanoPlot_",outputZip$geneSet)

      geneSetToHighlight <- outputZip$setGeneIDs
      geneToHighlight <- outputZip$id
      DF1$highlight <- "normal"
      DF1$highlight[geneSetToHighlight] <- "highlight"
      DF1$size <- 100
      DF1$highlight[geneToHighlight] <- "hit"
      DF1$size[geneToHighlight] <- 500
      geneToHighlight.df <- DF1[geneToHighlight, ]
      mycolours <- c("highlight"="red", "normal"="black", "hit"="deepskyblue")

      DF1_layer1 <- DF1[DF1$highlight == "normal",]
      DF1_layer2 <- DF1[DF1$highlight != "normal",]
      p <- ggplot(DF1, aes(x=x, y=y)) + geom_point(data=DF1_layer1, aes(x=x,y=y, colour = highlight, size = size))
      p <- p + geom_point(data=DF1_layer2, aes(x=x,y=y, colour = highlight, size = size))
      p <- p + scale_color_manual("Status", values = mycolours)
      p <- p + labs(x="log2FoldChange", y="-log10(pvalue)", title=strwrap(outputZip$geneSet, width=60),subtitle=paste0("(",input$tabGSEA,")"))
      p <- p + geom_vline(aes(xintercept=-log2(2)), color="gray")
      p <- p + geom_vline(aes(xintercept=log2(2)), color="gray") 
      p <- p + geom_hline(aes(yintercept=0.05), color="gray")
      p <- p + geom_label(data = geneToHighlight.df, aes(x = x, y = y * 0.85, label = symbol), alpha=0.3, size=3)
      p <- p + theme(legend.position = "none", plot.title = element_text(size=14,hjust = 0.5), 
      plot.subtitle = element_text(size=12,hjust = 0.5), plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "points")) + theme()

      if(input$fformat_genset=="png") png(paste0(outfile,".png"), height=fheight, width=fwidth, res=fres, units="cm")
      if(input$fformat_genset=="tiff") tiff(paste0(outfile,".tiff"), height=fheight, width=fwidth, res=fres, units="cm",compression="lzw")
      if(input$fformat_genset=="jpeg") jpeg(paste0(outfile,".jpeg"), height=fheight, width=fwidth, res=fres, units="cm",quality=100)
        plot(p)
      dev.off()

      if(input$fformat_genset=="pdf") ggsave( paste0(outfile,".pdf"), plot = p, device=input$fformat_genset, 
         height=fheight, width=fwidth, dpi=res)

      fs = c(paste0("beeswarm_",outputZip$geneSymbol,".",input$fformat_genset),paste0("volcanoPlot_",outputZip$geneSet,".",input$fformat_genset))

      zip(zipfile=fname, files=fs)

      if(file.exists(paste0("geneSet_",fname, ".zip"))) {
        file.rename(paste0("geneSet_",fname, ".zip"), fname)
      }

    },
    contentType = "application/zip"
  )

  #### Tab: Venn diagrams ####
  # We already have the Venn diagrams with the nicely formatted and linked HTMLs
  # See: http://stackoverflow.com/questions/22177974/r-shiny-using-iframe-for-local-files
  DF_Venn <- reactive({
    # check for the input variable
    comp = paste0("venn_",input$vennX,"_all.html")
    comp
  })

  # Observe whether the "Calculate" button is pressed and put the new values into the list
  # and make the new Venn diagrams with these new cut-offs.
  # NOTE: off course, we could also have the usr specify the colos for "up" and "down" -> future?!
  #       https://deanattali.com/2015/06/28/introducing-shinyjs-colourinput/
  #
  myVenn <- "14_libra_deg/venn/"
  observe({
   if (input$calculate > 0) {
     pVal_venn <- isolate(input$pVal_venn)
     lfc_venn <- isolate(input$lfc_venn)
     adjustM_venn <- isolate(ifelse(input$adjustM_venn=="Benjamini-Hochberg","BH","none"))
     colUp <- isolate(input$colUp)
     colDown <- isolate(input$colDown)

     if (lfc_venn != 0){
       showNotification("WARNING: decideTests is used in recalculating the Venn diagrams 


                          Although this function enables users to set p-value and lfc cutoffs simultaneously,

                          this combination criterion is not usually recommended. Unless the fold changes and p-values 

                          are very highly correlated, the addition of a fold change cutoff can increase the 

                          family-wise error rate or false discovery rate above the nominal level. 

                          Users wanting to use fold change thresholding are recommended to use treat instead of eBayes 

                          and to leave lfc at the default value when using decideTests",
                          duration=NULL, type="warning", closeButton = TRUE
       )
     }

     ### Check if calculation has not already been performed and a directory with results is already present

     # Make sure there is no directory backslash present in the venn_gs directory name...
     myVenn <- ifelse(substring("14_libra_deg/venn/", nchar("14_libra_deg/venn/")) == "/", substr("14_libra_deg/venn/", 1, nchar("14_libra_deg/venn/")-1),"14_libra_deg/venn/")

     values$vennDir <- paste(myVenn, gsub("\\.","_",pVal_venn), gsub("\\.","_",lfc_venn), adjustM_venn, sep="_") 
     myVenn <- values$vennDir
     
     if (dir.exists(values$vennDir)){
       print(paste0("This calculation has already been performed and results can be found in ", values$vennDir))
     } else {
       print(paste0("Results can be found in ",values$vennDir ))

       my.contrasts <- v$contrasts

       withProgress(message = 'Generating Venn diagrams', value = 0.1, {
         for (i in 1:length(vennList)){
           indices = unlist(vennList[i])
           vennReports(fit=fit,contrasts=my.contrasts, indices=indices, myTitle="vennDiagram_constrast_1",
                       reportsDir=values$vennDir, col=c(colUp,colDown), 
                       p.value=pVal_venn, adjust.method=adjustM_venn, lfc=lfc_venn, scriptDir="applied_scripts/"
           )
           incProgress(length(vennList)/i, message = paste("diagram ", i))
         } 
       })
            
       # Move the PNGs to the right place, i.e. not the values$vennDir, but one level up
       Sys.sleep(5)
       toBeMoved <- list.files(".", pattern="*.png", full.names=TRUE)
       file.copy(toBeMoved, to=paste0(values$vennDir,"/../"))
       file.remove(toBeMoved)
     }
   }
  })

  # Get the name of the directory holding the new Venn diagrams
  vennDir <- renderText({values$vennDir})

  getPage<-function(x) {
    myVenn <- gsub("^..\\/", "", gsub("^.\\/", "",gsub("\\/$", "", values$vennDir)))
    addResourcePath(myVenn, values$vennDir)
    #addResourcePath("vennHTML",vennDir())
    return(
      tags$iframe(
        #src = paste0("vennHTML","/venn_",x,".html"), seamless="seamless", width="100%", style="height:50em",  frameborder="0", scrolling="yes"
        src = paste0(myVenn,"//venn_",x,".html"), seamless="seamless", width="100%", style="height:50em",  frameborder="0", scrolling="yes"
      )
    )
  }
  
  output$venn <- renderUI({
    #    tags$iframe(
    #      seamless="seamless", align="center", frameborder="0", style="width:100%", height="100%", scrolling="yes",
    #      src=paste0("pwd/venn/",DF_Venn())
    #    )
    #x=paste0("venn_",input$vennX,".html")
    getPage(input$vennX)
  })

  #### Tab: Venn diagrams - Gene Sets ####
  # We already have the Venn diagrams with the nicely formatted and linked HTMLs
  # See: http://stackoverflow.com/questions/22177974/r-shiny-using-iframe-for-local-files
  DF_Venn_gs <- reactive({
    # check for the input variable
    comp = paste0("venn_",input$vennX_gs,"_all.html")
    comp
  })

  # Observe whether the "Calculate" button is pressed and put the new values into the list
  # and make the new Venn diagrams with these new cut-offs.
  # NOTE: off course, we could also have the user specify the colors for "up" and "down" -> future?!
  #       https://deanattali.com/2015/06/28/introducing-shinyjs-colourinput/
  #
  myVennGS <- "14_libra_deg/venn_gs/"

  observe({
    if (input$calculate_gs > 0) {
      pVal_venn_gs <- isolate(input$pVal_venn_gs)
      fdr_venn_gs <- isolate(ifelse(input$fdr_venn_gs=="Benjamini-Hochberg","BH","none"))
      colUp_gs <- isolate(input$colUp_gs)
      colDown_gs <- isolate(input$colDown_gs)

      if ( is.null(pVal_venn_gs)){ pVal_venn_gs = "0" }
      
      # Make sure there is no directory backslash present in the venn_gs directory name...
      myVennGS <- ifelse(substring("14_libra_deg/venn_gs/", nchar("14_libra_deg/venn_gs/")) == "/", substr("14_libra_deg/venn_gs/", 1, nchar("14_libra_deg/venn_gs/")-1),"14_libra_deg/venn_gs/")

      values$vennDir_gs <- paste(myVennGS, gsub("\\.","_",pVal_venn_gs), gsub("\\.","_",fdr_venn_gs), sep="_") 
      myVennGS <- values$vennDir_gs
      
      if (dir.exists(values$vennDir_gs)){
        print(paste0("This calculation has already been performed and results can be found in ", values$vennDir_gs))
      } else {
        print(paste0("Results can be found in ",values$vennDir_gs ))
      
        if ( fdr_venn_gs == "BH"){ 
          fdr_venn_gs = pVal_venn_gs
          pVal_venn_gs = NULL 
        }
        my.contrasts <- v$contrasts
      
        if ( fdr_venn_gs == "none"){ fdr_venn_gs = NULL } else { pVal_venn_gs = NULL }

        withProgress(message = 'Generating Venn diagrams', value = 0.1, {
          for (i in 1:length(vennList)){
            indices = unlist(vennList[i])
            vennReportsGeneSets(tab=paste0(dirname("14_libra_deg/v.AllComparisons.rds"),"/tt.genesets.rds"), filePrefix=paste0(pwd,"/14_libra_deg/reports/tab.camera.v2023.1.Hs_"),
              contrasts=my.contrasts, indices=indices, myTitle="vennDiagram_genesets_1",
              reportsDir=values$vennDir_gs, col=c(colUp_gs,colDown_gs), 
              p.value=pVal_venn_gs, fdr=fdr_venn_gs, scriptDir="applied_scripts/"
            )
            incProgress(length(vennList)/i, message = paste("diagram ", i))
          } 
         })
       
         # Move the PNGs to the right place, i.e. not the values$vennDir_gs, but one level up
         Sys.sleep(5)
         toBeMoved <- list.files(".", pattern="*.png", full.names=TRUE)
         file.copy(toBeMoved, to=paste0(values$vennDir_gs,"/../"))
         file.remove(toBeMoved)
      }
    }
  })

  # Get the name of the directory holding the new Venn diagrams
  vennDir_gs <- renderText({values$vennDir_gs})

  getPage_gs<-function(x) {
    myVennGS <- gsub("^..\\/", "", gsub("^.\\/", "", gsub("\\/$", "", values$vennDir_gs)))
    addResourcePath(myVennGS, values$vennDir_gs)
    # addResourcePath("venn_gs_HTML",vennDir_gs())
    return(
      tags$iframe(
        #src = paste0("venn_gs_HTML", "/venn_",x,".html"), seamless="seamless", width="100%", style="height:50em",  frameborder="0", scrolling="yes"
        src = paste0(myVennGS, "//venn_",x,".html"), seamless="seamless", width="100%", style="height:50em",  frameborder="0", scrolling="yes"
      )
    )
  }

  output$venn_gs <- renderUI({
  #    tags$iframe(
  #      seamless="seamless", align="center", frameborder="0", style="width:100%", height="100%", scrolling="yes",
  #      src=paste0("pwd/venn/",DF_Venn())
  #    )
  #x=paste0("venn_",input$vennX,".html")
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
    if(input$s != "") {
      selList = unlist(strsplit(input$s,","))
      print(selList)
      indxEntrez = grep("^entrez",colnames(DF))
      indxSYMBOL = grep("^hgnc_symbol",colnames(DF))
      indxEnsembl = grep("^ensembl",colnames(DF))
      data <- unique(DF[DF[,indxSYMBOL[1]] %in% selList,])
      data <- rbind(data,unique(DF[DF[,indxEntrez] %in% selList,]))
      data <- rbind(data,unique(DF[DF[,indxEnsembl] %in% selList,]))
      #      data <- subset(data,data[,indxSYMBOL]!=NA)
      
    }
  })
  
  # We order the labels based on the "treatment" factor: just assume that this has been constructed based on the
  # importance of classification
  # R doesn't like ";" in the columnnames of the design, so these have been replaced by "_", but these also have been
  # used to paste different subsets to each other...assuming the first "_" to be the most important, replace that by ";"
  # and split on this
  # Some designs do not contain a "splitter", though...
  indx.df = as.data.frame(strsplit(sub("_",";",as.character(colnames(v$design))),";"))
  nrFactors = length(indx.df[[1]])
  indx.df.t = as.data.frame(t(indx.df))
  indx.df.t$id = 1:nrow(indx.df.t)
  for (i in 1:nrFactors) {
    indx.lab = indx.df.t[order(indx.df.t[,i]),]
  }
  indx.lab = indx.lab$id
  
  # Get a filename for the directory that will hold the pictures (zipped dir)
  zipDir <- reactive({
    if (input$dir != ""){
      tmp = paste0(input$dir,".zip")
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

      indxSYMBOL = grep("^hgnc_symbol$",colnames(DF4.sel()))
      indxEnsembl = grep("^ensembl_gene_id$",colnames(DF4.sel()))

      fs = c()
      for (gene in DF4.sel()[,indxSYMBOL[1]]){
        newTitle = paste0(gene," expression" )
        if(grepl(rownames(v)[1], "ENSG000" )){
          geneOld = gene
          gene = rownames(v[grep(paste(DF4.sel()[,indxEnsembl[1]],collapse="|"), v$genes$ensembl_gene_id),])
        }
        png(paste0("filteredProbes",gene,".png") )
        par(mar=c(9.1, 4.1, 4.1, 17.1), xpd=TRUE)
        beeswarm(v$E[grep(paste0("^",gene,"$"), v$genes$hgnc_symbol),]~treatment,data=v$E,method="swarm",labels=levels(treatment), vertical=TRUE, pch=15, pwcol=myColors, main=newTitle, xlab=factorSel(), ylab="norm. log2 expressions")
        legend("topright", inset=c(-1.4,0),legend = myLevels, pch = 15, col = myColorsLevels, title = "Sample")
        dev.off()
        jpeg(paste0("filteredProbes",gene,"_barplot.jpg") )
        par(mar=c(15.1, 4.1, 4.1, 4.1), xpd=TRUE)
#        barplot(v[grep(paste0("^",gene,"$"), rownames(v)),indx.lab],names.arg=with(v,rownames(v$targets)[indx.lab]), main=newTitle, las=2, cex.names=0.8, ylab="norm. log2 expressions")
        barplot(v$E[grep(paste0("^",gene,"$"), v$genes$hgnc_symbol),],names.arg=with(v,rownames(v$targets)), main=newTitle, las=2, cex.names=0.8, ylab="norm. log2 expressions")
        dev.off()
        
        fs = c(fs,paste0("filteredProbes",gene,".png"), paste0("filteredProbes",gene,"_barplot.jpg"))
      }

      if(grepl("ENSG000" ,rownames(v)[1])){
        toCSV.dat <- data.frame(GeneName=rownames(v$E[grep(paste(DF4.sel()[,indxEnsembl[1]],collapse="|"), v$genes$hgnc_symbol),]),
                              v$E[grep(paste(DF4.sel()[,indxEnsembl[1]],collapse="|"), v$genes$hgnc_symbol),])     
      } else {
        toCSV.dat <- data.frame(GeneName=rownames(v$E[grep(paste(DF4.sel()[,indxSYMBOL[1]],collapse="|"), v$genes$hgnc_symbol),]),
                              v$E[grep(paste(DF4.sel()[,indxSYMBOL[1]],collapse="|"), v$genes$hgnc_symbol),])
      }

      write.table(toCSV.dat,file="selectedGenes.txt", row.names=FALSE, quote=FALSE, sep="	")
      
      fs = c(fs,"selectedGenes.txt")

      zip(zipfile=fname, files=fs)
      
      if(file.exists(paste0(fname, ".zip"))) {
        file.rename(paste0(fname, ".zip"), fname)
      }
      
    },
    contentType = "application/zip"
  )
   
})
  
