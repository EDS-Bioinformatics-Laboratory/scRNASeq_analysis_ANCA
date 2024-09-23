
# Check for packages and suppress updating packages
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
    logifySlider('log_slider', sci = false)
    logifySlider('pVal', sci = true)
  }, 5)})
"

# check existence of needed files - R 3.2.0
#if(!dir.exists("14_libra_deg/venn/")){
#  cat("No Venn information found, directory 'venn' does not exist!\n")
#}
                
# check existence of needed files - pre R 3.2.0
if(dir.exists("14_libra_deg/venn/")  == FALSE) {
  cat("No Venn information found, directory '14_libra_deg/venn/' does not exist!\n")
}

# See whether there is information about the experiment
if(file.exists("14_libra_deg/ExpInfo.txt")){
  expInfo <- 1
} else {
  cat("14_libra_deg/ExpInfo.txt could not be found\n")
}
                
# read in tt
if(file.exists("14_libra_deg/tt_all.AllComparisons.rds")){
  tt <- readRDS("14_libra_deg/tt_all.AllComparisons.rds")
} else {
  cat("14_libra_deg/tt_all.AllComparisons.rds could not be found\n")
}

# read in v
if(file.exists("14_libra_deg/v.AllComparisons.rds")){
  v <- readRDS("14_libra_deg/v.AllComparisons.rds")
  indx <- grep(c("group|lib.size|norm.factors"), colnames(v$targets))
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
    presentVenns <- c(dir("14_libra_deg/venn/",pattern="*.png"))
    presentVenns <- gsub("^venn_","",gsub("\\.png","",presentVenns))
    # Extract the names, count the "_vs_"s and split on 2nd, 4th etc.
    indx <- as.data.frame(strsplit2(presentVenns,"_vs_|_and_"))
    # Make combinations
    indx1 <- apply(indx,1,function(x) paste(x[seq(1,length(x),1)]))
    if ( dim(indx)[2] >= 2) {
      indx2 <- apply(indx,1,function(x) paste(x[seq(1,length(x),2)],x[seq(2,length(x),2)], sep="_vs_"))
      indxAll <- rbind(indx1,indx2)
      if ( dim(indx)[2] >= 3) {
         indx3 <- apply(indx,1,function(x) paste(x[seq(1,length(x),2)],x[seq(2,length(x),2)], x[seq(3,length(x),3)],sep="_vs_"))
         indxAll <- rbind(indxAll,indx3)
         if ( dim(indx)[2] >= 4) {
            indx4 <- apply(indx,1,function(x) paste(x[seq(1,length(x),2)],x[seq(2,length(x),2)], x[seq(3,length(x),3)], x[seq(4,length(x),4)],sep="_vs_"))
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
  cat("14_libra_deg/v.AllComparisons.rds could not be found\n")
}

# define title
title <- "RNASeq Analysis"

# Define UI for RNASeq application
shinyUI(fluidPage(
  # Theme (https://github.com/ClintWeathers/stuff_about_r/blob/master/datatable_links_plus_theme.md)
  theme = shinytheme("cosmo"), 
  
  # Set the color of the Warning/Error messages
  tags$head(
    tags$style(HTML("
      .shiny-output-error-validation { color: blue; }
    "))
  ),

  # Application title
  headerPanel(title),
  
  mainPanel(
    tabsetPanel(type="tabs",
                # Tab for Introduction of app
                tabPanel("Introduction",
                    fluidRow(column(
                      h2("Welcome to the RNASeq Analysis Shiny App!"),
                      br(),
                      tags$p(HTML(
                       "The RNASeq Analysis Shiny App takes the R objects resulting from an analysis performed on your RNASeq data and 
                        allows interactive visualization and inspection of the obtained results.
The results are organized using different tabs:"
                      )),
                      tags$p(HTML(
                        '<ol start=\'1\'> 
                         <li><font color="blue">"Exp. Info"</font> shows the experimental information regarding your data, the analysis steps performed, some observations regarding the data and possibly some conclusions.</li>
                         <li><font color="blue">"FastQC of final BAM files"</font> <i><b>(if present)</b></i> displays the results of the FastQC analysis of the BAM files obtained after alignment. Each colored tile is clickable and will show you the results of the analysis for that sample for that particular statistic. By clicking on the links in the rows starting with "Collated ..", you can obtain the information on that statistic for all samples. This way you can easily find possible outliers on sample level.</li>
                         <li><font color="blue">"After normalization & filtering"</font> shows you a multi-dimensional scaling (MDS) plot in which you can color samples based on 2 different factors to identify groups of samples. You can also click on samples to get their name. A static heatmap ("Heatmap" tab) is shown and a somewhat more dynamic one ("iHeatmap" tab) in which you can specify which factors you want to show..</li>
                         <li><font color="blue">"Differential Expression Plots"</font> gives you the results of the differentially expressed genes (DEG) analysis. Via a drop-down menu you can select which comparison you want to investigate. The corresponding "volcano plot" is then loaded and you can click on the dots to obtain a plot of the expression values of that particular gene in the various samples. The "Color by Factor" and "Group by Factor" menus allow you to change the appearance of this "beeswarm plot" accordingly. You can also display your pet gene by selecting it by name using different types of identiers. Finally you can save the plots in different formats and resolutions for use in presentations.
<b><font color="red">The generated figures are by no means publication ready and should therefore definitely not be used as such !!</font></b></li>
                         <li><font color="blue">"Differential Expression Table"</font> shows you the data underlying all plots. Again, you can specify the comparison you are interested in. You can sort on column information by clicking the small arrows displayed next to the column header. You can filter by clicking in the boxes just below the column headers. This will give you the filter options. You can also click on gene identifiers to obtain more information at the NCBI or Ensembl webpages.</li>
                         <li><font color="blue">"Venn diagrams - Genes"</font> allows you to see the intersections and differences of genes significantly expressed in different comparisons. The initial plot has been pre-calculated, but you can change this to your liking using different colors and/or cut-offs. The numbers in the plot are clickable and will show the corresponding genes. Genes have to show the same sign of expression up/down regulation in order to be counted as overlapping. For a statistically sound analysis it is advised to use the Benjamin-Hochberg multiple testing correction and to not use a cut-off on the log2 fold change  (that is, stick to the default value of 0).</li>
                         <li><font color="blue">"Geneset Enrichment"</font> shows the result of a CAMERA <a href="https://academic.oup.com/nar/article/40/17/e133/2411151"> (Wu & Smyth, 2012)</a> analysis using genesets obtained from the <a href="http://software.broadinstitute.org/gsea/msigdb/index.jsp">MSigDB</a>. You can again select the comparison of interest. By clicking on the name of the geneset, you will be led to the corresponding page at the MSigDB site. Here you can obtain much more information on the geneset. Just clicking the row of the geneset will highlight the genes in that particular geneset in the volcanoplot (in red). Information on these genes will be displayed below the geneset table. You can then again click on these genes to get the corresponding "beeswarm plot". Coloring and grouping can be performed as well. Selection of specific "Geneset Categories" is also an option. One can also provide the EntrezGene ID of gene(s) and get only the genesets these occur in.</li> 
                         <li><font color="blue">"Venn diagrams - Gene Sets"</font> allows you to see the intersections and differences of genesets significantly enriched in different comparisons. The initial plot has been pre-calculated, but you can change this to your liking using different colours and/or cut-offs. The numbers in the plot are clickable and will show the corresponding gene sets. Gene sets have to show the same sign of up/down regulation in order to be counted as overlapping.</li>
                         <li><font color="blue">"Selected Genes"</font> finally will allow you to obtain some information on a gene or set of genes you are interested in. It currently consists of a barplot and beeswarm plot of the expression.</li>
                         </ol>'
                      )),
                      br(),
                      tags$h4(HTML('<font color="red">This Shiny app is always under construction...!!</font>')),
                      br(),
                      tags$p(HTML(
                        'For questions, bug reporting, feature requests, remarks, comments, installation problems etc. please do not hesitate to
                                  <a href="mailto:a.jongejan@amc.uva.nl">contact me!</a>'
                      )),
                      width = 12)
                    )
                ),

                # Tab for Experimental Info
         tabPanel("Exp. Info",
                  fluidRow(
                    column(12,
                      pre(includeText("14_libra_deg/ExpInfo.txt"))     
                    )
                  )        
        ),

                # Tab for Fastqc files
                tabPanel("FastQC of final BAM files",
                    fluidRow(
                      column(12,
                        htmlOutput('fastqc')      
                      )
                    )                
                ),
                # Tab for MDS and heatmap
                tabPanel("After normalization\n& filtering",
                         tabsetPanel(
                           tabPanel("MDS plot",
                                    fluidRow(
                                      column(6,
                                        selectInput("f1MDSplot", "Color by factor 1:",factorNames)
                                      )
                                    ),
                                    fluidRow(
                                      column(6,
                                        selectInput("f2MDSplot", "Color by factor 2:",factorNames)
                                      )
                                    ),
                                    hr(),
                                    fluidRow(
                                      column(width = 7,
                                        plotOutput('mdsPlot', click="mdsPlot_click", width=550, height=550)
                                      ),    
                                      column(width = 5, offset=0.3,
                                        verbatimTextOutput("mdsPlot_clickinfo")
                                      )
                                    ),
                                    absolutePanel(
                                      bottom = 120, right = 20, width = paste0(legendWidth*9,"px"), height = paste0(legendHeight*18,"px"), draggable = TRUE,
                                      wellPanel(
                                        HTML(markdownToHTML(fragment.only=TRUE, text=c(  "Legend (`draggable`)" ))),
                                        plotOutput("legendMDSplot", height=paste0(legendHeight*18,"px"))
                                      ),
                                      style = "opacity: 0.92"
                                    )
                           ),
                           tabPanel("Heatmap",
                                    plotOutput('aheatmapPlot')
                           ),
                           tabPanel("iHeatmap",
                                    fluidRow(
                                      checkboxGroupInput("fD3Heatmap", "Color by factor:",names(nrLevels), selected = names(nrLevels)[1],inline = TRUE),
                                      checkboxInput("hD3Heatmap",label="Hide factors", value=FALSE, width=NULL)
                                    ),
                                    hr(),
                                    #  uiOutput("ui_heatmap")
                                    plotOutput("heatmap",height = 800),
                                    h3("Save pictures:"),
                                    fluidRow(
                                      column(
                                        width = 12,
                                        column(width = 2, selectInput("fres_heatmap", "Res", choices=c("100","200","300"), selected = "100")),
                                        column(width = 2, selectInput("fformat_heatmap", "File Type", choices=c("png","tiff","jpeg","pdf"), 
                                                                      selected = "png", multiple = FALSE, selectize = TRUE))
                                      ),
                                      tags$head(tags$style(".butt{background-color:#337ab7;} .butt{color: #fff;}")),
                                      column( 3, offset=3,
                                              downloadButton('ch_download', 'Download Plots', class = "butt")
                                      )
                                    ),
                                    hr()
                              )
                         )
                ),          
                
                # Tab for Differential Expression Plots
                tabPanel("Differential Expression Plots",
                         fluidRow(
                           column(6,   
                              selectInput("x", "Comparison:",gsub(")","",gsub("P.Value \\(","",grep("P.Value \\(",names(tt),value=TRUE))), width="600px")
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput("cDEP", "Color by factor:",factorNames)
                           ),
                           column(6,
                             sliderInput("log2FC", "Log2FC cutoff:", 0.5, 10, 1, step = NULL, round = FALSE, ticks = TRUE, animate = FALSE, width = NULL, sep = ",", pre = NULL, post = NULL, timeFormat = NULL, timezone = NULL, dragRange = TRUE)
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput("fDEP", "Group by factor:",factorNames)
                           ),
                           # https://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale
                           tags$head(tags$script(HTML(JS.logify))),
                           tags$head(tags$script(HTML(JS.onload))),
                           tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode", function(message) { eval(message.value); });'))),
                           column(6,
#                             sliderInput("pVal", "p-value cutoff:", 0.0000000001, 0.1, 0.05, step = NULL, round = FALSE, ticks = TRUE, animate = FALSE, width = NULL, sep = , pre = NULL, post = NULL, timeFormat = NULL, timezone = NULL, dragRange = TRUE)

                           sliderInput("pVal", "pVal cutoff (log10):", min = -10, max = 0, value = -0.5, step = 0.5)
                           )
                         ),
                         
                         hr(),
                         
                         fluidRow(
                           column(6,
                                  textInput("g", "Gene (SYMBOL, EntrezID or Ensembl Gene ID, separate by ','):", value = "1")
                           )
                         ),
                         
                         hr(),
                         
                         fluidRow(
                           verbatimTextOutput("volcano_sel"),                           
                           column(7,div(style = "height:500px;width:550px;",
                                        ggvisOutput("volcano_plot"))
                           ),
                           column(5,div(style = "height:500px;width:500px;",
                                        plotOutput("beeswarm_plot",width=400, height=400) )       
                           )
                           # needed for proper page layout
                           #    cellArgs = list(style="width: 400px;")    
                         ),
                         h3("Save pictures:"),
                         fluidRow(
                           column(
                             width = 12,
                             column(width = 2, numericInput("fheight_volcano", "Height (cm)", min=2, max=15, step=1, value = 10)),
                             column(width = 2, numericInput("fwidth_volcano", "Width (cm)", min=2, max=15, step=1, value = 10)),
                             column(width = 2, selectInput("fres_volcano", "Res", choices=c("100","200","300"), selected = "100")),
                             column(width = 2, selectInput("fformat_volcano", "File Type", choices=c("png","tiff","jpeg","pdf"), 
                                                    selected = "png", multiple = FALSE, selectize = TRUE))
                           ),
                           tags$head(tags$style(".butt{background-color:#337ab7;} .butt{color: #fff;}")),
                           column( 3, offset=3,
                             downloadButton('bn_download', 'Download Plots', class = "butt")
                           )
                           # needed for proper page layout
                           #    cellArgs = list(style="width: 400px;")    
                         ),
                         hr()
  
                ),
                
                # Tab for Differential Expression Table
                tabPanel("Differential Expression Table",
                         fluidRow(
                           column(6,   
                                  selectInput("tabX", "Comparison:",gsub(")","",gsub("P.Value \\(","",grep("P.Value \\(",names(tt),value=TRUE))), width='600px'),
                                  downloadButton('downloadData', 'Download')            
                           )
                         ),
                         
                         hr(),
                         
                         fluidRow(
                           # see https://stackoverflow.com/questions/34863466/data-table-column-filters-not-displaying-properly-in-shiny-using-dt-package
                           # The server=FALSE does not work ...
                           tags$head(tags$style(HTML( ".has-feedback .form-control { padding-right: 0px;}" ))),
                           column(12,
                                  DT::dataTableOutput("DE_table")
                                  
                           )
                         )
                ),
                
                # Tab for Venn diagrams - Genes
                tabPanel("Venn diagrams - Genes",
                         fluidRow(
                           column(2, numericInput("pVal_venn", "Cut off p-value:",0.05)),
                           column(2, numericInput("lfc_venn", "log Fold-change:",value=0.0, min=0.0)),
                           column(3, selectInput("adjustM_venn", "Multiple testing correction",c("none","Benjamini-Hochberg")))
                         ),
                         fluidRow(
                           column(
                             width = 7,
                             # https://stackoverflow.com/questions/44784083/shiny-label-position-textinput
                             tags$form(
                               class="form-horizontal",
                               tags$div(
                                 class="form-group",
                                 tags$label(class = "col-sm-4 control-label", `for` = "colUp", br(), "Color for genes:    Up"),
                                 column(width = 2, colourInput("colUp", "", "red", palette="limited", showColour="background")),
                                 tags$label(class = "col-sm-2 control-label", `for` = "colDown", br(), "Down"),
                                 column(width = 2, colourInput("colDown", "", "blue", palette="limited", showColour="background"))
                               )
                             )
                           )
                         ),
                         fluidRow(actionButton("calculate", "Calculate/Update",icon("sync"), 
                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                         ),

                         hr(),

                         fluidRow(
                           column(12,   
                             selectInput("vennX", "Comparison:",gsub("venn_|.html","",
                                         list.files(path="14_libra_deg/venn/", pattern="all")), width='600px')
                           )
                         ),
  
                         hr(),
  
                         fluidRow(
                           column(12,htmlOutput('venn'))
                         )                
                ),

                # Tab for Gene Set Enrichment Analysis
                tabPanel("GeneSet Enrichment",
                         fluidRow(
                           column(6,   
                                  selectInput("tabGSEA", "Comparison:",gsub("\\)","",gsub("P.Value \\(","",grep("P.Value \\(",names(tt),value=TRUE)))),
                                  # Update selection options via server.r, as this one reads in the BroadSets...
                                  selectInput("tabGSEA_Category", "Geneset Category:","All"),
                                  # Select a particular gene or multiple genes
                                  textInput("geneSetByGene", "Gene (EntrezID, separate by ','):", value=""),
                                  selectInput("cGSEA", "Color by factor:",factorNames),
                                  selectInput("fGSEA", "Group by factor:",factorNames)
                           ),
                           column(width=6, ggvisOutput("genSetVolcanoPlot"))
                         ),
                         
                         hr(),
                         
                         fluidRow(
                           tags$head(tags$style(HTML( ".has-feedback .form-control { padding-right: 0px;}" ))),
                           column(12,
                                  DT::dataTableOutput("genSetTab")
                           ),
                           hr(),
                           column(12,
                                  DT::dataTableOutput("genSetTab_genes")
                           )
                         ),
                         h3("Save pictures:"),
                         fluidRow(
                           column(
                             width = 12,
                             column(width = 2, numericInput("fheight_genset", "Height (cm)", min=2, max=15, step=1, value = 10)),
                             column(width = 2, numericInput("fwidth_genset", "Width (cm)", min=2, max=15, step=1, value = 10)),
                             column(width = 2, selectInput("fres_genset", "Res", choices=c("100","200","300"), selected = "100")),
                             column(width = 2, selectInput("fformat_genset", "File Type", choices=c("png","tiff","jpeg","pdf"), 
                                                    selected = "png", multiple = FALSE, selectize = TRUE))
                           ),
                           tags$head(tags$style(".butt{background-color:#337ab7;} .butt{color: #fff;}")),
                           column( 3, offset=3,
                             downloadButton('gs_download', 'Download Plots', class = "butt")
                           )
                           # needed for proper page layout
                           #    cellArgs = list(style="width: 400px;")    
                         ),
                         hr(),
  
                         absolutePanel(
                           top = 60, right = 20, width = 500, draggable = TRUE,
                           wellPanel(
                             HTML(markdownToHTML(fragment.only=TRUE, text=c(  "Beeswarm plot of the selected Gene (`draggable`)" ))),
                             plotOutput("beeswarm_plot2", height="400px")
                           ),
                           style = "opacity: 0.92"
                         ) #,
                         # absolutePanel(
                         #   top = 20, right = 60, width = 300, draggable = TRUE,
                         #   wellPanel(
                         #     HTML(markdownToHTML(fragment.only=TRUE, text=c("Volcano plot of the selected GeneSet (`draggable`)" ))),
                         #     plotOutput("genSetVolcanoPlot", height="200px")
                         #   ),
                         #   style = "opacity: 0.92"
                         # )
                ),
                
                # Tab for Venn diagrams of Gene Sets
                tabPanel("Venn diagrams - Gene Sets",
                         fluidRow(
                           column(2, numericInput("pVal_venn_gs", "Cut off p-value:",0.05, min=0.0, max=1.0)),
                           column(4, selectInput("fdr_venn_gs", "Multiple testing correction",c("none","Benjamini-Hochberg")))
                         ),
                         fluidRow(
                           column(width = 7,
                           # https://stackoverflow.com/questions/44784083/shiny-label-position-textinput
                             tags$form(
                               class="form-horizontal",
                               tags$div(
                                 class="form-group",
                                 tags$label(class = "col-sm-4 control-label", `for` = "colUp_gs", br(), "Color for genes:    Up"),
                                 column(width = 2, colourInput("colUp_gs", "", "red", palette="limited", showColour="background")),
                                 tags$label(class = "col-sm-2 control-label", `for` = "colDown_gs", br(), "Down"),
                                 column(width = 2, colourInput("colDown_gs", "", "blue", palette="limited", showColour="background"))
                               )
                             )
                           )
                         ),
                         fluidRow(actionButton("calculate_gs", "Calculate/Update",icon("sync"), 
                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                         ),

                         hr(),

                         fluidRow(
                           column(12,   
                             selectInput("vennX_gs", "Comparison:",gsub("venn_|.html","",
                                                    list.files(path="14_libra_deg/venn_gs/", pattern="all")), width='600px')
                             )
                         ),

                         hr(),

                         fluidRow(
                           column(12,htmlOutput('venn_gs'))
                         )                
                ),

                # Generate Beeswarm and Probe expression barplots for selected Genes
                tabPanel("Selected Genes",
                         fluidRow(
                           column(6,   
                                  textInput("s", "Gene (SYMBOL, EntrezID or Ensembl Gene ID, separate by ','):", value="")
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput("cSel", "Color by factor:",factorNames)
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput("fSel", "Group by factor:",factorNames)
                           )
                         ),
                         fluidRow(
                           column(6,   
                                  textInput("dir", "Files will be saved into <dir>.zip:",value=""),
                                  downloadButton('downloadPictures', 'Download')                          
                           )
                         )
                )   
    )
  )
  
))
  
