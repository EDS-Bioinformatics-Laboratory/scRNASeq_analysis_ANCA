# Script excerpt of the analysis run for Pete Quinn
# in order to generate clickable Venn Diagrams

# Function to generate the list of HTML pages to be used in creating the clickable Venn diagram
makeHtmlList <- function (vennMatrix,htmlFiles, reportDir="venn")
{
  # Make lists of the individual HTML files for the Up and Down
  # Just one contrast
  if(ncol(vennMatrix) == 1){
    myVennUp = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up"))  ) )

    myVennDown = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down"))  ) )
  }
  # - for two contrasts
  if(ncol(vennMatrix) == 2){
    myVennUp = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up"))
    ) )

    myVennDown = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down"))
    ) )
  }
  # For three comparisons
  if(ncol(vennMatrix) == 3){
    myVennUp = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("up_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up"))
    ) )

    myVennDown = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("down_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down"))
    ) )
  }
  # And for four contrasts
  if(ncol(vennMatrix) == 4){
    myVennUp = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("up_x_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up"))
    ) )

    myVennDown = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("down_x_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_down_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_x_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_x_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_down_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_x_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_x_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_down_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_down_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_x_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_down_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_down_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_x_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("x_down_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                    HTMLReport(shortName=htmlFiles[grep("down_down_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down"))
    ) )
  }

  # And for five contrasts
  if(ncol(vennMatrix) == 5){
    myVennUp = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("up_x_x_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_x_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_up_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_x_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_x_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),								  
                                  HTMLReport(shortName=htmlFiles[grep("up_up_x_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_up_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_x_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_x_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_up_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_x_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_x_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_up_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_up_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_x_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_up_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),								  
                                  HTMLReport(shortName=htmlFiles[grep("up_up_x_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_x_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_up_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_x_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_up_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_up_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_up_up_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_up_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_up_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_x_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_up_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_up_up_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_up_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_up_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("x_up_x_up_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up")),
                                  HTMLReport(shortName=htmlFiles[grep("up_x_up_x_up",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_up"))
    ) )
    myVennDown = list(vennout=list(HTMLReport(shortName=htmlFiles[grep("down_x_x_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_down_x_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_down_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_x_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_x_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),								  
                                  HTMLReport(shortName=htmlFiles[grep("down_down_x_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_x_down_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_x_x_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_x_x_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_down_down_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_down_x_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_down_x_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_down_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_down_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_x_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_down_down_x_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),								  
                                  HTMLReport(shortName=htmlFiles[grep("down_down_x_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_down_x_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_x_down_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_x_x_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_down_down_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_x_down_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_down_down_down_x",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_down_down_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_x_down_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_down_x_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_down_down_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_down_down_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_down_down_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_x_down_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("x_down_x_down_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down")),
                                  HTMLReport(shortName=htmlFiles[grep("down_x_down_x_down",htmlFiles)],reportDirectory=paste0("./",reportDir,"/venn_down"))
    ) )
  }
  htmlList = list(up=myVennUp, down=myVennDown)
  return(htmlList)
}

# Function to generate the coordinates and write the HTML pages
#
# You can obtain image coordinates from fx. http://imagemap-generator.dariodomi.de/
#
makeVennHTML <- function (lst, page, mapname, image, reportDirectory, originalDir, condition) 
{
  fun <- function(x, y) paste0("<area shape=\"rect\" coords=", x, " href=", y, "></area>")
  if(length(lst$vennout)==1){
    if(condition=="up"){
	    loclst <- list("370,375,425,405")
	} else {
	    loclst <- list("370,405,425,430")
	} 
  }
  if(length(lst$vennout)==3){
    if(condition=="up"){
        loclst <- list("200,350,300,390","500,350,600,390","350,350,450,390")
	} else {
        loclst <- list("200,395,300,435","500,395,600,435","350,395,450,435")
    }    
  }
  if(length(lst$vennout)==7){
    if(condition=="up"){
        loclst <- list("230,280,290,310","510,280,570,310","370,520,430,550",
                       "370,290,430,320","300,390,360,430","440,390,500,430","370,370,430,400")
	} else {
        loclst <- list("230,310,290,340","510,310,570,340","370,550,430,580",
                       "370,320,430,350","300,430,360,460","440,430,500,460","370,400,430,430" )
    }    
  }
  if(length(lst$vennout)==15){
    if(condition=="up"){
        loclst <- list("110,285,170,320", "285,185,345,215", "480,185,540,215", 
                       "655,295,725,325", "200,240,260,270", "200,510,260,550", "385,615,445,645", 
                       "385,225,445,265", "545,510,605,550", "565,240,625,280", "270,315,330,355", 
                       "455,555,515,595", "310,555,370,595", "500,315,560,355", "385,440,445,480")
	} else {
        loclst <- list("110,325,170,360", "285,220,345,255", "480,220,540,255", 
                       "655,330,725,365", "200,275,260,310", "200,555,260,590", "385,650,445,685", 
                       "385,270,445,305", "545,555,605,590", "565,285,625,320", "270,360,330,395", 
                       "455,600,515,635", "310,600,370,635", "500,360,560,395", "385,485,445,520")
    }
  }
  if(length(lst$vennout)==32){
    if(condition=="up"){
        loclst <- list("165,340,190,360","380,165,410,190","580,315,610,335","530,580,560,600",
		               "275,600,305,620","300,305,330,325","265,405,295,425","515,465,545,485",
					   "290,475,320,495","455,275,485,295","350,275,380,295","370,520,400,540",
					   "525,400,555,420","495,330,525,350","425,540,455,560","300,370,330,390",
					   "320,300,350,320","345,480,375,500","525,430,540,440","475,475,505,495",
					   "395,300,425,320","490,370,520,390","340,335,370,355","450,335,480,355",
					   "495,415,515,435","410,480,440,500","315,435,345,455","395,405,425,425",
					   "480,310,510,330","525,425,555,445","405,520,435,540","290,445,310,465")
	} else {
        loclst <- list("165,360,190,380","380,195,410,220","580,335,610,355","530,600,560,620",
		               "275,620,305,640","290,325,320,345","265,425,295,445","515,485,545,505",
					   "290,495,320,515","455,295,485,315","350,295,400,315","370,540,400,560",
					   "525,420,555,440","495,350,525,370","425,560,455,580","300,390,330,410",
					   "320,310,350,330","345,500,375,520","525,440,540,455","475,495,505,515",
					   "395,320,425,340","490,390,520,410","340,355,370,375","450,355,480,375",
					   "495,440,515,455","410,500,440,520","315,455,345,475","395,425,425,445",
					   "480,330,510,355","525,445,555,465","405,540,435,560","300,465,330,485")
    }
  }
  nam2 <- paste0(image, ".png")
  #Original code: urlst <- gsub(reportDirectory, ".", sapply(lst$vennout, path))
  urlst <- gsub(originalDir, ".", sapply(lst$vennout, path))
  strng <- do.call("c", mapply(fun, loclst, urlst, SIMPLIFY = FALSE))
  strng <- c(paste0("<map name=\"", sub("#", "", mapname),  "\">"), strng, "</map>")
  
  # print(getwd())
  # print(paste0("Opening ",paste0(reportDirectory,"/",page, ".html")))
  
  hpage <- openPage(paste0(page, ".html"), dirname = reportDirectory)
  hwrite(paste("The Venn diagrams all contain clickable links. Click on the counts", 
               "to see a table of the genes involved.", 
               "Also please note that the tables are sortable - simply click", 
               "on any header to sort on that column."), hpage, br = TRUE)
  
  hwrite(hmakeTag("img", border = 0, width = 800, height = 800, 
                  src = nam2, alt = nam2, usemap = mapname), hpage)
  #hwriteImage(paste("Venn Diagram", num), hpage)
  cat("\n",file = hpage)
  cat(strng, file = hpage, sep = "\n")
  
  closePage(hpage)
  paste0(reportDirectory, "/", page, ".html")
}

# Here I should also have parameters determining the column in the 'tt' object that holds the 'entrezgene' information,
# as that one gets an HTML link. The other annotation columns should be made more generic, they should be able to keep
# the heading they have in the original 'tt' object
makeClickableVenn <- function(vennMatrix, tt, indxEntrez, indxDesc, indxEnsembl, 
                              reportDir="venn", pValString="p.value = 0.01 (not adjusted)",
                              myCol=c("red","blue"))  
  # In order to shorthen the name of the contrasts, I tried to replace the name with the index of the comparison:
  # see fx. vennReportsGeneSets.r and vennReports.r where we have the indices variable, but I forgot to add this to the
  # final call to makeClickableVenn.... 
  # Now, I have added the indices to the colnames of the testMatrix and perform a "string split" to get the indices in the
  # final name of the comparisons.
  # Another option would be to add the indices variable to the call to makeClickableVenn and then add ",indices=1" in the
  # function declaration above.
  

{
  if (!require("pacman")){ 
    if ( paste0(R.Version()[c("major","minor")], collapse = ".") >= "3.6.0" ){
      BiocManager::install("pacman", dependencies = TRUE)
      
      # install packages using 'pacman'
      pacman::p_load(hwriter, cwhmisc, stringr, update = FALSE)
      
    } else {
      if ( paste0(R.Version()[c("major","minor")], collapse = ".") >= "3.5.0" ){
         source("http://bioconductor.org/biocLite.R")
         biocLite("pacman", suppressUpdates=TRUE)
         
         # install packages using 'pacman'
         pacman::p_load(hwriter, cwhmisc, update = FALSE)
         
      } else {
        source("http://bioconductor.org/biocLite.R")
        if ( require("hwriter") == FALSE ) { biocLite("hwriter", suppressUpdates=TRUE) }
        if ( require("cwhmisc") == FALSE ) { biocLite("cwhmisc", suppressUpdates=TRUE) }
      }
    }
  }
  
  # Biobase, limma and edgeR are packages from Bioconductor
  if(require("limma") == FALSE) {
    if ( paste0(R.Version()[c("major","minor")], collapse = ".") >= "3.6.0" ){
      BiocManager::install("limma", dependencies = TRUE)
    } else {
      source("http://bioconductor.org/biocLite.R")
      biocLite("limma", suppressUpdates=TRUE)
    }
  }
  if(require("affycoretools") == FALSE) {
    if ( paste0(R.Version()[c("major","minor")], collapse = ".") >= "3.6.0" ){
      BiocManager::install("affycoretools", dependencies = TRUE)
    } else {
      source("http://bioconductor.org/biocLite.R")
      biocLite("affycoretools", suppressUpdates=TRUE)
    }
  }
  if(require("ReportingTools") == FALSE) {
    if ( paste0(R.Version()[c("major","minor")], collapse = ".") >= "3.6.0" ){
      BiocManager::install("ReportingTools", dependencies = TRUE)
    } else {
      source("http://bioconductor.org/biocLite.R")
      biocLite("ReportingTools", suppressUpdates=TRUE)
    }
  }

  require(affycoretools)
  require(cwhmisc)
  require(limma)
  require(hwriter)
  require(ReportingTools)
  require(stringr)
  
#  if(!exists("indxEntrez")){
#    indxEntrez = grep("^entrezgene$",colnames(tt))
#  }
#  if(!exists("indxDesc")){
#    indxDesc = grep("^description",colnames(tt))
#  }
#  if(!exists("indxEnsembl")){
#    indxEnsembl = grep("^ensembl_gene_id",colnames(tt))
#  }
  
  # split the testMatrix into a matrix for "up" and "down" and convert to binair
  # For UP convert all 'down' interactions to 'no' (i.e zero)
  testMatrix.up = vennMatrix
  testMatrix.up[testMatrix.up == -1] = 0

  # And for DOWN, convert all 'up' interaction into zero and then all 'down' interactions to '1'
  testMatrix.down = vennMatrix
  testMatrix.down[testMatrix.down == 1] = 0
  testMatrix.down[testMatrix.down == -1] = 1

  # It seems as if vennCounts2 can not handle 4 contrasts...so, resort to Limma vennCount
  vennCountsResults_up = vennCounts(testMatrix.up)
  vennCountsResults_down = vennCounts(testMatrix.down)

  # Put the two count objects into a list, so you can loop over them
  matrices <- list(up=testMatrix.up, down=testMatrix.down)

  # Initiate a list to hold all generated HTML files
  htmlFiles = list()

  # Get rid of redundant information with the topTags/topTable, tt
  # This could be made more generic..but later on you again need these columns..
  for (col in 1:ncol(tt)){
    if ( class(tt[,col]) == "character" | class(tt[,col]) == "factor"){
      tt[,col] = as.character(sapply(tt[,col],function(x) ifelse(x != "",unique(unlist(strsplit(as.character(x),"//"))), NA)))
    }
  }
#  tt[,indxEntrez] = as.character(sapply(tt[,indxEntrez],function(x){if (x != ""){unique(unlist(strsplit(as.character(x),"//")))} else {NA}}))
#  tt[,indxDesc] = as.character(sapply(tt[,indxDesc],function(x){if (x != ""){unique(unlist(strsplit(as.character(x),"//")))} else {NA}}))
#  tt[,indxEnsembl] = as.character(sapply(tt[,indxEnsembl],function(x){if (x != ""){unique(unlist(strsplit(as.character(x),"//")))} else {NA}}))

  baseDir = getwd()
  
  for(i in 1:length(matrices)){
  
    testMatrix <- matrices[[i]]
    cat("\n",names(matrices)[i],"\n")

    possibleValues=c()
	
	  for(col in ncol(testMatrix):1){
      possibleValues = append(possibleValues, 2^(col-1))
    }
  
    allPossibleValues = unique(testMatrix %*% possibleValues)
  
    direction = c( "x",names(matrices)[i])
  
    for(val in allPossibleValues){
      cat(val, "\n")
    
      # Select all row (i.e. genes) that sum up to the desired value. Contains all combinations
      # of 1 and 0 that sum up to 'val'
      indx = which(testMatrix %*% possibleValues == val); 
      cat("Total counts:",length(indx), "\n")
    
      tt.subset = tt[indx,]
    
      contrastNameLong=paste(str_split(colnames(testMatrix),"_",n=2, simplify=TRUE)[,2], collapse="_v_")
      # The contrastNames cab become very long ...
      # Replace by index of the comparison
      # contrastName=paste0("comp_", paste(str_split(colnames(testMatrix),"_",n=2, simplify=TRUE)[,1], collapse="_"),"_",paste(indices, collapse="v"))
      contrastName=paste0("comp_", paste(str_split(colnames(testMatrix),"_",n=2, simplify=TRUE)[,1], collapse="_"))
      
      if (package.version("cwhmisc") < 6.1){
        signs = paste(direction[1 + as.numeric(strsplit2(formatC(as.numeric(intToBase(as.numeric(val),base=2),""), width = ncol(testMatrix), format = "d", flag = "0"), split=""))], collapse="_")
      } else {
        signs = paste(direction[1 + as.numeric(strsplit2(formatC(as.numeric(int2B(as.numeric(val),B=2)[[1]],""), width = ncol(testMatrix), format = "d", flag = "0"), split=""))], collapse="_")
      }
      contrastNameLong = paste(contrastNameLong, signs, sep="_")
      contrastName = paste(contrastName, signs, sep="_")
      
      
      htmlFiles = append(htmlFiles,paste0(contrastName,".html"))
    
      cat("\t",contrastNameLong,"\n")
      #print(paste0("./",reportDir,"/venn_",names(matrices)[i],"/"))
      dir.create(paste0("./",reportDir,"/venn_",names(matrices)[i],"/"), recursive=TRUE, showWarnings = FALSE)
      setwd(paste0("./",reportDir,"/venn_",names(matrices)[i],"/"))
    
      vennReport = HTMLReport(shortName=contrastName,title=paste0("venn (contrasts: ",contrastNameLong,")"),reportDirectory=".")
      csvReport = CSVFile(shortName=contrastName,reportDirectory=".")

      publish(hwrite("<FORM><INPUT Type=\"button\" VALUE=\"Back\" onClick=\"history.go(-1);return true;\"></FORM>"),vennReport)
	  
      addEGIDLink <- function(object, ...){        
        indxEntrez = grep("^entrezgene$",colnames(tt))
        indxEnsembl = grep("^ensembl_gene_id$",colnames(tt))
        indxGeneSet = grep("^GeneSet$",colnames(tt))
                
        for ( colName in colnames(object[1:ncol(object),-c(indxEnsembl,indxEntrez)]) ){
            colClass <- class(object[,colName])
            object[,colName] <- ifelse(!is.na(as.character(object[,colName])),
                                       hwrite(as.character(object[,colName]), table = FALSE), 
                                       hwrite("",table=FALSE)
            )
            if (colClass == "numeric"){ object[,colName] <- as.numeric(object[,colName]) }
        }
        
        if (exists("indxEntrez")){
          object[,indxEntrez] <- ifelse(!is.na(as.character(object[,indxEntrez])),
                                     hwrite(as.character(object[,indxEntrez]),link = paste0("http://www.ncbi.nlm.nih.gov/gene/",as.character(object[,indxEntrez])), table = FALSE),
                                     hwrite("NA",table=FALSE)
          )
        }

        if (exists("indxEnsembl")){  
          object[,indxEnsembl] <- ifelse(!is.na(as.character(object[,indxEnsembl])),
                                       hwrite(as.character(object[,indxEnsembl]),link = paste0("http://www.ensembl.org/Gene/Summary?g=",as.character(object[,indxEnsembl])), table = FALSE), 
                                       hwrite("NA",table=FALSE)
          ) 
        }
         
        if (exists("indxGeneSet")){  
          object[,indxGeneSet] <- ifelse(!is.na(as.character(object[,indxGeneSet])),
                                         hwrite(as.character(object[,indxGeneSet]),
                                                link = paste0("http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=",as.character(gsub("^C[1-7]_|^H_","",as.character(object[,indxGeneSet])))), table = FALSE), 
                                         hwrite("NA",table=FALSE)
          ) 
        }
        return(object)
      }
    
      # Flatten some list for CSV output
      csv.subset = as.character(tt.subset)
      #csv.subset$entrezgene = as.character(csv.subset$entrezgene)
      #csv.subset$ensembl_gene_id = as.character(csv.subset$ensembl_gene_id)
    
      # 'publish' won't work if there are no rows..
      # Do not output if not selected in the Venn Diagram (i.e. 'val' == 0)
      if(dim(tt.subset)[1] > 0 & val > 0){ 
        publish(tt.subset,vennReport,.modifyDF=list(addEGIDLink))
        #        publish(as.data.frame(csv.subset),csvReport)
        my.df <- data.frame(lapply(tt.subset, as.character), stringsAsFactors=FALSE)
        colnames(my.df) <- colnames(tt.subset)
        # replace all "," with "//", because we're going to export a comma separated file and Excel will choke on extra comma's...
        #        my.df <- data.frame(apply(apply(my.df,2,gsub, patt=",",replace="//"),2,as.character))
        my.df <- lapply(my.df, function(x){gsub(",","//",x)})
        
        # print("Writing CSV file ...")
        write.table(my.df,file=paste0("./", contrastName,".csv"),sep=",",quote=FALSE,row.names=FALSE)
      }
      else{
        publish(hwrite("No overlap found for this comparison"), vennReport)
      }
      # And add a Download button
	    # http://stackoverflow.com/questions/11620698/how-to-trigger-a-file-download-when-clicking-an-html-button-or-javascript
      publish(hwrite(paste0("<FORM method=\"get\" action=\"",paste0(contrastName,".csv"),"\"><button type=\"submit\">Download</button></FORM>")),vennReport)
      
      # Finish
      finish(vennReport)
      # print("HTML written ...")
      
      # And return to base ...
      setwd(baseDir)
    }
  }

  # Generate a base filename based on the contrasts being compared
  contrastName=paste0("venn_",paste(str_split(colnames(testMatrix),"_",n=2, simplify=TRUE)[,2], collapse="_and_"))

  # Recreate the Venn Diagram using Limma (with both Up and Down) to overwrite the Venn Diagram produced by affycoretools:::vennPage
  dir.create(paste0("./",reportDir), showWarnings = FALSE)
  png(paste0("./",reportDir,"/",contrastName,".png"),width = 800, height = 800)
    colnames(vennMatrix) <- str_split(colnames(vennMatrix),"_",n=2, simplify=TRUE)[,2]
  vennDiagram(vennMatrix, include=c("up","down"),counts.col=myCol, cex = 1.5 )
#  txtLine <- ifelse(ncol(testMatrix)==2, 1,ifelse(ncol(testMatrix)==3,3, 4))
#  atLine <- ifelse(ncol(testMatrix)==2, -2.0,ifelse(ncol(testMatrix)==3, -2.0, 50 ))
#  mtext(pValString,side=1,line=txtLine,cex=1.2, at=atLine)
  mtext(pValString,side=1,line=3,cex=1.2, adj=0)
  dev.off()

  # Generate the list of HTML files created above
  # print("Create list of HTML files")
  vennList <- makeHtmlList(testMatrix,htmlFiles, reportDir=reportDir)

  # Make HTML files for 'up' and 'down', respectively
  pwd <- getwd()  # Get the present working directory
  setwd(reportDir) # Go to the report directory, filenames get too long if you paste them together with the filepath!!
  #makeVennHTML(vennList$up, paste0(contrastName,".up"), "#Venn diagrams Up", contrastName, paste0("./",reportDir), "up")
  #makeVennHTML(vennList$down, paste0(contrastName,".down"), "#Venn diagrams Down", contrastName, paste0("./",reportDir), "down")
  makeVennHTML(vennList$up, paste0(contrastName,".up"), "#Venn diagrams Up", contrastName, reportDirectory = "./", originalDir = paste0("./",reportDir),  "up")
  makeVennHTML(vennList$down, paste0(contrastName,".down"), "#Venn diagrams Down", contrastName, reportDirectory = "./", originalDir = paste0("./",reportDir), "down")

  # Combine the individual HTML Files for Up and Down into one general HTML and take care that both Up and Down numbers can be clicked upon separately
  
  #conDown = file(paste0("./",reportDir,"/",contrastName,".down.html"))
  conDown = file(paste0("./", contrastName,".down.html"))
  htmlFileDown = readLines(conDown)
  htmlToInsert = grep("rect",htmlFileDown)
  close(conDown)

  #conUp = file(paste0("./",reportDir,"/",contrastName,".up.html"))
  conUp = file(paste0("./",contrastName,".up.html"))
  htmlFileUp = readLines(conUp)
  htmlFileAfter = grep("rect",htmlFileUp)
  close(conUp)

  # Combine Up and Down
  combinedHtml <- append(htmlFileUp, htmlFileDown[htmlToInsert], after=tail(htmlFileAfter,n=1))
  combinedHtml <- gsub("venn1",contrastName, combinedHtml)
#  combinedHtml <- gsub(paste0("./venn/",contrastName,".up.html"),contrastName, combinedHtml)
# 20201217 - AJ
# As I made changes to vennReports.r and vennGeneSets_Reports.r it now seems I have to outcomment this line below:
  #combinedHtml <- gsub(reportDir,"", combinedHtml)
  
  # Write into new Html
  #writeLines(combinedHtml,con=paste0("./",reportDir,"/",contrastName,"_all.html") )
  writeLines(combinedHtml,con=paste0("./",contrastName,"_all.html") )
  
  # Unload libraries
  #detach("package:affycoretools", unload=TRUE)
  #detach("package:ReportingTools", unload=TRUE)
  #detach("package:hwriter", unload=TRUE)
  #detach("package:cwhmisc", unload=TRUE)
  
  # and go back to the working directory
  setwd(pwd)
}
 
