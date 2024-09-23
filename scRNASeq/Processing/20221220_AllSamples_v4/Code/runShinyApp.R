
    # Script to start Shiny app
    # - Detects whether we are on the CDW
    # - Creates a local directory for the R libraries to be installed
    # - Makes sure to use the Chrome browser if detected
    
    # https://www.r-bloggers.com/deploying-desktop-apps-with-r/
    dir.create("RLibs")
    .libPaths("RLibs")
    
    # Check availability of BiocManager package
    if (!requireNamespace("BiocManager") && tools:::.BioC_version_associated_with_R_version() >= 3.8 ){
      install.packages("BiocManager")
    }
    
    #
    if (!require("shiny")){ 
      if ( version$major >= 4 || version$minor > 5.1 ){
        BiocManager::install("shiny")
      } else {
        source("http://bioconductor.org/biocLite.R")
        biocLite("shiny", suppressUpdates=TRUE) 
      }
    }
    
    # Check where we are, i.e. local or CDW
    checkLocation <- system("WHERE /F /R C:\\ DirectFlex.exe", intern=TRUE)
    if ( is.null(attributes(checkLocation)) ){
      # We are on the CDW...
      myBrowser <- checkLocation
    } else {
      #Set default browser to Chrome
      myBrowser <- system("WHERE /F /R C:\\ chrome.exe", intern=TRUE)
    }
    
    options(browser = myBrowser)
    
    # You can probably also use IE, but maybe not all figure/pages will work - CHECK!
    if ( !is.null(attributes(myBrowser)) ){
      # Set default browser to IE
      myBrowser <- system("WHERE /F /R C:\\ iexplore.exe", intern=TRUE)
      # # Search returns multiple hits, take the one in "Program File" and "Internet Explorer" directory
      myBrowser <- myBrowser[grep("Program Files\\",myBrowser)]
      options(browser = myBrowser[1])
    }
    
    message('library paths:
', paste('... ', .libPaths(), sep='', collapse='
'))
    
    launch.browser = function(appUrl, browser.path=myBrowser) {
      message('Browser path: ', browser.path)
      #  shell(sprintf('"%s %s"', browser.path, appUrl), intern=TRUE)
      system("cmd.exe", input=paste(browser.path, appUrl, sep=" "), intern=TRUE)
    }
    
    shiny::runApp('./', launch.browser=launch.browser)
    
    # Perhaps remove temporary directory that holds the packages?
    # unlink("RLibs")
    
  
