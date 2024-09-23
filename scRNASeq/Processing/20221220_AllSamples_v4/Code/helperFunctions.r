#### my own publish_gosttable function ####
my_publish_gosttable <- function (gostres, highlight_terms = NULL, use_colors = TRUE, 
                                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                                  sortCol = NULL, filename = NULL)
{
  term_id <- p_values <- query <- p_value <- NULL
  if (class(gostres) == "list") {
    if (!("result" %in% names(gostres))) 
      stop("Name 'result' not found from the input")
    df <- gostres$result
  }
  else if (class(gostres) == "data.frame") {
    df <- gostres
  }
  else {
    stop("The input 'gostres' should be a data frame or a list from the gost() function.")
  }
  if (!"term_id" %in% colnames(df)) 
    stop("The required column 'term_id' is missing from the input")
  if (!any(grepl("p_value", colnames(df)))) 
    stop("Column 'p_value(s)' is missing from the input")
  if (is.null(highlight_terms)) {
    highlight_terms = df
  }
  if (is.data.frame(highlight_terms)) {
    message("The input 'highlight_terms' is a data.frame. The column 'term_id' will be used.")
    if ("term_id" %in% colnames(highlight_terms)) {
      highlight_terms <- highlight_terms$term_id
    }
    else {
      stop("No column named 'term_id'.")
    }
  }
  subdf <- base::subset(df, term_id %in% highlight_terms)
  if (nrow(subdf) == 0) {
    stop("None of the term IDs in the 'highlight_terms' were found from the results.")
  }
  subdf$id <- match(subdf$term_id, highlight_terms)
  show_columns <- unique(append(show_columns, c("id","term_id", "p_value")))
  gp_colnames <- c("id", "source", "term_id", 
                   "term_name", "term_size", "query_size", 
                   "intersection_size", "p_value")
  colnames <- gp_colnames[which(gp_colnames %in% show_columns)]
  if (length(setdiff(show_columns, gp_colnames)) > 0) {
    colnames <- append(colnames, setdiff(show_columns, gp_colnames))
  }
  if ("p_values" %in% colnames(subdf)) {
    if ("meta" %in% names(gostres)) {
      meta <- gostres$meta
      subdf$query <- list(names(meta$query_metadata$queries))
    }
    else {
      qnames = paste("query", seq(1, length(subdf$p_values[[1]])), sep = "_")
      subdf$query <- list(names(qnames))
    }
    subdf <- tidyr::unnest(data = subdf, p_values, query)
    subdf <- dplyr::rename(subdf, p_value = p_values)
    subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 3)
    showdf <- subdf[, stats::na.omit(match(c(colnames, "query"),  names(subdf)))]
    showdf <- tidyr::spread(showdf, query, p_value)
    idx <- which(!is.na(match(names(showdf), unique(subdf$query))))
  }
  else {
    if ("query" %in% names(subdf) & length(unique(subdf$query)) > 1) {
      subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 3)
      showdf <- subdf[, stats::na.omit(match(c(colnames, "query"), names(subdf)))]
      showdf <- tidyr::spread(showdf, query, p_value)
      idx <- which(!is.na(match(names(showdf), unique(subdf$query))))
    }
    else {
      subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 3)
      showdf <- subdf[, stats::na.omit(match(colnames, names(subdf)))]
      idx <- which(!is.na(match(names(showdf), "p_value")))
    }
  }
  if ( !is.null(sortCol) ){ 
    if ( sortCol == "p_value" ){ 
      showdf$p_value <- as.numeric(showdf$p_value)
    }
    showdf <- showdf[order(showdf[sortCol]),]
    message(paste0("The table will be sorted on ",sortCol, "."))
  }
  
  if (is.null(filename)) {
    if ( nrow(showdf) > 100 ){
      message(paste0("The table contains more than 100 rows and therefore hardly be legible. Please provide a filename"))
    } else {  
      colours <- matrix("white", nrow(showdf), ncol(showdf))
      if (use_colors) {
        cols <- sapply(showdf[, idx], function(x) mapViridis(-log10(as.numeric(x))))
      }
      else {
        cols <- sapply(showdf[, idx], function(x) "white")
      }
      colours[, idx] <- cols
      fontcolours <- matrix("black", nrow(showdf), ncol(showdf))
      fontcolours[, idx] <- ifelse(use_colors, "white", "black")
      fontfaces <- matrix("plain", nrow(showdf), ncol(showdf))
      fontfaces[, idx] <- "bold"
      showdf[is.na(showdf)] <- ""
      #hj <- matrix(c(1, 1, 0, 0, 1, 1, 1), ncol=7, nrow=nrow(showdf), byrow=TRUE)
      th <- gridExtra::ttheme_default(base_size = 10, 
                                      padding = grid::unit(c(4,  4), "mm"), 
                                      core = list(padding.h = grid::unit(c(15,15), "mm"), 
                                                  padding.v = grid::unit(c(15,15), "mm"), 
                                                  bg_params = list(fill = colours, col = "black", lwd = 0.5), 
                                                  fg_params = list(hjust = 0, x = 0.01, col = fontcolours, fontface = fontfaces)), 
                                      colhead = list(bg_params = list(fill = "gray99", lwd = 0.5, col = "black"), 
                                                     fg_params = list(col = "gray39",  fontface = "bold")), 
                                      rowhead = list(fg_params = list(col = "black",fontface = "bold")))
      tb <- gridExtra::tableGrob(showdf, theme = th, rows = NULL)
      h <- grid::unit.c(sum(tb$heights))
      w <- grid::unit.c(sum(tb$widths))
      tg <- gridExtra::grid.arrange(tb, ncol = 1, widths = w, heights = h, 
                                    newpage = TRUE, bottom = grid::textGrob("g:Profiler (biit.cs.ut.ee/gprofiler)", 
                                                                            x = 0.95, hjust = 1, gp = grid::gpar(fontsize = 10, 
                                                                                                                 font = 8, col = "cornflowerblue")))
      p <- ggplot2::ggplot() + 
        ggplot2::annotation_custom(tg) + 
        ggplot2::geom_blank() + 
        ggplot2::theme_void()
      return(p)
    }
  }
  else {
    imgtype <- strsplit(basename(filename), split = "\\.")[[1]][-1]
    if (length(imgtype) == 0) {
      filename = paste0(filename, ".pdf")
    }
    if (tolower(imgtype) %in% c("png", "pdf", "jpeg", "tiff", "bmp")) {
      if ( nrow(showdf) > 25 ){
        myDF <- split(showdf, (seq(nrow(showdf))-1) %/% 25)
        message(paste0("The table contains more than 25 entries (i.e. ", nrow(showdf),") and will be split in multiple (",length(myDF),") plots."))
        for ( i in 1:length(myDF) ){
          colours <- matrix("white", nrow(myDF[[i]]), ncol(myDF[[i]]))
          if (use_colors) {
            cols <- sapply(myDF[[i]][, idx], function(x) mapViridis(-log10(as.numeric(x))))
          }
          else {
            cols <- sapply(myDF[[i]][, idx], function(x) "white")
          }
          colours[, idx] <- cols
          fontcolours <- matrix("black", nrow(myDF[[i]]), ncol(myDF[[i]]))
          fontcolours[, idx] <- ifelse(use_colors, "white", "black")
          fontfaces <- matrix("plain", nrow(myDF[[i]]), ncol(myDF[[i]]))
          fontfaces[, idx] <- "bold"
          myDF[[i]][is.na(myDF[[i]])] <- ""
          
          #hj <- matrix(c(1, 1, 0, 0, 1, 1, 1), ncol=7, nrow=nrow(myDF[[i]]), byrow=TRUE)
          th <- gridExtra::ttheme_default(base_size = 10, 
                                          padding = grid::unit(c(4,  4), "mm"), 
                                          core = list(padding.h = grid::unit(c(15,15), "mm"), 
                                                      padding.v = grid::unit(c(15,15), "mm"), 
                                                      bg_params = list(fill = colours, col = "black", lwd = 0.5), 
                                                      fg_params = list(hjust = 0, x = 0.01, 
                                                                       col = fontcolours,
                                                                       fontface = fontfaces)), 
                                          colhead = list(bg_params = list(fill = "gray99", lwd = 0.5, col = "black"), 
                                                         fg_params = list(col = "gray39",  fontface = "bold")), 
                                          rowhead = list(fg_params = list(col = "black",fontface = "bold")))
          
          tb <- gridExtra::tableGrob(myDF[[i]], theme = th, rows = NULL)
          h <- grid::unit.c(sum(tb$heights))
          w <- grid::unit.c(sum(tb$widths))
          tg <- gridExtra::grid.arrange(tb,ncol = 1, widths = w, heights = h, 
                                        newpage = TRUE, 
                                        bottom = grid::textGrob("g:Profiler (biit.cs.ut.ee/gprofiler)", 
                                                                x = 0.95, hjust = 1, 
                                                                gp = grid::gpar(fontsize = 10,font = 8, col = "cornflowerblue"))
          )
          p <- ggplot2::ggplot() + 
            ggplot2::annotation_custom(tb) + 
            ggplot2::geom_blank() + 
            ggplot2::theme_void()
          width = grid::convertWidth(sum(tg$widths), "in", 
                                     TRUE) + 0.5
          height = grid::convertHeight(sum(tg$heights), "in", 
                                       TRUE) + 0.5
          filename <- strsplit(filename, split = "\\.")[[1]][1]
          filename_p <- paste0(filename, "_p",i,".",tolower(imgtype))
          ggplot2::ggsave(filename = filename_p, plot = p + theme(plot.margin = unit(c(0,2,0,2), "cm")), height = height, width = width)
          message("The image is saved to ", filename_p)
        } 
      } 
      else {
        colours <- matrix("white", nrow(showdf), ncol(showdf))
        
        if (use_colors) {
          cols <- sapply(showdf[, idx], function(x) mapViridis(-log10(as.numeric(x))))
        }
        else {
          cols <- sapply(showdf[, idx], function(x) "white")
        }
        
        colours[, idx] <- cols
        fontcolours <- matrix("black", nrow(showdf), ncol(showdf))
        fontcolours[, idx] <- ifelse(use_colors, "white", "black")
        fontfaces <- matrix("plain", nrow(showdf), ncol(showdf))
        fontfaces[, idx] <- "bold"
        showdf[is.na(showdf)] <- ""
        #hj <- matrix(c(1, 1, 0, 0, 1, 1, 1), ncol=7, nrow=nrow(showdf), byrow=TRUE)
        th <- gridExtra::ttheme_default(base_size = 10, 
                                        padding = grid::unit(c(4,  4), "mm"), 
                                        core = list(padding.h = grid::unit(c(15,15), "mm"), 
                                                    padding.v = grid::unit(c(15,15), "mm"), 
                                                    bg_params = list(fill = colours, col = "black", lwd = 0.5), 
                                                    fg_params = list(hjust = 0, x = 0.01, col = fontcolours, fontface = fontfaces)), 
                                        colhead = list(bg_params = list(fill = "gray99", lwd = 0.5, col = "black"), 
                                                       fg_params = list(col = "gray39",  fontface = "bold")), 
                                        rowhead = list(fg_params = list(col = "black",fontface = "bold")))
        tb <- gridExtra::tableGrob(showdf, theme = th, rows = NULL)
        h <- grid::unit.c(sum(tb$heights))
        w <- grid::unit.c(sum(tb$widths))
        tg <- gridExtra::grid.arrange(tb, ncol = 1, widths = w, heights = h, 
                                      newpage = TRUE, bottom = grid::textGrob("g:Profiler (biit.cs.ut.ee/gprofiler)", 
                                                                              x = 0.95, hjust = 1, gp = grid::gpar(fontsize = 10, 
                                                                                                                   font = 8, col = "cornflowerblue")))
        p <- ggplot2::ggplot() + 
          ggplot2::annotation_custom(tg) + 
          ggplot2::geom_blank() + 
          ggplot2::theme_void()
        
        width = grid::convertWidth(sum(tg$widths), "in", TRUE) + 0.5
        height = grid::convertHeight(sum(tg$heights), "in", TRUE) + 0.5
        ggplot2::ggsave(filename = filename, plot = p + theme(plot.margin = unit(c(0,2,0,2), "cm")), height = height, 
                        width = width)
        message("The image is saved to ", filename)
        return(p)
      }
    } 
    else {
      stop("The given file format is not supported.\nPlease use one of the following extensions: .png, .pdf, .jpeg, .tiff, .bmp")
    }
  }
}