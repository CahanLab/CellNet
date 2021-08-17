#Adaptation to ccn_hmClass to prevent re-ordering the groups by name
acn_queryClassHm <- function (classMat, grps = NULL, isBig = FALSE, cRow = FALSE, 
                         cCol = FALSE, fontsize_row = 4, fontsize_col = 4, main = NA, 
                         scale = FALSE, customAnnoColor = NULL, ...) 
{
   cools <- colorRampPalette(c("black", "limegreen", "yellow"))(100)
   if (isBig) {
      bcol <- NA
   } else {
      bcol <- "grey60"
   }
   if (is.null(grps)) {
      if (scale) {
         mymin <- min(classMat)
         mymax <- max(classMat)
      }
      else {
         mymin <- 0
         mymax <- 1
      }
      return(pheatmap::pheatmap(classMat, col = cools, breaks = seq(from = mymin, 
                                                                    to = mymax, length.out = 100), cluster_rows = cRow, 
                                cluster_cols = cCol, main = main, fontsize_row = fontsize_row, border_color=bcol,
                                fontsize_col = fontsize_col, ...))
   }
   else {
      #grps <- grps[order(grps)]
      cells <- names(grps)
      #groupNames <- sort(unique(grps))
      groupNames <- unique(grps)
      xcol <- colorRampPalette(rev(brewer.pal(n = 12, name = "Set3")))(length(groupNames))
      names(xcol) <- groupNames
      anno_colors <- list(group = xcol)
      xx <- data.frame(group = as.factor(grps))
      rownames(xx) <- cells
      if (scale) {
         mymin <- min(classMat)
         mymax <- max(classMat)
      }
      else {
         mymin <- 0
         mymax <- 1
      }
      if (is.null(customAnnoColor) == FALSE) {
         if (!all(groupNames %in% names(customAnnoColor))) {
            stop(paste0("Not all group name has a color in custom color vetor.", 
                        "\n"))
         }
         customAnnoColor <- customAnnoColor[groupNames]
         names(customAnnoColor) <- groupNames
         anno_colors <- list(group = customAnnoColor)
         xx <- data.frame(group = as.character(grps))
         rownames(xx) <- cells
      }
      return(pheatmap::pheatmap(classMat, col = cools, breaks = seq(from = mymin, 
                                                                    to = mymax, length.out = 100), cluster_rows = cRow, 
                                cluster_cols = cCol, show_colnames = FALSE, annotation_names_row = FALSE, 
                                annotation_col = xx, main = main, annotation_names_col = FALSE, 
                                annotation_colors = anno_colors, fontsize_row = fontsize_row, border_color=bcol,
                                fontsize_col = fontsize_col, ...))
   }
}


heatmapPlotly <- function(c_scoresMatrix, sampTab) {
   melted_tab = data.frame("classificationScore" = NULL, 
                           "sampleName" = NULL, 
                           "tissueType" = NULL)
   
   for(sampleName in colnames(c_scoresMatrix)) {
      temp_cscore = c_scoresMatrix[, sampleName]
      samp_descrip <- sampTab[sampleName, "description1"]
      tempTab = data.frame("classificationScore" = temp_cscore, 
                           "sampleName" = sampleName,
                           "description" = samp_descrip,
                           "tissueType" = rownames(c_scoresMatrix))
      melted_tab = rbind(melted_tab, tempTab)
   }
   melted_tab$description <- as.character(melted_tab$description)
   melted_tab[which(is.na(melted_tab$description)),"description"] <- "random"
   melted_tab$description <- as.character(lapply(melted_tab$description, 
                                          FUN = function(descrip){
                                             return(paste(strwrap(descrip, width=15), collapse="\n"))
                                          }))
   
   melted_tab$tissueType <- factor(melted_tab$tissueType,
                                  levels = sort(unique(melted_tab$tissueType), decreasing=TRUE))
   cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
   
   p <- ggplot(data = melted_tab, aes(sampleName, tissueType, fill = classificationScore))+
      geom_tile(color = "white") +
      facet_grid(.~description, scales="free_x") +
      scale_fill_gradient2(low = cools[1], high = cools[length(cools)], mid = cools[length(cools)/2], 
                           midpoint = 0.5, limit = c(0,1), space = "Lab", 
                           name="Classification Score") +
      xlab("Samples")+
      ylab("Tissue Types")+
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 8, hjust = 1)) +
      theme(axis.title.x = element_blank()) +
      theme(panel.spacing = unit(0, "lines"))
   
   return(ggplotly(p) %>% layout(margin=c(0,0,0.05,0.05)))
}


heatmapPlotlyRef <- function(c_scoresMatrix, sampTab) {
   melted_tab = data.frame("classificationScore" = NULL, 
                           "sampleName" = NULL, 
                           "citation" = NULL,
                           "tissueType" = NULL,
                           "description" = NULL)
   
   for(sampleName in colnames(c_scoresMatrix)) {
      temp_cscore = c_scoresMatrix[, sampleName]
      samp_citation <- sampTab[sampleName, "citation"]
      samp_descrip <- sampTab[sampleName, "description1"]
      tempTab = data.frame("classificationScore" = temp_cscore, 
                           "sampleName" = sampleName,
                           "citation" = samp_citation,
                           "tissueType" = rownames(c_scoresMatrix),
                           "description" = samp_descrip
                           )
      melted_tab = rbind(melted_tab, tempTab)
   }
   melted_tab$citation <- as.character(melted_tab$citation)
   melted_tab[which(is.na(melted_tab$citation)),"citation"] <- "random"
   melted_tab$citation <- as.character(lapply(melted_tab$citation, 
                                                 FUN = function(descrip){
                                                    return(paste(strwrap(descrip, width=15), collapse="\n"))
                                                 }))
   melted_tab$description <- as.character(melted_tab$description)
   melted_tab[which(is.na(melted_tab$description)),"description"] <- "random"
   melted_tab$tissueType <- factor(melted_tab$tissueType,
                                   levels = sort(unique(melted_tab$tissueType), decreasing=TRUE))
   cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
   
   p <- ggplot(data = melted_tab, aes(sampleName, tissueType, fill = classificationScore, text=description))+
      geom_tile(color = "white") + 
      #facet_grid(.~citation, scales="free_x") +
      facet_wrap(~citation, scales="free_x") + 
      scale_fill_gradient2(low = cools[1], high = cools[length(cools)], mid = cools[length(cools)/2], 
                           midpoint = 0.5, limit = c(0,1), space = "Lab", 
                           name="Classification Score") +
      xlab("Samples")+
      ylab("Tissue Types")+
      theme_minimal()+ 
      theme(axis.text.x=element_blank()) +
      theme(panel.spacing = unit(0, "lines"))
   
   return(ggplotly(p, tooltip="all") %>% layout(margin=c(0,0,0.05,0.02)))
}
