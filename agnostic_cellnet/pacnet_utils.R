#Adaptation to ccn_hmClass to prevent re-ordering the groups by name
acn_queryClassHm <- function (classMat, grps = NULL, isBig = FALSE, cRow = FALSE, 
                         cCol = FALSE, fontsize_row = 4, fontsize_col = 4, main = NA, 
                         scale = FALSE, customAnnoColor = NULL, ...) 
{
   cools <- colorRampPalette(c("black", "limegreen", "yellow"))(100)
   bcol <- "white"
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


heatmapPlotly <- function(c_scoresMatrix) {
   melted_tab = data.frame("classificationScore" = NULL, 
                           "sampleName" = NULL, 
                           "tissueType" = NULL)
   
   for(sampleName in colnames(c_scoresMatrix)) {
      temp_cscore = c_scoresMatrix[, sampleName]
      
      tempTab = data.frame("classificationScore" = temp_cscore, 
                           "sampleName" = sampleName, 
                           "tissueType" = rownames(c_scoresMatrix))
      melted_tab = rbind(melted_tab, tempTab)
   }
   
   melted_tab$tissueType = factor(x = melted_tab$tissueType,levels = rev(levels(melted_tab$tissueType)))
   cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
   
   p <- ggplot(data = melted_tab, aes(sampleName, tissueType, fill = classificationScore))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = cools[1], high = cools[length(cools)], mid = cools[length(cools)/2], 
                           midpoint = 0.5, limit = c(0,1), space = "Lab", 
                           name="Classification Score") +
      xlab("Query Samples")+
      ylab("Tissue Types")+
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))
   return(ggplotly(p))    
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
   if (length(unique(melted_tab$citation)) >= 15) {
      col_num <- 5
   } else if (length(unique(melted_tab$citation)) >= 10) {
      col_num <- 3
   } else {
      col_num <- 2
   }
   
   p <- ggplot(data = melted_tab, aes(sampleName, tissueType, fill = classificationScore, text=description))+
      geom_tile(color = "white") + 
      #facet_grid(.~citation, scales="free_x") +
      facet_wrap(~citation, scales="free_x", ncol=col_num) + 
      #facet_wrap(~citation, scales="free_x") + 
      scale_fill_gradient2(low = cools[1], high = cools[length(cools)], mid = cools[length(cools)/2], 
                           midpoint = 0.5, limit = c(0,1), space = "Lab", 
                           name="Classification Score") +
      xlab("Samples")+
      ylab("Tissue Types")+
      theme_minimal()+ 
      theme(axis.text.x=element_blank()) +
      theme(axis.title.y=element_blank()) + 
      theme(panel.spacing = unit(0.5, "lines"))
   
   ggplotly(p, tooltip="all") %>% layout(margin=c(0,0,0.02,0.01))
}


heatmapRef <- function(c_scoresMatrix, sampTab) {
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
      facet_wrap(~citation, scales="free_x") + 
      scale_fill_gradient2(low = cools[1], high = cools[length(cools)], mid = cools[length(cools)/2], 
                           midpoint = 0.5, limit = c(0,1), space = "Lab", 
                           name="Classification Score") +
      xlab("Samples")+
      ylab("Tissue Types")+
      theme_minimal()+ 
      theme(axis.text.x=element_blank()) +
      theme(axis.title.y=element_blank()) + 
      theme(panel.spacing = unit(0.5, "lines"))
   return(p)
}


subsetGRNall <- function(grnAll, iGenes) {
   # Subset grnTable based on iGenes
   allTargets <- grnAll$overallGRN$grnTable$TG
   newGRNTable <- grnAll$overallGRN$grnTable[which(allTargets %in% iGenes),]
   newTFsAll <- newGRNTable$TF
   newGRNTable <- newGRNTable[which(newTFsAll %in% iGenes),]
   grnAll$overallGRN$grnTable <- newGRNTable
   
   ## Subset overallGRN graph based on iGenes
   vertex_names <- V(grnAll$overallGRN$graph)$name
   graph_iGenes <- which(vertex_names %in% iGenes)
   newGraph <- induced_subgraph(graph=grnAll$overallGRN$graph, vids=graph_iGenes, impl="copy_and_delete")
   grnAll$overallGRN$graph <- newGraph
   
   # Subset specGenes based on iGenes and tissue type
   tissueTypes <- names(grnAll$specGenes$context$general)
   newGeneral <- grnAll$specGenes$context$general
   for (tissueType in tissueTypes) {
      tissueSpecGenes <- newGeneral[[tissueType]]
      tissueSpecGenes <- tissueSpecGenes[which(names(tissueSpecGenes) %in% iGenes)]
      newGeneral[[tissueType]] <- tissueSpecGenes
   }
   grnAll$specGenes$context$general <- newGeneral
   
   # Subset ctGRN geneLists, graphLists, and tfTargets  based on iGenes and tissue type
   grnAll$ctGRNs$geneLists <- newGeneral
   
   newGraphLists <- grnAll$ctGRNs$graphLists
   for (tissueType in tissueTypes) {
      tissueGRN <- newGraphLists[[tissueType]]
      iVertices <- vertex_attr(tissueGRN, name="name")
      iVertices <- iVertices[which(iVertices %in% iGenes)]
      tissueGRN <- induced_subgraph(graph=tissueGRN, vids=iVertices, impl="copy_and_delete")
      newGraphLists[[tissueType]] <- tissueGRN
   }
   grnAll$ctGRNs$graphLists <- newGraphLists
   
   newTFTargets <- grnAll$ctGRNs$tfTargets
   for (tissueType in tissueTypes) {
      tissueTFTargets <- newTFTargets[[tissueType]]
      tissueTFTargets <- tissueTFTargets[which(names(tissueTFTargets) %in% iGenes)]
      for (TF in names(tissueTFTargets)) {
         newTargets <- tissueTFTargets[[TF]]
         newTargets <- newTargets[which(newTargets %in% iGenes)]
         tissueTFTargets[[TF]] <- newTargets
      }
      newTFTargets[[tissueType]] <- tissueTFTargets
   }
   grnAll$ctGRNs$tfTargets <- newTFTargets
   
   return(grnAll)
}

subsetTrainNormParam <- function(trainNormParam, grnAll, iGenes) {
   tissueTypes <- names(grnAll$specGenes$context$general)
   newTVals <- trainNormParam$tVals
   for (tissueType in tissueTypes) {
      newIndices <- which(names(newTVals[[tissueType]][["mean"]]) %in% iGenes)
      newTVals[[tissueType]][["mean"]] <- newTVals[[tissueType]][["mean"]][newIndices]
      newTVals[[tissueType]][["sd"]] <- newTVals[[tissueType]][["sd"]][newIndices]
   }
   trainNormParam$tVals <- newTVals
   return(trainNormParam)
}


### compute zscore
zscore<-function(x,### numeric vector
 meanVal, ### mean of distribution to compute zscore of x against 
 sdVal ### standard deviation of distribution to compute zscore of x agains
){ 
  (x-meanVal)/sdVal;
  ### zscore
}


### Compute the mean zscore of given genes in each sample
cn_zscoreVect<-function(genes, ### genes
 xvals, ### named vector
 tVals, ### tvals
 ctt ### ctt
){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
  ### zscore vector
}



#' network influence score
#'
#' Computes network influence score (NIS).
#' @param expDat ranked query expression matrix, output of running logRank()
#' @param stQuery query sample metadata table
#' @param iGenes character vector of gene names found in both training data and query data 
#' @param grnAll GRN object, output of running ccn_makeGRN()
#' @param trainNormParam normalization parameters, output of running ccn_trainNorm()
#' @param subnet name of cell/tissue type subnetwork to evaluate
#' @param ctt string indicating the cell/tissue type to compare against
#' @param colname_sid colname of stQuery with unique sample labels
#' @param relaWeight whether to weight by overall expression such that TFs with higher expression in ctt are more important (1=do the weighting)
#'
#' @return numeric matrix where rows are TFs in CT GRN and columns are query samples. Negative score indicates TF should be upregulated to acquire a ctt fate. Positive score indicates TF should be downregulated to acquire a ctt fate.
pacnet_nis <- function(expDat, stQuery, iGenes, grnAll, trainNormParam, subnet, ctt, colname_sid="sra_id", relaWeight=0) {
  tfTargList<-grnAll$ctGRNs$tfTargets
  nTargets<-vector();
  targetScore<-vector();
  tfScore<-vector();
  totalScore<-vector();
  tfWeights<-vector();
  
  tfs<-names(tfTargList[[subnet]]);
  netGenes<-names(grnAll$ctGRNs$geneLists[[subnet]])
  netGenes<-intersect(netGenes, iGenes)
  
  sids<-as.vector(stQuery[,colname_sid])
  
  ans<-matrix(0, nrow=length(tfs), ncol=nrow(stQuery))
  rownames(ans)<-tfs;
  colnames(ans)<-sids;
  
  tVals <- trainNormParam[["tVals"]]
  
  # compute a matrix of zscores.
  zzzMat<-matrix(0, nrow=length(netGenes), ncol=nrow(stQuery));
  
  for(i in seq(length(sids))){
    sid<-sids[i];
    print(paste0("Computing zscores for ", sid))
    xvals<-as.vector(expDat[netGenes, sid])
    names(xvals)<-netGenes
    zzzMat[,i]<-cn_zscoreVect(netGenes, xvals, tVals, ctt)
  }
  
  rownames(zzzMat)<-netGenes
  colnames(zzzMat)<-stQuery[,colname_sid]
  
  for(sid in sids) {
    print(paste0("TF scoring for ", sid));
    xvals<-as.vector(expDat[,sid]);
    names(xvals)<-rownames(expDat);
    
    # assign weights
    
    ### # meanVect<-unlist(tVals[[ctt]][['mean']][netGenes]);
    meanVect<-unlist(tVals[[subnet]][['mean']][netGenes]);
    meanVect <- meanVect / min(meanVect) # Added to make function work with rank-based tVals
    weights<-(2**meanVect)/sum(2**meanVect);
    
    for(i in seq(length(tfs))){
      tf<-tfs[i];
      # zscore of TF relative to target C/T
      ##      tfScore[i]<-zscore(xvals[tf], tVals[[ctt]][['mean']][[tf]], tVals[[ctt]][['sd']][[tf]]);
      
      tfScore[i]<-zzzMat[tf,sid];
      
      targs<-tfTargList[[subnet]][[tf]];
      targs<-intersect(targs, iGenes);
      
      # Zscores of TF targets, relative to C/T
      ##      tmp<-cn_zscoreVect(targs, xvals, tVals, ctt );
      tmp<-zzzMat[targs,sid];
      targetScore[i]<-sum(tmp*weights[targs]);
      
      ## new one:
      totalScore[i]<-targetScore[i] + (length(targs)*tfScore[i]*weights[tf]);
      
      if(relaWeight!=1){ # don't weight by expression
        meanW<-mean(weights)
        totalScore[i]<- sum(tmp)*meanW + (length(targs)*tfScore[i])*meanW
      }
      nTargets[i]<-length(targs) ;
      tfWeights[i]<-weights[tf];
    }
    
    xxx<-data.frame(tf=tfs, tfScore=tfScore, targetScore=targetScore, nTargets=nTargets,tfWeight=tfWeights, totalScore=totalScore);
    xxx<-xxx[order(xxx$totalScore),]; # puts the worst ones at top when plotting
    xxx$tf<-factor(xxx$tf, as.vector(unique(xxx$tf)));
    ans[as.vector(xxx$tf),sid]<-as.vector(xxx$totalScore);
  }
  
  return(ans) # returns network influence score.
}

