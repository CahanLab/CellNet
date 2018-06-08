newHm <- function
### heatmap of the classification result
(cnRes, 
 ### cellnet result
 isBig=FALSE
 ### is this a big heatmap
){
  classMat<-cnRes$classRes
  ## SORT BY TISSUE TYPE
  classMat = classMat[,order(colnames(classMat))]
  
  ## MAP COLOR TO TISSUE TYPE
  columns = unique(colnames(classMat))
  counts = list()
  for(i in 1:length(columns)) {
    counts[[i]] = length(which(colnames(classMat) == columns[[i]]))
  }
  annotation_col = data.frame(TissueType = factor(rep(columns, counts)))
  ann_colors = list(c(columns = rainbow(length(columns))))
  ## 
  colnames(classMat)  =  rownames(annotation_col)
  
  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
  bcol<-'white';
  if(isBig){
    bcol<-NA;
  }
  
  # GENERATE HEATMAP
  pheatmap(classMat,
           col=cools,
           border_color=bcol,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_colnames = FALSE,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           )
  ##library(grid)
  ##grid.ls(grid.force())
  ##grid.gedit("col_annotation", gp = gpar(col="grey70")) 
 
}


