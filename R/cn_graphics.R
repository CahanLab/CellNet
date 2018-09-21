# CellNet
# (C) Patrick Cahan 2012-2016


#' make a rainbow colored dot plot
#'
#' make a rainbow colored dot plot
#' @param expDat expression data matrix
#' @param stAll sample table
#' @param gene gene name
#' @param dLevel column name to group samples by
#'
#' @return ggplot object
#'
#' @examples 
#' mp_rainbowPlot(expDat,sampTab, "Actb", "description1")
#'
#' @export
mp_rainbowPlot<-function### make a rainbow colored dot plot
(expDat,
 stAll,
 gene,
 dLevel="description1"
){
  xDat<-cbind(stAll, gene=expDat[gene,]);
  xi<-which(colnames(xDat)==dLevel);
  colnames(xDat)[xi]<-'type';
  xplotb<-ggplot(xDat, aes(x=type, y=gene, color=type)) + 
    geom_point(position='jitter', shape=19,alpha=5/8, size=.7,show.legend=F) + 
    theme_bw() + coord_flip() + ylab(gene) + xlab("");
  xplotb;
}

#' boxplot of network influence scores
#'
#' NIS<0 means that the expression of the TF is lower in query samples as compared to target cell type, and magnitude determined by number and dysregualtion of target genes
#' @param tfScores result of running cn_nis_all
#' @param targetCT what is the target cell type? 
#' @param sampTab sample table to select subset to plot
#' @param group which subset to plot
#' @param dLevel level to group on
#' @param limitTo number of TFs to show
#'
#' @return ggplot object of box and whiskers
#'
#' @examples
#' plot_nis(tfScores, "hspc", sampTab, "day21")
#'
#' @export
plot_nis<-function#### boxplot of NIS scores, requires plyr, tidyr
(tfScores,
 targetCT,
 sampTab,
 group,
 limitTo=20,
 dLevel="description1"
){

  stTmp<-sampTab[sampTab[,dLevel]==group,]
  tfScores<-tfScores[[targetCT]][,rownames(stTmp)];

  xx<-as.data.frame(t(tfScores))
  cnames<-colnames(xx)
  newX<-gather_(xx, "gene", "expression", cnames)
  newx2<-transform(newX, gene=reorder(gene, -expression))
  newx3<-ddply(newx2, "gene", transform, medVal=median(expression, na.rm=TRUE))
  
  nTFs<-length(unique(newx3$gene))
  if(limitTo==0 | limitTo>nTFs){
    limitTo<-nTFs
  }

  genes <- unique(as.vector(newx3[order(abs(newx3$medVal), decreasing=TRUE),]$gene))[1:limitTo]
  tmpAns <- newx3[which(newx3$gene==genes[1]),]
  for(gene in genes[2:length(genes)]){
    tmpAns<-rbind(tmpAns,newx3[which(newx3$gene==gene),])
  }

  ggplot(tmpAns, aes(x=gene, y=expression)) + 
    geom_boxplot(aes(fill=medVal)) + coord_flip() + theme_bw() + 
    scale_fill_gradient2(low='purple', mid='white', high='orange') + 
    ylab("Network influence score") + xlab("Transcriptional regulator") + theme(legend.position="none", axis.text=element_text(size=8))
}


#' Plot GRN status
#'
#' wrapper to barplot secific GRN
#' @param cnObj result of analyzing query data with cn_apply
#' @param cnProc result of creating a cn_make_processor
#' @param snName subnet name of which to plot establishment or status level 
#' @param ctrlSamps names of samples in training data
#' @param bOrder order of bars
#' @param dLevel which stquery level to plot on
#' @param sidCol sample id colname
#'
#' @return
#'
#' @examples
#' cn_barplot_grnSing(cnRes, cnProc, "hspc", c("esc", "hspc"), bOrder=c("esc_train", "day0", "day5", "day10", "day20", "hspcxs_train"))
#'
#' @export
cn_barplot_grnSing<-function ### wrapper to barplot secific GRN
(cnObj,
 cnProc,
 snName,
 ctrlSamps,
 bOrder,
 dlevel='dLevelQuery',
 sidCol="sample_id"
){
  
  if(dlevel=='dLevelQuery'){
    dlevel<-cnObj[['dLevelQuery']]
  }

  qScores<-cnObj[['normScoresQuery']];
  ctrlScores<-cnProc[['trainingScores']]
  
  ###  .cn_barplot_grnSing(qScores, ctrlScores, cnObj[['stQuery']], cnObj[['dLevelQuery']], snName, ctrlSamps, bOrder);
  cn_barplot_grnSing_base(qScores, ctrlScores, cnObj[['stQuery']], dlevel, snName, ctrlSamps, bOrder, sidCol=sidCol);
}

#' barplot this specific GRN
#'
#' barplot this specific GRN
#' @param qScores queryScores
#' @param control scores
#' @param stQuery sample table
#' @param dLevelQ dLevel of query samples
#' @param snName name of subnet to plot establishment level 
#' @param ctrSamps names of samples in training data
#' @param bOrder order of bars
#' @param sidCol sample id colname
#'
#' @return gpplot barplot
cn_barplot_grnSing_base<-function### 
(qScores, 
 ctrlScores,
 stQuery,
 dLevelQ,
 snName,
 ctrSamps,
 bOrder,
 sidCol="sample_id"
 ){
    
  # convert into a data.frame
  aa<-cn_extract_SN_DF(qScores, stQuery,dLevelQ, rnames=snName, sidCol=sidCol);
  aa3<-cn_reduceMatLarge(aa, "score", "description", "subNet");
  aa3<-cbind(aa3, src=rep('query', nrow(aa3)));
  tmpAns<-data.frame();
  for(ctrSamp in ctrSamps){
    xxx<-ctrlScores[ctrlScores$grp_name==ctrSamp & ctrlScores$subNet==snName,];
    xxx$grp_name<-paste(xxx$grp_name, "_train", sep='');
    tmpAns<-rbind(tmpAns, xxx);
  }
  tmpAns<-cbind(tmpAns, src=rep("train", nrow(tmpAns)));
  
  if(is.null(bOrder)){
    bOrder<-c(as.vector(tmpAns$grp_name), as.vector(aa3$grp_name));
    ##bOrder<-aa3$grp_name[order(aa3$mean, decreasing=TRUE)];    
  }
  
  aa3<-rbind(aa3, tmpAns);
  aa3$grp_name<-factor(aa3$grp_name, bOrder);
  
  ##
  # convert is.na(stdev) -> 0
  xi<-which(is.na(aa3$stdev));
  if(length(xi)>0){
    aa3[xi,'stdev']<-0;
  }
  # # 
  ans<-  ggplot(na.omit(aa3), aes(x=grp_name, y=mean, fill=src)) +
    geom_bar(width=.75,position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev),width=.2,position=position_dodge()) +
    scale_fill_brewer(palette = "Paired")  +
    theme_bw() +
    theme(text = element_text(size=8), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
    ggtitle(snName) + theme(axis.title.x = element_blank())+ ylab("GRN status")
  ans;
}


#' heatmap of the classification result
#'
#' Heatmap of the classification result.
#' @param cnRes returned from cn_sapply
#' @param isBig is this a big heatmap? TRUE or FALSE
#'
#' @return nothing
#'
#' @examples
#' cn_HmClass(cnRes, isBig=TRUE)
#'
#' @export
cn_HmClass<-function
(cnRes, 
 isBig=FALSE
){
 
  classMat<-cnRes$classRes;
  cools<-colorRampPalette(c("black", "limegreen", "yellow"))( 100 )
  bcol<-'white';
  if(isBig){
    bcol<-NA;
  }
  pheatmap(classMat,
    col=cools,
    breaks=seq(from=0, to=1, length.out=100),
    border_color=bcol,
    cluster_rows = FALSE,
    cluster_cols = FALSE)
  # classification heatmap
}

#' heatmap of gene expression
#'
#' heatmap of gene expression
#' @param expDat expression matrix
#' @param sampTab sample table
#' @param isBig boolean
#' @param dist distance metricx used, determines also whether to row-scale or not
#'
#' @return nothing
#'
#' @examples
#' selectedGenes<-c("Actb", "Gapdh", "Hlf", "Sox2")
#' cn_HmVars(expDat[selectedGenes,], sampTab)
#'
#' @export
cn_HmVars<-function
(expDat, 
 sampTab,
 isBig=FALSE,
 dist='correlation'
){

 if(dist=='correlation'){
    scale<-'row'
  }
  else{
    scale<-'none'
  }

  ##cools<-colorRampPalette(c("blue", "white", "red"))( 100 )
  bcol<-'white';
  if(isBig){
    bcol<-NA;
  }
  pheatmap(expDat,
   # annotation_col=sampTab,
    annotation_legend = FALSE,
    show_colnames = F,
   ### col=cools,
    border_color=bcol,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    dist=dist,
    scale=scale)
  # 
}


#' heatmap of gene expression for disease modeling studies
#' Added 6-4-18
#' @param cnResQuery cellNet query results object
#' @param isBig boolean
#'
#' @return nothing
#'
#' @examples
#' cnRes_iPSC <- utils_loadObject("cn_iPSC_pooled_results_nat_prot_Jul_14_2017.R")
#' col_annotations_iPSC <- data.frame(Disease_Status = factor(sampTab_iPSC$disease_status), Study = factor(sampTab_iPSC$study_id))
#' set.seed(0) # reseed random number generator, makes column coloring reproducible
#' cn_HmClass_disease(cnRes_iPSC, col_annotations_iPSC)
#'
#' @export
cn_HmClass_disease<-function (cnResQuery, col_annotations, isBig = TRUE) 
{
   library(pheatmap)
   library(RColorBrewer)
   classMat <- cnResQuery$classRes
   newSampTab <- cnResQuery$stQuery
   colnames(classMat)<-make.unique(colnames(classMat))
   rownames(col_annotations)<-colnames(classMat)
   
   study_names<-unique(newSampTab$study_id)
   breaks<-match(study_names, newSampTab$study_id)-1
   study_colors<-sample(rainbow(length(study_names)))
   #study_colors<-sample(brewer.pal(length(study_names), "Set3"))
   names(study_colors)<-study_names
   ann_colors <- list(
      Disease_Status = c(control="#a6cee3",disease="#e78ac3",treated="#d699ff"), 
      Study = study_colors)
   hm_colors <- colorRampPalette(c("black", "limegreen", "yellow"))(100)
   if (isBig) {
      bcol <- NA
   }
   bcol <- "white"
   pheatmap(classMat,
            annotation_col=col_annotations,
            gaps_col = breaks,
            col = hm_colors, 
            annotation_colors=ann_colors,
            breaks = seq(from = 0, to = 1, length.out = 100), 
            border_color = bcol, 
            cluster_rows = FALSE, 
            cluster_cols = FALSE,
            show_colnames = FALSE)
}

#' Output ordered classifier performance heatmap by description1 label of validation samples
#' Added 6-4-18
#' @param cnProc cellNet processor object from training
#' @param classifierPerformance cellNet classifier performance object from training
#'
#' @return nothing
#'
#' @examples
#' ordered_class_perf_hm(cnProcFile="cnProc_HS_RS_Jun_20_2017.rda", classPerfFile="classifierPerformance_June-20-2017.rda")
#'
#' @export
ordered_class_perf_hm<-function(cnProc, classifierPerformance) {
   
   # Order stTrain by description1, then pull out the ordered SRRs
   stTrain<-cnProc$stTrain
   stTrain_desc1_ordered<-stTrain[order(stTrain$description1),]
   stTrain_srrs<-stTrain_desc1_ordered$sra_id
   
   # Pull out the validation SRRs
   validation_srrs<-colnames(classifierPerformance$classRes)
   
   # Get the ordered indices of the validation SRRs using the ordered SRRs from stTrain
   ordered_srr_indices<-match(stTrain_srrs, validation_srrs)
   ordered_srr_indices<-ordered_srr_indices[!is.na(ordered_srr_indices)]
   
   # Make a new classifier Performance object with an ordered classRes
   newClassPerf<-classifierPerformance
   newClassPerf$classRes<-classifierPerformance$classRes[,ordered_srr_indices]
   
   # Make heatmap
   cn_HmClass(newClassPerf, isBig = TRUE)
   
}



if(FALSE){
cn_barplot_exp<-function
### barplot gene expression
(cnObj,
 ### result of analyzing query data with CN
 cnProc,
 ### result of creating a cellnet processor
 gName,
 ### name of gene to plot
 ctrSamps,
 ### names of samples in training data. e.g. c("hspc", "liver")
 bOrder=NULL,  
 ### order of bars
 logRatio=NULL
 ### if !null, expression expr compared to mean value of logRatio values 
 ){
  
  ddl<-cnProc[['dLevelTrain']];
  qScores<-cnObj[['expQuery']];
  stTrain<-cnProc[['stTrain']];
  
  x<-data.frame();
  for(ct in ctrSamps){
      x<-rbind(x, stTrain[which(stTrain[,ddl]==ct),]);  
  }

  expTrain<-cnProc[['expTrain']];
  expTrain<-expTrain[,rownames(x)];
  stTrain<-x;
  stTrain<-stTrain[,c("sample_id", "sample_name", ddl)];
  colnames(stTrain)[3]<-"grp_name";
  stTrain[,3]<-paste(stTrain[,"grp_name"], "_train", sep='');
  expX<-cbind(stTrain, expTrain[gName,]);
  colnames(expX)[4]<-gName;
  
  ddq<-cnObj[['dLevelQuery']];
  stQuery<-stQuery[,c("sample_id", "sample_name", ddq)];
  nQs<-length(unique(stQuery[,ddq]));
  
  colnames(stQuery)[3]<-"grp_name";
  expQ<-cbind(stQuery, expQuery[gName,]);
  colnames(expQ)[4]<-gName
  expNew<-rbind(expX, expQ);
  cat("2\n")
  expQ<-utils_reduceMat(expNew, gName, "grp_name");
  newExp<-cbind(expQ, src=c( rep("train", length(ctrSamps)), rep("query",nQs)));
  
  
  if(is.null(bOrder)){
      bOrder<-as.vector(newExp$grp_name);
  }
  newExp$grp_name<-factor(newExp$grp_name, bOrder);

  

  #### 
  
  ##
  # convert is.na(stdev) -> 0
  xi<-which(is.na(newExp$stdev));
  if(length(xi)>0){
    newExp[xi,'stdev']<-0;
  }

  ans<-  ggplot(na.omit(newExp), aes(x=grp_name, y=mean, fill=src)) +
    geom_bar(width=.75,position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev),width=.2,position=position_dodge()) +
    scale_fill_brewer(palette = "Paired")  +
    theme_bw() +
    theme(text = element_text(size=8), axis.text.x = element_text(angle=90, vjust=1)) +
    ggtitle(gName) + theme(axis.title.x = element_blank())+ ylab("Expression")
  ans;
  newExp;
  ### single gene barplot
}
}

#' Plot results of cn_classAssess
#'
#' Plot one precision recall curve per CT
#' @param assessed result of runnung cn_classAssess
#'
#' @return ggplot pbject
#'
#' @examples
#' testAssTues<-cn_splitMakeAssess(stTrain, expTrain, ctGRNs, prop=.5)
#' plot_class_PRs(testAssTues$ROCs)
#'
#' @export
plot_class_PRs<-function
(assessed
  ){
  ctts<-names(assessed);
  df<-data.frame();
  for(ctt in ctts){
    tmp<-assessed[[ctt]];
    tmp<-cbind(tmp, ctype=ctt);
    df<-rbind(df, tmp);
  }

  prsAll<-transform(df, TP = as.numeric(as.character(TP)), 
    TN = as.numeric(as.character(TN)), 
    FN = as.numeric(as.character(FN)), 
    FP = as.numeric(as.character(FP)));

    precfunc<-function(df){
      ans<-vector();
      for(i in 1:nrow(df)){
        ans<-append(ans, df[i,"TP"]/(df[i,"TP"]+df[i,"FP"]));
      }
      ans;
    }

    sensfunc<-function(df){
      ans<-vector();
      for(i in 1:nrow(df)){
        ans<-append(ans, df[i,"TP"]/(df[i,"TP"]+df[i,"FN"]));
      }
      ans;
    }

  precs<-precfunc(prsAll)
  sens<-sensfunc(prsAll)
  prsAll2<-cbind(prsAll, data.frame(recall=sens, precision=precs));

  ggplot(data=prsAll2, aes(x=as.numeric(as.vector(recall)), y=as.numeric(as.vector(precision)))) + geom_point(size = .5, alpha=.5) +  geom_path(size=.5, alpha=.75) +
  theme_bw() + xlab("Recall") + ylab("Precision") + facet_wrap( ~ ctype, ncol=4) +
  theme(axis.text = element_text(size=5)) + ggtitle("Classification performance")
}


#' Order classifier performance heatmap by description1 label of validation samples
#' 
#' @param cnProc CellNet Processor object from training
#' @param classifierPerformance Classifier performance object generated with cn_splitMakeAssess
#' 
#' @export
cn_class_perf_hm_ordered<-function(cnProc, classifierPerformance) {
   
   # Order stTrain by description1, then pull out the ordered SRRs
   stTrain<-cnProc$stTrain
   stTrain_desc1_ordered<-stTrain[order(stTrain$description1),]
   stTrain_srrs<-stTrain_desc1_ordered$sra_id
   
   # Pull out the validation SRRs
   validation_srrs<-colnames(classifierPerformance$classRes)
   
   # Get the ordered indices of the validation SRRs using the ordered SRRs from stTrain
   ordered_srr_indices<-match(stTrain_srrs, validation_srrs)
   ordered_srr_indices<-ordered_srr_indices[!is.na(ordered_srr_indices)]
   
   # Make a new classifier Performance object with an ordered classRes
   newClassPerf<-classifierPerformance
   newClassPerf$classRes<-classifierPerformance$classRes[,ordered_srr_indices]
   
   # Make heatmap
   cn_HmClass(newClassPerf, isBig = TRUE)
   
}



########### Plot pdfs
#' Added 6-4-18


#' Create pdf of classification heatmap, one column for every sample.
#' 
#' @param cnResQuery CellNet query results object
#' @param study_name desired study name
#' 
#' @example plot_classification_hm(cnResQuery, "Study_1")
#' 
#' @export
pdf_classification_hm <- function
(cnResQuery, 
 study_name
){
   mydate<-utils_myDate()
   newSampTab <- cnResQuery$stQuery 
   height = 8
   width = 4 + length(newSampTab$sample_name) / 3
   
   fname<-paste("hm_classification_query_", study_name, "_rna-seq_", mydate,".pdf", sep="")
   pdf(file=fname, width=width, height=height, onefile=FALSE)
   cn_HmClass(cnResQuery)
   dev.off()
}

#' Create pdf of GRN status of starting and target cell types, grouping samples by dlevel
#' 
#' @param cnResQuery CellNet query results object
#' @param cnProc CellNet Processor object from training
#' @param study_name desired study name
#' @param cell_type_1 starting cell/tissue type
#' @param cell_type_2 ending or target cell/tissue type
#' @param dlevel column of sample table used for sample grouping
#' 
#' @example plot_grn_status(cnResQuery, cnProc, "esc", "heart", "Study_1", dlevel="description1")
#' 
#' @import ggplot2
#' 
#' @export
pdf_grn_status <- function
(cnResQuery, 
 cnProc, 
 cell_type_1, 
 cell_type_2, 
 study_name,
 dlevel="description2"
){
   mydate<-utils_myDate()
   newSampTab <- cnResQuery$stQuery  
   height = 5
   width = 4 + length(unique(newSampTab$description2)) / 2
   
   print("Plotting GRN Status...")
   bOrder<-c(paste(cell_type_1, "_train", sep=""), 
             unique(as.vector(newSampTab$description2)), 
             paste(cell_type_2, "_train", sep=""))
   cn_barplot_grnSing(cnResQuery, cnProc, cell_type_1, c(cell_type_1, cell_type_2), 
                      bOrder, sidCol="sra_id", dlevel=dlevel)
   fname <- paste("grnStatus_", study_name, "_", cell_type_1, "_", mydate, ".pdf", sep="");
   ggplot2::ggsave(fname, width=width, height=height)
   
   
   bOrder<-c(paste(cell_type_1, "_train", sep=""), 
             unique(as.vector(newSampTab$description2)), 
             paste(cell_type_2, "_train", sep=""))
   cn_barplot_grnSing(cnResQuery, cnProc, cell_type_2, c(cell_type_1, cell_type_2), 
                      bOrder, sidCol="sra_id", dlevel=dlevel);
   fname <- paste("grnStatus_", study_name, "_", cell_type_2, "_", mydate, ".pdf", sep="");
   ggplot2::ggsave(fname, width=width, height=height)
   
}

#' Create multiplot - adapted from Bioconductor package "scater"
#' 
#' @import grid
#' 
#' @example multiplot(plotlist = grn_plot_list, layout=plot_layout)
#' 
multiplot <- function (plotlist = NULL, layout = NULL, cols = 1) 
{
   library(grid)
   plots <- plotlist
   numPlots = length(plots)
   if (is.null(layout)) {
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), 
                       ncol = cols, nrow = ceiling(numPlots/cols))
   }
   if (numPlots == 1) {
      print(plots[[1]])
   }
   else {
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), 
                                                 ncol(layout))))
      for (i in 1:numPlots) {
         matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
         print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, 
                                         layout.pos.col = matchidx$col))
      }
   }
}

#' Create pdf of GRN status plots, one plot with all dlevel groups for every C/T type
#' 
#' @param cnResQuery CellNet query results object
#' @param cnProc CellNet Processor object from training
#' @param study_name desired study name
#' @param dlevel column of sample table used for sample grouping
#' 
#' @example plot_grn_status_by_CT(cnResQuery, cnProc, "Study_1", dlevel="description1")
#' 
#' 
#' @export
pdf_grn_status_by_CT <- function
(cnResQuery, 
 cnProc, 
 study_name,
 dlevel="description2"
){
   mydate<-utils_myDate()
   newSampTab<-cnResQuery$stQuery
   tissue_types <- rownames(summary(cnProc$ctGRNs$ctGRNs$geneLists))
   tissue_types<-tissue_types[order(tissue_types)]
   
   
   grn_plot_list<-vector("list", length(tissue_types))
   i=1
   for (tissue in tissue_types) {
      cat(paste("Plotting GRN status for ", tissue, "...\n", sep=""))
      bOrder<-c(paste(tissue, "_train", sep=""), unique(as.vector(newSampTab[[dlevel]])))
      grn_plot<-cn_barplot_grnSing(cnResQuery, cnProc, tissue, tissue, bOrder, sidCol="sra_id", dlevel=dlevel)
      grn_plot_list[[i]]<-grn_plot
      i=i+1
   }
   grn_plot_list<-as.list(grn_plot_list)
   
   fname <- paste("grn_status_by_CT_", study_name, "_", mydate, ".pdf", sep="");
   n_plot_rows<-ceiling(length(tissue_types) / 4)
   height = 4 * n_plot_rows
   width = 16
   
   cat("Finishing plot...\n")
   pdf(fname, height=height, width=width)
   num_empty<-n_plot_rows * 4 - length(tissue_types)
   plot_layout<-matrix(c(1:length(tissue_types), rep(0, (num_empty))), ncol=4, byrow=TRUE)
   multiplot(plotlist = grn_plot_list, layout=plot_layout)
   dev.off()
   cat("Done!\n")
}

#' Create pdf of GRN status plots, one plot with all C/T types for every dlevel group
#' 
#' @param cnResQuery CellNet query results object
#' @param cnProc CellNet Processor object from training
#' @param study_name desired study name
#' @param target_cell_type ending or target cell/tissue type
#' @param dlevel column of sample table used for sample grouping
#' 
#' @example plot_grn_status_by_dlevel(cnResQuery, cnProc, "Study_1", "heart", dlevel="description1")
#' 
#' @import tibble
#' @import grid
#' @import gridExtra
#' @import reshape2
#' @import ggplot2
#' 
#' @export
pdf_grn_status_by_dlevel <- function
(cnResQuery, 
 cnProc,
 study_name, 
 target_cell_type,
 dlevel="description2"
){
   mydate<-utils_myDate()
   newSampTab<-cnResQuery$stQuery
   tissue_types <- rownames(summary(cnProc$ctGRNs$ctGRNs$geneLists))
   tissue_types<-tissue_types[order(tissue_types)]
   
   qScores <- as.data.frame(cnResQuery$normScoresQuery)
   
   dlevel_names<-unique(newSampTab[,dlevel])
   descrip_grn_plot_list<-vector("list", length(dlevel_names))
   i = 1
   for (descrip in dlevel_names) {
      cat(paste("Plotting GRN status for ", descrip, " ...\n", sep=""))
      descrip_indices<-which(newSampTab[[dlevel]] == descrip)
      descrip_grn_scores<-get_grn_scores(newSampTab, descrip_indices, qScores)
      plot_df<-data.frame(matrix(ncol=2, nrow=length(tissue_types)))
      for (j in 1:nrow(plot_df)){
         plot_df[j,1]<-mean(as.numeric(descrip_grn_scores[j,]))
         plot_df[j,2]<-sd(as.numeric(descrip_grn_scores[j,]))
      }
      rownames(plot_df)<-rownames(descrip_grn_scores)
      
      target_train_scores<-cnProc$trainingScores[cnProc$trainingScores$grp_name == target_cell_type,]
      plot_df<-cbind(plot_df, target_train_scores$mean, target_train_scores$stdev)
      
      colnames(plot_df)<-c(paste(descrip, "mean"), paste(descrip, "st dev"), 
                           paste(target_cell_type, "train mean"), paste(target_cell_type, "train stdev"))
      means_df<-plot_df[,c(1,3)]
      stdevs_df<-plot_df[,c(2,4)]
      
      stdevs_df<-tibble::rownames_to_column(stdevs_df, var="tissue_type")
      melted_stdevs_df<-melt(stdevs_df)
      
      means_df<-tibble::rownames_to_column(means_df, var="tissue_type")
      melted_means_df<-melt(means_df)
      
      y_maxes<-melted_means_df$value + melted_stdevs_df$value
      y_mins<-melted_means_df$value - melted_stdevs_df$value
      
      descrip_plot<-ggplot(melted_means_df, aes(x=reorder(tissue_type, -value), y=value, fill=variable)) + 
         geom_bar(stat="identity", position = "dodge") + 
         geom_errorbar(mapping=aes(ymax = y_maxes, ymin = y_mins), position = "dodge", width = 0.9) +
         ylab("GRN Status Score") + scale_fill_manual(values = c("#a6cee3", "#2952a3")) +
         xlab("Tissue Type GRN Status") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
         theme(legend.position = "bottom") + ggtitle(descrip)
      descrip_grn_plot_list[[i]]<-ggplotGrob(descrip_plot)
      i = i + 1
   }
   
   fname <- paste("grn_status_by_", dlevel, "_", target_cell_type,"_", study_name, "_", mydate, ".pdf", sep="");
   n_plot_rows<-ceiling(length(dlevel_names) / 4)
   height = 6 * n_plot_rows
   width = 24
   
   cat("Finishing plot...\n")
   pdf(fname, height=height, width=width)
   grid.arrange(grobs=descrip_grn_plot_list, ncol=4)
   dev.off()
   cat("Done!\n")
}

#' Internal function, subset GRN scores for given sample indices
#' 
#' @param sampTab sample metadata table
#' @param sample_indices integer index vector of which samples to subset
#' @param qScores matrix of study GRN scores
#' 
get_grn_scores<-function
(sampTab, 
 sample_indices, 
 qScores
){
   samptab_subset<-sampTab[sample_indices,]
   subset_srrs<-samptab_subset$sra_id
   subset_grn_scores<-qScores[subset_srrs]
   #Check if you got the right columns
   #subset_check<-filter(newSampTab, sra_id %in% colnames(subset_grn_scores))
   return(subset_grn_scores)
}

#' Create pdf of Network Influence Scores (NISs) for given C/T type, grouping samples by dlevel
#' @param cnResQuery CellNet query results object
#' @param cnProc CellNet Processor object from training
#' @param study_name desired study name
#' @param target_cell_type ending or target cell/tissue type
#' @param dlevel column of sample table used for sample grouping
#' 
#' @export
pdf_nis_plots <- function
(cnResQuery, 
 cnProc, 
 target_cell_type, 
 study_name,
 dlevel="description2"
){
   mydate<-utils_myDate()
   newSampTab<-cnResQuery$stQuery
   rownames(newSampTab)<-as.vector(newSampTab$sra_id)
   tfScores<-cn_nis_all(cnResQuery, cnProc, target_cell_type);
   
   fname <- paste("cn_tfScores_", target_cell_type, "_", study_name, "_", mydate, ".R", sep="")
   save(tfScores, file=fname)
   
   groups <- unique(newSampTab$description2)
   fname<-paste("NIS_", study_name, "_", target_cell_type, "_", mydate, ".pdf", sep="")
   pdf(fname)
   for (group in groups) {
      print(paste("Plotting NIS for ", group, sep=""))
      group_plot <- plot_nis(tfScores, target_cell_type, newSampTab, group, dLevel=dlevel, limitTo=0) + ggtitle(group)
      print(group_plot)
   }
   dev.off()
}


