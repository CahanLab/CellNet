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



if(FALSE){
cn_stdout<-function
### save pdfs of standard output
(cnObj,
 ### cnRes object
 cnProc,
 ### CellNet object used to produce cnObj
 tfScores,
 ### result of running cn_nis
 fname_prefix,
 ### what to put at the front of the pdffile name
 targetType
 ### what is the target cell type, must be one of the cell types listed in CellNet object
){
  
  ##<<note function to produce a 'standard' output consisting of:
  ##<< (1) classification heatmap
  ##<< (2) starting and target TCT GRN establishments
  ##<< (3) Aberrant TCT GRN establishments
  ##<< (4) all TCT GRN establishments
  findWidth<-function(panels){
    #  panels<-tmp[['npanels']];
    height<-ceiling(panels/4)*3;
    if(panels<4){
      width<-panels*3;
    }
    else{
      width<-12;
    }
    list(height=height, width=width)
  }
  
  stQuery<-cnObj[['stQuery']];
  cttBest<-cnProc[['grnList']];
  dLevel<-cnObj[['dLevelQuery']];
  
  # Classification heatmap
  myWidth<-nrow(stQuery)*.4;
  myHeight<-length(cttBest)*.3;
  
  xtime<-format(Sys.time(), "%Y_%b_%d_%H_%M_%S");
  fname<-paste(fname_prefix, "_", xtime, "_plots.pdf",sep='');
  tempTitle<-paste("Classification heatmap ", fname_prefix, sep='');
  pdf(fname, width=8.5, height=11);
  cn_hmClass(cnObj, main=tempTitle);
  
  # set this up so that only the training data from the 'grn' cell type is also shown,
  # along with the query samples
  grnNames<-rownames(cnObj[['normScoresQuery']]);
  for(grnName in grnNames){    
    tSamp<-paste(grnName, "_train", sep='');
    print(cn_barplot_grnSing(cnObj, cnProc, grnName, grnName, bOrder=NULL, norm=T));    
  }
  # plot transcriptional regulator scores
  # for now, this only looks at the target cell/tissue type and is a heatmap
##  tfS<-tfScores[which(tfScores$grn==targetType),];
##  tfS<-tfS[,1:(ncol(tfS)-1)];
  print(cn_plotnis(tfScores[[targetType]]));
  dev.off();
  fname;
  ### return the file name of the plot pdf file.
}
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


