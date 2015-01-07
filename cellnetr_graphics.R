
mp_rainbowPlot<-function### make a rainbow colored dot plot
(expDat,### expression data matrix
 stAll, ### sample table
 gene, ### gene name
 dLevel="description2" ### stAll column name to group samples by
){
  xDat<-cbind(stAll, gene=expDat[gene,]);
  xi<-which(colnames(xDat)==dLevel);
  colnames(xDat)[xi]<-'type';
  xplotb<-ggplot(xDat, aes(x=type, y=gene, color=type)) + 
    geom_point(position='jitter', shape=19,alpha=5/8, size=.7,show_guide=F) + 
    theme_bw() + coord_flip() + ylab(gene) + xlab("");
  xplotb;
  # rainbow plot
}

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

cn_barplot_grnSing<-function
### barplot this specific GRN
(cnObj,
 ### result of analyzing query data with CN
 cnProc,
 ### result of creating a cellnet processor
 snName,
 ### name of subnet to plot establishment level 
 ctrSamps,
 ### names of samples in training data
 bOrder,  
 ### order of bars
 norm=FALSE 
 ### normalize?
){
  
  qScores<-cnObj[['queryScores']];
  if(norm){  
    nVals<-cnProc[['normVals']];  
    qScores<-cnObj[['normScoresQuery']];
    #ctrlScores<-cnProc[['ctrlScores']];
    ctrlScores<-cnProc[['trainingScores']]
  }
  else{
    ctrlScores<-cnProc[['raw_scores']];
    xx<-cn_extract_SN_DF(ctrlScores, cnProc[['stTrain']], cnProc[['dLevelTrain']],rnames=snName);
    xx3<-cn_reduceMatLarge(xx, "score", "description", "subNet");
    ctrlScores<-xx3;
  }
  
  # convert into a data.frame
  aa<-cn_extract_SN_DF(qScores, cnObj[['stQuery']], cnObj[['dLevelQuery']], rnames=snName);
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
  
  #### 
  
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
    theme(text = element_text(size=8), axis.text.x = element_text(angle=90, vjust=1)) +
    ggtitle(snName) + theme(axis.title.x = element_blank())+ ylab("GRN status")
  ans;
  ### single GRN barplot
}


cn_hmClass<-function
### plot classification heatmap
(cn, 
 ### CellNet Object
 main=NULL, 
 ### title
 dLevel=NULL,
 ### dLevel default=NULL
 isBig=FALSE
 ### is big? =FALSE
 ){
  
  classRes<-cn$classRes;
  sampTab<-cn$stQuery;
  
  if(!is.null(dLevel)){
    colnames(classRes)<-as.vector(sampTab[,dLevel]);
  }
  .cn_HmClass(classRes,main, isBig=isBig);  
  ### classification heatmap
}

.cn_HmClass<-function
### heatmap of the classification result
(classMat, 
 ### aMat columns are samples, rows are reference cell types
 main=NULL, 
 ### title
 sampTab=NULL,  
 ### sample table
 dLevel="description1", 
 ### level at which to extract column labels
 scale='none',
 ### scale?
 clusterR=FALSE,
 ### cluster the rows?
 clusterC=FALSE,
 ### cluster the columns?
 margin=c(6,6),
 ### margin
 cexR=.75,
 ### cexR
 cexC=.5,
 ### cexC
 aDist=dist,
 ### distance metric
 isBig=FALSE
 ### is Big? false
){
  
  if(is.null(sampTab)){
    labCol=colnames(classMat);
  }
  else{
    labCol<-as.vector(sampTab[,dLevel]);
  }
  if(isBig){
    heatmap.2(classMat,
              col=colorpanel(99, "black","limegreen","yellow"),
              breaks=seq(from=0, to=1,length.out=100 ),
              #col=colorpanel(100, "black","yellow"),
              scale=scale,
              trace='none',
              key=T,
              Rowv=clusterR,
              Colv=clusterC,
              labCol=labCol,
              density.info='none',
              cexRow=cexR,
              cexCol=cexC,
              margin=margin,
              dist=aDist,
              main=main);  
  }
  else{
    heatmap.2(classMat,
              col=colorpanel(99, "black","limegreen","yellow"),
              breaks=seq(from=0, to=1,length.out=100 ),
              #col=colorpanel(100, "black","yellow"),
              scale=scale,
              trace='none',
              key=T,
              Rowv=clusterR,
              Colv=clusterC,
              labCol=labCol,
              density.info='none',
              cexRow=cexR,
              cexCol=cexC,
              margin=margin,
              colsep=seq(ncol(classMat)),
              rowsep=seq(nrow(classMat)),
              sepcol='white',
              sepwidth=c(0.001,0.001),
              dist=aDist,
              main=main);  
  }
  # classification heatmap
}


# generic red/blue heatmap
mp_hmVars<-function# basic heatmap
(expDat,
 ### numeric matrix
 genes,
 ### rownames to include in the heatmap
 main='',
 ### optional title for teh top of the HM
 clusterR=T,
 ### whether to cluster the Rows
 clusterC=F,
 ### whether to cluster the columns
 scale='row',
 ### normalize the data: 'row', 'column', or 'none'
 big=FALSE,
 ### if not big then add cell separators
 dist=utils_myDist,
 ###
 margin=c(12,6),
 ###
 RowSideColors=NULL,
 ###
 ColSideColors=NULL,
 ###
 cexCol=1,
 cexRow=1,
 ccol=''
){
  
  
  
  require(gplots);
  genes<-intersect(rownames(expDat), genes);
  
  if(length(ccol)==4){
    ccol<-colorpanel(ccol$n, ccol$low,ccol$mid, ccol$high);
  }
  else{
    ccol<-bluered(100);
  }
  if(is.null(RowSideColors)){
    RowSideColors<-rep('white', length(genes));
  }
  if(is.null(ColSideColors)){
    ColSideColors<-rep('white', ncol(expDat));
  }
  if(!clusterR & ! clusterC){
    dendrogram='none';
  }
  if(clusterR & clusterC){
    dendrogram='both';
  }
  if(clusterR & !clusterC){
    dendrogram='row';
  }
  if(!clusterR & clusterC){
    dendrogram='column';
  }
  if(!big){
    heatmap.2(expDat[genes,],
              #col=bluered(100),
              col=ccol,
              scale=scale,
              trace='none',
              key=T,
              dendrogram=dendrogram,
              Rowv=clusterR,
              Colv=clusterC,
              density.info='none',
              margin=margin,
              colsep=seq(ncol(expDat)),
              rowsep=seq(length(genes)),
              sepcol='white',
              sepwidth=c(0.001,0.00),
              main=main,
              dist=dist,
              RowSideColors=RowSideColors,
              ColSideColors=ColSideColors,
              cexCol=cexCol,
              cexRow=cexRow);
  }
  else{
    heatmap.2(expDat[genes,],col=ccol, scale=scale, trace='none', key=T,dendrogram=dendrogram,Rowv=clusterR,Colv=clusterC,density.info='none',margin=margin,main=main,dist=dist,labRow='',labCol='',RowSideColors=RowSideColors,ColSideColors=ColSideColors,cexCol=cexCol,cexRow=cexRow);
  }
}



cn_plotnis<-function
### plot Network influence scores
(scoresDF,
 ### a numeric matrix of tfscores, rownames=genes
 main=NULL, 
 ### title
 limit=50
 ### limit output to top 50 TFs
 ){

  xmeans<-apply(scoresDF, 2, mean);
  worst<-which.min(xmeans);
  #tfs<-rownames(scoresDF);
  tfs<-rownames(scoresDF)[order(scoresDF[,worst], decreasing=F)];
  #tfs<-tfs[order(xmeans, decreasing=F)];
  scoresDF<-scoresDF[tfs,];
  if(nrow(scoresDF)>limit){
    scoresDF<-scoresDF[1:limit,];
  }
  
   scale='none';
   clusterR=FALSE;
   clusterC=FALSE;
   margin=c(12,6);
   cexR=.75;
   cexC=.5;
   aDist=dist;
   col=colorpanel(99, "blue","white","red");
   breaks=seq(from=-5, to=5,length.out=100);
  heatmap.2(scoresDF,
              col=col,
              breaks=breaks,
              scale=scale,
              trace='none',
              key=T,
              Rowv=clusterR,
              Colv=clusterC,
              density.info='none',
              cexRow=cexR,
              cexCol=cexC,
              margin=margin,
              colsep=seq(ncol(scoresDF)),
              rowsep=seq(nrow(scoresDF)),
              sepcol='white',
              sepwidth=c(0.001,0.001),
              dist=aDist);  
    
  ##scores1a<-tsco1[tsco1$tfScore<0 ,]#& tsco1$totalScore>0,]
  #  ggplot(data=scoresDF, aes(x=tf, y=totalScore))+ geom_bar(stat="identity", width=.8, fill=color) +
  #    theme_bw() + coord_flip() + xlab("") +  theme(axis.text.y = element_text(size=4)) + ylab("Network influence score") + ggtitle(main);
  # barplot of NIS
  #mp_hmVars(scoresDF, rownames(scoresDF), scale='none', clusterR=F, clusterC=F, cexR=.75)
  # heatmap of TF scores
}



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
  
  ddq<-cnRes[['dLevelQuery']];
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


