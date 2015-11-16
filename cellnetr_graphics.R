# CellNet
# (C) Patrick Cahan 2012-2014

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

cn_barplot_grnSing<-function ### wrapper to barplot secific GRN
(cnObj,
 ### result of analyzing query data with CN
 cnProc,
 ### result of creating a cellnet processor
 snName,
 ### name of subnet to plot establishment level 
 ctrlSamps,
 ### names of samples in training data
 bOrder
 ### order of bars
){
  
  qScores<-cnObj[['normScoresQuery']];
  ctrlScores<-cnProc[['trainingScores']]
  
  .cn_barplot_grnSing(qScores, ctrlScores, cnObj[['stQuery']], cnObj[['dLevelQuery']], snName, ctrlSamps, bOrder);
}


.cn_barplot_grnSing<-function### barplot this specific GRN
(qScores, ### queryScores
 ctrlScores, #### control scores
 stQuery,
 dLevelQ,
 snName,### name of subnet to plot establishment level 
 ctrSamps,### names of samples in training data
 bOrder  ### order of bars
 ){
    
  # convert into a data.frame
  aa<-cn_extract_SN_DF(qScores, stQuery,dLevelQ, rnames=snName);
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


plot_grn_summary<-function# plot basics stats of grnTable produced by cn_grnDoRock
(grnRes### result of running cn_grnDoRock
  ){
  x<-grnRes$grnStuff$grnTable;
  nTargetsPerTF<-table(x$TF);
  ntp<-data.frame(TF=names(nTargetsPerTF), count=as.vector(nTargetsPerTF));
  nRegPerTarg<-table(x$TG)
  npp<-data.frame(TF=names(nRegPerTarg), count=as.vector(nRegPerTarg));

  x1<-ggplot(ntp, aes(x=count)) + geom_histogram(colour="black", fill="white") + 
    ggtitle("Number of targets per TF") + labs(x="Number of targets");

  x2<-ggplot(npp, aes(x=count)) + geom_histogram(colour="black", fill="white") + 
    ggtitle("Number of regulators per target") + labs(x="Number of regulators");

  ###multiplot(x1, x2,cols=2);
  list(x1, x2);
}

plot_commList_summary<-function# plot basics stats of graphList produced by cn_grnDoRock
(grnRes### result of running cn_grnDoRock
  ){
    x<-grnRes$ctGRNs$graphLists;
    funsize<-function(agraph, type="nodes"){
      if(type=="nodes"){
        ans<-length(V(agraph));
      }
      else{
        ans<-length(E(agraph));
      }
      ans;
    }

    nnodes<-unlist(lapply(x, funsize, type='nodes'));
    df1<-data.frame(sn_name=names(x), size=nnodes);
    xplot1<-ggplot(df1, aes(x=size)) + geom_histogram(colour="black", fill="white") + 
    ggtitle("Node distribution") + labs(x="Number of nodes") + labs(y='Number of subnets');

    nedges<-unlist(lapply(x, funsize, type='edges'));
    df2<-data.frame(sn_name=names(x), size=nedges);
    xplot2<-ggplot(df2, aes(x=size)) + geom_histogram(colour="black", fill="white") + 
    ggtitle("Edge distribution") + labs(x="Number of edges") + labs(y='Number of subnets');
    ####multiplot(xplot1, xplot2, cols=2);
    list(xplot1, xplot2);
  }

plot_ctSN_general_summary<-function# plot basics stats of ct communities produced by cn_grnDoRock
(grnRes### result of running cn_grnDoRock
  ){
    ctGRNs<-grnRes$ctGRNs$graphLists; # list of ct -> one graph

    nnodes<-unlist(lapply(ctGRNs, vcount));
    df1<-data.frame(sn_name=names(ctGRNs), size=nnodes);
    xplot1<-ggplot(df1, aes(x=sn_name, y=size, width=.7)) + geom_bar(stat="identity") + #, colour="black", fill="white") + 
    ggtitle("GRN size - nodes") + labs(y="Number of nodes") + labs(x='') + coord_flip() + theme_bw() + theme(axis.text.y = element_text(size=5))

    nedges<-unlist(lapply(ctGRNs, ecount));
    df2<-data.frame(sn_name=names(ctGRNs), size=nedges);
    xplot2<-ggplot(df2, aes(x=sn_name, y=size, width=.7)) + geom_bar(stat="identity") + #, colour="black", fill="white") + 
    ggtitle("GRN size - edges") + labs(y="Number of edges") + labs(x='') + coord_flip() + theme_bw() + theme(axis.text.y = element_text(size=5))
    ###multiplot(xplot1, xplot2, cols=1);
    list(xplot1, xplot2);
}

plot_ctSN_subnet_summary<-function# plot basics stats of subnet ct communities produced by cn_grnDoRock
(grnRes### result of running cn_grnDoRock
  ){

    nodesize<-function(agraph, ntype='Target'){
      length(V(agraph)[V(agraph)$type==ntype]);
    }
    

    ctGRNs<-grnRes$ctGRNs$subnets$graphs; # list of ct -> one graph
    ctts<-names(ctGRNs);
    ntargs<-vector();
    nregs<-vector();
    snnames<-vector();
    df1<-data.frame();
    for(ctt in ctts){
      ntargs<-append(ntargs, unlist(lapply(ctGRNs[[ctt]], nodesize, ntype='Target')));
      nregs<-append(nregs, unlist(lapply(ctGRNs[[ctt]], nodesize, ntype='Regulator')));
      snnames<-append(snnames, names(ctGRNs[[ctt]]));
    }
    df1<-data.frame(sn_name=snnames, size=ntargs, type="Target")
    df1<-rbind(df1, data.frame(sn_name=snnames, size=nregs, type="Regulator"));

    xplot1<-ggplot(df1, aes(x=sn_name, y=size, width=.7, fill=type)) + geom_bar(stat="identity") + #, colour="black", fill="white") 
      ggtitle("Sub-network size - nodes") + labs(y="Number of genes") + 
      labs(x='') + coord_flip() + theme_bw() + 
      theme(axis.text.y = element_text(size=5))

    
    if(FALSE){
    nnodes<-unlist(lapply(ctGRNs, vcount));
    df1<-data.frame(sn_name=names(ctGRNs), size=nnodes);
    
    nedges<-unlist(lapply(ctGRNs, ecount));
    df2<-data.frame(sn_name=names(ctGRNs), size=nedges);
    xplot2<-ggplot(df2, aes(x=sn_name, y=size, width=.7)) + geom_bar(stat="identity") + #, colour="black", fill="white") + 
    ggtitle("GRN size - edges") + labs(y="Number of edges") + labs(x='') + coord_flip() + theme_bw()
    multiplot(xplot1, xplot2, cols=1);
  }
  xplot1
}



# FROM::: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

grn_report<-function# make a one page pdf of GRN results
(ctGRNs ### result of running cn_grnDoRock
){  
  a1<-plot_grn_summary(ctGRNs);
  a2<-plot_commList_summary(ctGRNs);
  a3<-plot_ctSN_general_summary(ctGRNs)
  multiplot(a1[[1]], a2[[1]], a3[[1]], a1[[2]], a2[[2]], a3[[2]],cols=2);
}

plot_class_rocs<-function# plot results of cn_classifications
(assessed # result of runnung cn_classAssess
  ){
  ctts<-names(assessed);
  df<-data.frame();
  for(ctt in ctts){
    tmp<-assessed[[ctt]];
    tmp<-cbind(tmp, ctype=ctt);
    df<-rbind(df, tmp);
  }

  rocsAll<-transform(df, TP = as.numeric(as.character(TP)), 
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

  precs<-precfunc(rocsAll)
  sens<-sensfunc(rocsAll)
  rocsAll2<-cbind(rocsAll, data.frame(recall=sens, precision=precs));

  ggplot(data=rocsAll2, aes(x=as.numeric(as.vector(recall)), y=as.numeric(as.vector(precision)))) + geom_point(size = 1) + 
  theme_bw() + xlab("Recall") + ylab("Precision") + facet_wrap( ~ ctype, ncol=4) +
  theme(axis.text = element_text(size=5)) + ggtitle("Classification performance")
}


