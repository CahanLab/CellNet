
# Patrick Cahan (C) 2017
# patrick.cahan@gmail.com

# find most variable genes
#' @export
findVarGenes<-function
(expNorm,
  geneStats,
  zThresh=2, 
  meanType="overall_mean"){
  allGenes<-rownames(expNorm)
  sg<-binGenesAlpha(geneStats)
  zscs<-rep(0, nrow(sg));
  names(zscs)<-rownames(sg);
  bbins<-unique(sg$bin);
  for(bbin in bbins){
    xx<-sg[sg$bin==bbin,];
    tmpZ<-scale(xx$fano);
    zscs[ rownames(xx) ]<-tmpZ[,1];
  }

  zByAlpha<-names(which(zscs>zThresh))
  
  # by overall mean
  sg<-binGenes(geneStats, meanType=meanType)
  zscsM<-rep(0, nrow(sg));
  names(zscsM)<-rownames(sg);
  bbins<-unique(sg$bin);
  for(bbin in bbins){
    xx<-sg[sg$bin==bbin,];
    tmpZ<-scale(xx$fano);
    zscsM[ rownames(xx) ]<-tmpZ[,1];
  }
  zByMean<-names(which(zscsM>zThresh))
  union(zByMean, zByAlpha)
}


#' weighted subtraction from mapped reades, applied to all
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
weighted_down<-function
(expRaw,
 total
 ){
    expCountDnW<-apply(expRaw, 2, downSampleW, total)
    #log(1+expCountDnW)
    expCountDnW
  }

#' @export
trans_prop<-function 
(expDat,
 xFact=1e5
){
  ans<-matrix(0, nrow=nrow(expDat), ncol=ncol(expDat));
  for(i in seq(ncol(expDat))){
    ans[,i]<-expDat[,i]/sum(expDat[,i]);    
  }
  ans<-ans*xFact;
  colnames(ans)<-colnames(expDat);
  rownames(ans)<-rownames(expDat);
  log(1+ans)
}


#' @export
plot_tsne<-function
(sampTab,
 tsneRes,
 cName="group"
 ){

  datTab<-tsneRes[rownames(sampTab),]
  colnames(datTab)[1:2]<-c("dim.1", "dim.2")

  xres<-cbind(sampTab, datTab)
  xi<-which(colnames(xres)==cName)
  colnames(xres)[xi]<-"groupX"

  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 12,name = "Paired")))(length(unique(xres$groupX)))
  ggplot(xres, aes(x=dim.1, y=dim.2, colour=groupX) ) + geom_point(pch=19, alpha=3/4, size=1) + theme_bw() + scale_colour_manual(values=ColorRamp) #+ facet_wrap( ~ k, nrow=3)
}




#' @export
to_tsne<-function
(datMat, #cols are vars, rows are cells/samples
 perplexity=30,
 theta=0.30){
  tres<-Rtsne(datMat, pca=FALSE, perplexity=perplexity, theta=theta)
  xres<-tres$Y
  colnames(xres)<-c("TSNE.1", "TSNE.2")
  rownames(xres)<-rownames(datMat)
  xres
}


#' @export
sc_statTab<-function# make a gene stats table (i.e. alpha, mu, etc)
(expDat, # expression matrix
 dThresh=0 # threshold for detection
 ){
  
  statTab<-data.frame()
  muAll<-sc_compMu(expDat, threshold=dThresh);
  alphaAll<-sc_compAlpha(expDat,threshold=dThresh);
  meanAll<-apply(expDat, 1, mean);
  covAll<-apply(expDat, 1, sc_cov);
  fanoAll<-apply(expDat,1, sc_fano);
  maxAll<-apply(expDat, 1, max);
  sdAll<-apply(expDat, 1, sd);
  
  statTabAll<-data.frame(gene=rownames(expDat), mu=muAll, alpha=alphaAll, overall_mean=meanAll, cov=covAll, fano=fanoAll, max_val=maxAll, sd=sdAll)
  statTabAll;
}


# compute alpha given detection threshold
#' @export
sc_compAlpha<-function 
(expMat,
 threshold=0,
  pseudo=FALSE){
  
  indexFunction<-function(vector, threshold){
    names(which(vector>threshold));
  }
  
  indexes<-apply(expMat, 1, indexFunction, threshold);
  alphas<-unlist(lapply(indexes, length));
  ans<-alphas/ncol(expMat)
  if(pseudo){
    ans<-(alphas+1)/(ncol(expMat)+1)
  }
  ans
}

# compute Mu given threshold
#' @export
sc_compMu<-function 
(expMat,
 threshold=0){
  
  afunct<-function(vector, threshold){
    mean( vector[which(vector>threshold)] );
  } 
  
  mus<-unlist(apply(expMat, 1, afunct, threshold))
  mus[is.na(mus)]<-0;
  mus;
}

# replavce NAs with 0
repNA<-function
(vector){
  vector[which(is.na(vector))]<-0;
  vector;
}

#compute fano factor on vector
sc_fano<-function
(vector){
  var(vector)/mean(vector);
}

# compute coeef of variation on vector
sc_cov<-function
(vector){
  sd(vector)/mean(vector);
}


#' @export
tsneMult<-function# facet tsne plot by gene
(tsneDat, # cols TSNE.1, TNSE.2, and genes
 genesToPlot, # genes to plot
 colorPal="BuPu",
 revCol=TRUE
 ){
  require(tidyr)
  tsne<-tsneDat[,c("TSNE.1", "TSNE.2", genesToPlot)]
  tsneLong<-gather_(tsne, "gene", "expression", genesToPlot)

  if(revCol){
    ColorRamp <- rev(colorRampPalette(rev(brewer.pal(n = 7,name = colorPal)))(100))[10:100]
  }
  else{
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = colorPal)))(100)
  }
  ggplot(tsneLong, aes(x=TSNE.1, y=TSNE.2, colour=expression) ) + 
  geom_point(pch=19, alpha=2/4, size=.25) + 
  theme_bw() + 
  scale_colour_gradientn(colours=ColorRamp) +
  facet_wrap( ~ gene)
 
}

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
# See: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' @export
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

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



### bins genes into x groups based on overallmean
binGenesAlpha<-function
(geneStats,
 nbins=20){
  max<-max(geneStats$alpha);
  min<-min(geneStats$alpha);

  binGroup<-rep(nbins, length=nrow(geneStats));
  names(binGroup)<-rownames(geneStats);
  rrange<-max-min;
  inc<-rrange/nbins
  borders<-seq(inc, max, by=inc)
  for(i in length(borders):1){
    xnames<-rownames(geneStats[which(geneStats$alpha<=borders[i]),]);
    binGroup[xnames]<-i;
  }
  cbind(geneStats, bin=binGroup);
}

### bins genes into x groups based on overallmean
binGenes<-function
(geneStats,
 nbins=20,
 meanType="overall_mean"){

##  max<-max(geneStats$overall_mean);
##  min<-min(geneStats$overall_mean);

  max<-max(geneStats[,meanType])
  min<-min(geneStats[,meanType])

  binGroup<-rep(nbins, length=nrow(geneStats));
  names(binGroup)<-rownames(geneStats);
  rrange<-max-min;
  inc<-rrange/nbins
  borders<-seq(inc, max, by=inc)
  for(i in length(borders):1){
    xnames<-rownames(geneStats[which(geneStats$overall_mean<=borders[i]),]);
    binGroup[xnames]<-i;
  }
  cbind(geneStats, bin=binGroup);
}


