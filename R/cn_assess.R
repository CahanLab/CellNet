# CellNet
# (C) Patrick Cahan 2012-2016


#' Assess classifiers based on validation data
#'
#' Assess classifiers based on validation data
#' @param ct_scores matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
#' @param stVal sample table
#' @param classLevels column name of stVal to use as ground truth to assess classifiers
#' @param resolution increment at which to evalutate classification
#' @param dLevelSID column to indicate sample id
#'
#' @return list of data frames with threshold, sens, precision
#' @export 
cn_classAssess<-function# make ROCs for each classifier
(ct_scores,# matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
 stVal, # sampletable
 classLevels="description2",
 dLevelSID="sample_id",
 resolution=0.005 # increment at which to evalutate classification
){
  allROCs<-list();
  evalAll<-matrix(0, nrow=nrow(ct_scores),ncol=2);
  classifications<-rownames(ct_scores);
  rownames(stVal)<-as.vector(stVal[,dLevelSID]);
  i<-1;
  for(xname in classifications){
    classification<-classifications[i];
    tmpROC <- cn_eval(ct_scores[xname,],
                           stVal,
                           classLevels,
                           xname,threshs=seq(0,1, by=resolution), dLevelSID=dLevelSID);
    allROCs[[xname]]<-tmpROC;
    i<-1+i;
  }
  allROCs;
}

#' return the sens at given FPR
#'
#' return the sens at given FPR
#' @param rocRes result of running cn_eval
#' @param what is sens at this FPR?
#'
#' @return what is sens at selected FPR?
cn_findSensAt<-function
(rocRes,
 at=0.05){
  xi<-which(rocRes[,"FPR"]<=at);
  tmp<-rocRes[xi,];
  cat(nrow(tmp),"\n")
  as.vector(tmp[which.max(tmp[,"FPR"]),'TPR']);
}

#' run cn_clPerf across thresholds
#'
#' run cn_clPerf across thresholds
#' @param vect named vector
#' @param dLevel description level
#' @param classification classification matrix
#' @param threshs seq of pval cutoffs
#' @param dLevelSID column to indicate sample id
#'
#' @return return a data frame of the number of TP, FN, FP, and TN, and pval cutoff
cn_eval<-function# return a data frame of the number of TP, FN, FP, and TN, and pval cutoff
(vect, # named vector
 sampTab,
 dLevel, # description level)
 classification,
 threshs=seq(0,1,by=0.05), # pval cutoffs
 dLevelSID="sample_id"
){
  ans<-matrix(0,nrow=length(threshs), ncol=7);
  for(i in seq(length(threshs))){
    thresh<-threshs[i];
    ans[i,1:4]<-cn_clPerf(vect, sampTab, dLevel, classification, thresh, dLevelSID=dLevelSID);
  }
  ans[,5]<-threshs;
  colnames(ans)<-c("TP", "FN", "FP", "TN", "thresh","FPR", "TPR");
  TPR<-ans[,'TP']/(ans[,'TP']+ans[,'FN']);
  FPR<-ans[,'FP']/(ans[,'TN']+ans[,'FP']);
  ans[,'TPR']<-TPR;
  ans[,'FPR']<-FPR;
  ans;
}

#' determine performance of classification at given threshold
#'
#' determine performance of classification at given threshold
#' @param vect vector of values
#' @param sampTab sample table
#' @param dLevel colname
#' @param classification actual classification
#' @param thresh threshold above which to make a call
#' @param dLevelSID column to indicate sample id
#'
#' @return vector of TP FN FP TN 
cn_clPerf<-function # assumes rownames(sampTab) == sampTab identifier used as colname for vect
(vect,
 sampTab,
 dLevel,
 classification, # actual classification
 thresh,
 dLevelSID="sample_id"){
  TP<-0;
  FN<-0;
  FP<-0;
  TN<-0;
  sampIDs<-names(vect);  
  classes<-as.vector(sampTab[sampIDs,dLevel]);
  
  ###actualPos<-as.vector(sampTab[sampTab[,dLevel]==classification,]$sample_id);#which(classes==classification));
  actualPos<-as.vector(sampTab[sampTab[,dLevel]==classification,dLevelSID])
  actualNeg<-setdiff(sampIDs, actualPos);
  
  calledPos<-names(which(vect>thresh));
  calledNeg<-names(which(vect<=thresh));
  
  TP <- length(intersect(actualPos, calledPos));
  FP <- length(intersect(actualNeg, calledPos));
  FN <- length(actualPos)-TP;
  TN <- length(actualNeg)-FP;
  c(TP, FN, FP, TN);  
}

### GRN STUFF ###

#' assess performance of GRN predictions based on zscores. this function runs cn_calcPRs on all TFs in gol_standard
#'
#' Recall (Sensitivity): fraction of gold standard regulatory relationships detected by GRN. Precision (1-FPR): Proportion of called relationships that are in the gold standard
#' @param zscores matrix of zscores
#' @param goldStandard gold stanard gene list
#' @param universe genes tested
#' @param tRange threshold sequence
#' @param random boolean, whether to compute a random selection of calls, too, based on the number of genes actually predicted
#' @param randomIter numer of random selects to make, returns only the average
#' @param funType how to treat >1 TF, intersect or union
#'
#' @return List of performance results of form: data frame of Score (cutoff), Precision, Recall, Predictions (n), obs = real GRN, rand = list of random GRNs, one per 
#'
#' @examples
#' prList<-cn_assessGRN(zscores, gs, random=TRUE)
#'
#' @export
cn_assessGRN<-function
(zscores,
 goldStandard,
 tRange=seq(3,6,by=0.25),
 random=FALSE,
 randomIter=10,
 funType='intersect'
){
  prs<-list();
  testTFs<-names(goldStandard)
  testTFs<-intersect(testTFs, colnames(zscores))
  for(ttf in testTFs){
    prs[[ttf]]<-cn_calcPRs(zscores, ttf, goldStandard[[ttf]], rownames(zscores), tRange=tRange, random=random, randomIter=randomIter)
  }
  prs
}

#' compile results of assessing GRNs for multiple TFs
#'
#' result can be immediately plotted
#' @param prs result of running cn_assessGRN
#'
#' @return data frame
#'
#' @examples
#' prList<-cn_assessGRN(zscores, gs, random=TRUE)
#' assessedGRN<-cn_compileGRNassessments(prList)
#' ggplot(assessedGRN, aes(x=factor(Score), y=Fold)) +
#'  geom_boxplot(alpha=.75,colour="black", outlier.colour="#BDBDBD", outlier.size=.5, fill="#FF9999") +
#'  theme_bw() + ggtitle("AUPR Fold increase") + ylim(0,15) + ylab("AUPR fold improvement over random") + 
#'  xlab("GRN Z-score threshold")
#'
#' @export
cn_compileGRNassessments<-function
(prs){
  x3<-data.frame();
  nnames<-names(prs);
  for(nname in nnames){
    a<-cn_addAUPRs(prs[[nname]]);
    a2<-cn_summ_AUPRs(a, nname);
    x3<-rbind(x3, pool_auPRs(a2[[1]]));
  }
  x3;
}


#' assess performance of GRN predictions based on zscores
#'
#' Recall (Sensitivity): fraction of gold standard regulatory relationships detected by GRN. Precision (1-FPR): Proportion of called relationships that are in the gold standard
#' @param zscores matrix of zscores
#' @param tfs vector of tf(s) perturbed
#' @param goldStandard gold stanrd gene list
#' @param universe genes tested
#' @param tRange threshold sequence
#' @param random boolean, whether to compute a random selection of calls, too, based on the number of genes actually predicted
#' @param randomIter numer of random selects to make, returns only the average
#' @param funType how to treat >1 TF, intersect or union
#'
#' @return List of performance results of form: data frame of Score (cutoff), Precision, Recall, Predictions (n), obs = real GRN, rand = list of random GRNs
#'
#' @examples
#' grnPR<-cn_calcPRs(zscs, "Sox2", aGoldStandard[["Sox2"]], allGenes, tRange=seq(0,8, by=.5), random=TRUE, randomIter=3);
#'
cn_calcPRs<-function
(zscores,
 tfs,
 goldStandard,
 universe=NULL,
 tRange=seq(3,6,by=0.25),
 random=FALSE,
 randomIter=10,
 funType='intersect'
){
  if(!is.null(universe)){
    zscores<-zscores[universe,];
    goldStandard<-intersect(universe, goldStandard);
  }
  
  allGenes<-rownames(zscores);
  ans<-matrix(0, nrow=length(tRange), ncol=4);
  randList<-list();
  # setup list to store performance of random GRNs
  if(random){
    for(xi in seq(randomIter)){
      randList[[xi]]<-matrix(0, nrow=length(tRange), ncol=4);
      colnames(randList[[xi]])<-c("Score", "Precision", "Recall", "Predictions");
    }
  }  
  
  x<-zscores[,tfs];
  numTFs<-length(tfs);
  for(i in seq(length(tRange))){
    threshold<-tRange[i];
    if(length(tfs)>1){
      tmpPos<-vector();
      for(ai in seq(length(tfs))){            
        tmpPos<-append(tmpPos, rownames(zscores[ which(x[,ai]>threshold) ,]));
      }
      if(funType=='intersect'){
        aTab<-table(tmpPos);
        positives<-names(which(aTab==numTFs));
      }
      else{
        positives<-unique(tmpPos);
      }
    }
    else{
      positives<-rownames(zscores[ which(x>threshold) ,]);
    }
    nPos<-length(positives);
    TP<-intersect(goldStandard, positives);
    FN<-setdiff(goldStandard, positives);    
    precision <- length(TP)/nPos;
    recall <- length(TP)/(length(TP)+length(FN));
    ans[i,]<-c(threshold, precision, recall, nPos);
  }
  if(random){
    for(j in seq(randomIter)){
      # to randomly select TF, uncomment next line
      # xTF<-sample(colnames(zscores),1);
      # re-order zscores
      newX<-sample(zscores[,tfs],length(x));
      for(i in seq(length(tRange))){
        threshold<-tRange[i];
        positives<-rownames(zscores[ which(newX>threshold) ,]);          
        nPos<-length(positives);
        TP<-intersect(goldStandard, positives);
        FN<-setdiff(goldStandard, positives);    
        precision <- length(TP)/nPos;
        recall <- length(TP)/(length(TP)+length(FN));
        randList[[j]][i,]<-c(threshold, precision, recall, nPos);
      }
    }
  }
  colnames(ans)<-c("Score", "Precision", "Recall", "Predictions");
  
  list(obs=ans, rand=randList);
}

#' convert a result of cn_calcPRs to something easily plot-able
#'
#' convert a result of cn_calcPRs to something easily plot-able
#' @param prRes result of running cn_calcPRs
#'
#' @return something
cn_convertPR<-function
(prRes
){
  obs<-prRes[['obs']];
  obs<-cbind(obs, type=rep('obs', nrow(obs)));
  
  rrand<-prRes[['rand']];
  for(i in seq(length(rrand))){
    tmp<-rrand[[i]];
    tmp<-cbind(tmp, type=rep('random', nrow(tmp)));
    obs<-rbind(obs, tmp);
  }
  ans<-data.frame(obs);
  ans<-transform(ans, Precision = as.numeric(as.vector(Precision)));
  ans<-transform(ans, Recall = as.numeric(as.vector(Recall)));
  ans<-transform(ans, Predictions = as.numeric(as.vector(Predictions)));
  ans;  
}

#' utility function to compile cn_calcPRs
#'
#' utility function to compile cn_calcPRs
#' @param bigList alist
#'
#' @return something
cn_addAUPRs<-function### 
(bigList
){
  
  tmpX<-cn_computeAUCPR(bigList[[1]]);
  bigList[[1]]<- cbind(bigList[[1]], AUPR=tmpX);
  
  xList<-bigList[[2]];
  for(i in seq(length(xList))){
    tmpX<-cn_computeAUCPR(xList[[i]]);
    xList[[i]]<- cbind(xList[[i]], AUPR=tmpX);      
  }
  bigList[[2]]<-xList;
  bigList;
}

#' processes result of cn_addAUPRs: p-value of fold improvements,df for visualization
#'
#' processes result of cn_addAUPRs: p-value of fold improvements,df for visualization
#' @param bigList result of cn_addAUPRs
#' @param name name of analysis
#'
#' @return something
cn_summ_AUPRs<-function ### 
(bigList,
 name
){
  
  # compute p-vals at each zscores as 1 +  (number of times in which random AUPR >= observed AUPR) / 1 + number of random iterations 
  obs<-bigList[[1]];
  randoms<-bigList[[2]];
  numers<-rep(1, nrow(obs));
  denoms<-rep(1+length(randoms), nrow(obs));
  for(i in seq(nrow(obs))){
    auprs<-unlist(lapply(randoms, "[", i, "AUPR") );
    numers[i]<-1 + length(which(auprs>=obs[i,"AUPR"]));
  }
  pvals<-data.frame(Score=obs[,"Score"],Pval=numers/denoms, GS=rep(name, nrow(obs)));  
  
  ans<-data.frame();
  for(i in seq(length(randoms))){
    newF<-obs;
    newF<-cbind(newF, Fold=obs[,"AUPR"]/randoms[[i]][,"AUPR"]);
    ans<-rbind(ans, newF);
  }
  ans<-cbind(ans, gs=name);
  list(bigTable=ans, pvals=pvals);
}


#' compute AUPCR
#'
#' compute AUPCR
#' @param perfDF a df
#' @param precisionCol "Precision"
#' @param recallCol "Recall"
#' @param predCol "Predictions"
#'
#' @return auprs
cn_computeAUCPR<-function
(perfDF,
 precisionCol="Precision",
 recallCol="Recall",
 predCol="Predictions"
){
  
  ### Notes: starts at top left, and accumulates to max
  str<-(nrow(perfDF)-1);
  stp<-2;
  
  # sometimes you can still get NA in 
  areas<-rep(0, nrow(perfDF));
  
  pts<-seq(str,stp,-1);
  for(i in pts){
    a<-(i+1);
    cc<-(i-1);
    ptA<-c(perfDF[a,recallCol], perfDF[a,precisionCol]);
    ptB<-c(perfDF[i,recallCol], perfDF[i,precisionCol]);
    ptC<-c(perfDF[cc,recallCol], perfDF[cc,precisionCol]);
    tmpArea<-cn_rectArea(ptA, ptB, ptC);
    if(is.nan(tmpArea)){
      #cat(perfDF[i,]$Score,"\n");
      tmpArea<-0;
    }
    areas[i]<-areas[(i+1)]+tmpArea;
  }
  
  # far right point
  a<-2;
  b<-1;
  cc<-1;
  ptA<-c(perfDF[a,recallCol], perfDF[a,precisionCol]);
  ptB<-c(perfDF[i,recallCol], perfDF[i,precisionCol]);
  ptC<-c(perfDF[cc,recallCol], perfDF[cc,precisionCol]);
  areas[b]<-areas[2]+cn_rectArea(ptA, ptB, ptC);
  
  # far left point
  a<-nrow(perfDF);
  b<-nrow(perfDF);
  cc<-(b-1);
  ptA<-c(perfDF[a,recallCol], perfDF[a,precisionCol]);
  ptB<-c(perfDF[i,recallCol], perfDF[i,precisionCol]);
  ptC<-c(perfDF[cc,recallCol], perfDF[cc,precisionCol]);
  areas[b]<-cn_rectArea(ptA, ptB, ptC);
  areas;
}

#' compute area of rect given by 
#'
#' compute area of rect given by 
#' @param ptA point a
#' @param ptB point b
#' @param ptC point C
#'
#' @return area
cn_rectArea<-function
(ptA,
 ptB,
 ptC
){
  xRight <- ptB[1]+( (ptC[1]-ptB[1])/2);
  xLeft  <- ptA[1]+( (ptB[1]-ptA[1])/2);
  width<-abs(xRight - xLeft);
  width*ptB[2];
}

# it's really annoying that I used auto-select for the zscore cutoff for gs assessment because now ti si hard to compare across cell types/gs
# deal with this by pooling based on a range
# -3:10,inc=1 -3->-2,

#' pool auprs 
#'
#' pool auprs 
#' @param perfDF perfDF
#' @param rng seq -4 to 10
#'
#' @return df
pool_auPRs<-function###
(perfDF,
 ###
 rng=seq(-4,10,by=1)
 ###
){
  ans<-data.frame();
  nBins<-length(rng)
  for(i in seq( nBins-1)){
    val<-rng[i];
    #cat(val,"\n");
    nextVal<-rng[i+1];
    xi<-which( perfDF[,'Score']>=val & perfDF[,'Score']<nextVal);
    perfDF[xi,'Score']<-val;
  }
  # last one
  xi<-which(perfDF[,'Score']>=rng[nBins]);
  perfDF[xi,'Score']<-rng[nBins];
  perfDF;
}

######################################
#
# Unused
#
######################################

if(FALSE){
CN3_AUCs<-function# make ROCs for each classifier, and return AUC as data frame
(ct_scores,# matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids
 stVal, # sampletable
 classLevels="description2", #
 resolution=0.005 # increment at which to evalutate classification
){
  allROCs<-list();
  evalAll<-matrix(0, nrow=nrow(ct_scores),ncol=2);
  classifications<-rownames(ct_scores);
  rownames(stVal)<-as.vector(stVal$sample_id);
  i<-1;
  for(xname in classifications){
    classification<-classifications[i];
    tmpROC <- CN3_eval(ct_scores[xname,],
                           stVal,
                           classLevels,
                           xname,threshs=seq(0,1, by=resolution));
    allROCs[[xname]]<-tmpROC;
    evalAll[i,1]<-CellNet_auc(tmpROC);
    evalAll[i,2]<-classification;
    i<-i+1;
  }
  list(aucs=evalAll, rocs=allROCs);
}}

if(FALSE){
  CN3_findSensAt<-function # return the sens at max FPR)<5%
(rocRes,
 at=0.05){
  xi<-which(rocRes[,"FPR"]<=at);
  tmp<-rocRes[xi,];
  cat(nrow(tmp),"\n")
  as.vector(tmp[which.max(tmp[,"FPR"]),'TPR']);
}}

if(FALSE){
CN3_eval<-function# return a data frame of the number of TP, FN, FP, and TN, and pval cutoff
(vect, # named vector
 sampTab,
 dLevel, # description level)
 classification,
 threshs=seq(0,1,by=0.05) # pval cutoffs
){
  ans<-matrix(0,nrow=length(threshs), ncol=7);
  for(i in seq(length(threshs))){
    thresh<-threshs[i];
    ans[i,1:4]<-celnet_clPerf(vect, sampTab, dLevel, classification, thresh);
  }
  ans[,5]<-threshs;
  colnames(ans)<-c("TP", "FN", "FP", "TN", "thresh","FPR", "TPR");
  TPR<-ans[,'TP']/(ans[,'TP']+ans[,'FN']);
  FPR<-ans[,'FP']/(ans[,'TN']+ans[,'FP']);
  ans[,'TPR']<-TPR;
  ans[,'FPR']<-FPR;
  ans;
}}

if(FALSE){
  CellNet_auc<-function # basic AUC calc given ROC obj (as returned by celnet_eval)
(rocDF # cols: TP, FP, FN, TN, thresh, TPR, FPR
)
{
  # note: calculate as SUM( average(height (TPRx+1, TPRx) * width ( FPRx+1 - FPRx )
  llen<-nrow(rocDF);
  mySum<-0
  for(i in seq(llen-1)){
    myWidth<-rocDF[i,'FPR']-rocDF[(i+1),'FPR'];
    myHeight<- mean( c(rocDF[(i+1),'TPR'], rocDF[i,'TPR']));
    # cat(i,"\t", myHeight,"\t", myWidth,"\n");
    mySum<-mySum + (myHeight*myWidth);
  }
  as.vector(mySum);
}}

if(FALSE){
  celnet_addErr<-function # average random values and add error bars
(prRes){
  tmp<-celnet_convertPR(prRes);
  scores<-unique(tmp$Score);
  types<-unique(tmp$type);
  nper<-length(scores);
  
  ans<-data.frame();
  for(type in types){
    median_precision<-rep(0, nper);
    median_recall<-rep(0, nper);
    precision_stdvs<-rep(0, nper);
    recall_stdvs<-rep(0, nper);
    mean_predictions<-rep(0, nper);
    types<-rep(type, nper);
    tmpDF<-tmp[tmp$type==type,];
    for(i in seq(nper)){
      tmpX<-tmpDF[tmpDF$Score==scores[i],];
      median_precision[i]<-median(tmpX[,"Precision"]);
      median_recall[i]<-median(tmpX[,"Recall"]);
      mean_predictions[i]<-mean(tmpX[,"Predictions"]);
      x<-utils_SD(tmpX[,"Precision"]);
      if(is.na(x)){
        x<-0;
      }
      precision_stdvs[i]<-x;
      x<-utils_SD(tmpX[,"Recall"]);
      if(is.na(x)){
        x<-0;
      }
      recall_stdvs[i]<-x;
    }
    tmpAns<-data.frame(Score=scores, median_precision=median_precision,
                       median_recall=median_recall,
                       mean_predictions=mean_predictions,
                       precision_stdvs=precision_stdvs,
                       recall_stdvs=recall_stdvs,
                       type=types);
    ans<-rbind(ans, tmpAns);
  }
  ans;
}
}

if(FALSE){

cn_makePRtab<-function###
(aList
){
  ans<-data.frame();
  nnames<-names(aList);
  for(i in seq(length(aList))){
    wReal<-aList[[i]][[1]];
    tmpPRobs<-cn_computeAUCPR(wReal);
    tmpPRr<-cn_computeAUCPR(wReal, "randomPrecision", "randomRecall");
    tmpAns<-data.frame(name=rep(titles[i],2), type=c("CLR", "Random"), AUPR=c(tmpPRobs, tmpPRr));
    ans<-rbind(ans, tmpAns);
  }
  ans;
}
}

if(FALSE){
cn_makeRandNets<-function # make random grns as baseline comparison
(allGenes,
 cttComms){
  newComms<-cttComms;
  for(nname in names(newComms)){
    x<-newComms[[nname]];
    lens<-sapply(x, length);
    newNets<-list();
    xNames<-names(x);
    for(xName in xNames){    
      newNets[[xName]]<-sample(allGenes, lens[xName]);
    }
    newComms[[nname]]<-newNets;
  }
  newComms;  
}
}

