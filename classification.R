# CellNet
# (C) Patrick Cahan 2012-2014

# making and applying classifiers

# 
# APPLYING
#
cn_classify<-function
### run binary classifiers on expression data
(classList, 
 ### result of running CN3_makeClassifier
 expDat, 
 ### expression data 
 cttComms 
 ### list of tct->list(grn)
){
  ctts<-names(classList);
  ans<-matrix(0, nrow=length(ctts), ncol=ncol(expDat));
  rownames(ans)<-ctts;
  for(ctt in ctts){
    #cat(ctt,"\n")
    ##    xgenes<-cttComms[[ctt]][[1]];
    xgenes<-cttComms[[ctt]];
    myProbs<-predict(classList[[ctt]], t(expDat[xgenes,]), type='prob');
    ans[ctt,]<-myProbs[,ctt];
  }
  colnames(ans)<-colnames(expDat);
  ans;
  # classification matrix
}


#
# MAKING CLASSIFIERS
#

samp_for_class<-function# return a sampTab for training a test classifier
(sampTab,
  prop=0.5,
  dLevel="description1"){

  ctts<-unique(as.vector(sampTab[,dLevel]));
  stTrain<-data.frame();
  stVal<-data.frame();
  for(ctt in ctts){
    stTmp<-sampTab[sampTab[,dLevel]==ctt,];
    stTrain<-rbind(stTrain, .samp_for_class(stTmp,prop));

    idsval<-setdiff( rownames(stTmp), rownames(stTrain) );
    stVal<-rbind(stVal, stTmp[idsval,]);
  }
  list(stTrain=stTrain, stVal=stVal);
}

.samp_for_class<-function
(sampTab,
  prop=0.5){

  expIDcounts<-sort(table(sampTab$exp_id));
  expIDs<-names(expIDcounts);
  total<-sum(expIDcounts);

  runningTotal<-0;
  i<-1;
  xi<-floor( length(expIDcounts)/2 );
  while(i<=length(expIDcounts)){
    runningTotal<-sum(expIDcounts[1:i]);
    if(runningTotal/total > prop){
      xi<- i-1;
      break;
    }
    i<-i+1;
  }
  expIDs<-expIDs[1:xi];
  
  stTrain<-data.frame();
  for(expID in expIDs){
    stTrain<-rbind(stTrain, sampTab[sampTab$exp_id==expID,]);
  }
  stTrain;
}

cn_makeRFs<-function# Make one Random Forest classification per TCT
(expTrain,
 stTrain,
 geneLists, # list of tct=>gene vector
 dLevel='description1'){
  ans<-list();
  cnames<-names(geneLists);
  for(cname in cnames){
    cat(cname,"\n");
    xgenes<-intersect(rownames(expTrain), geneLists[[cname]]);  
    ans[[cname]]<-cn_makeClassifier(expTrain[xgenes,],cname,stTrain, dLevel=dLevel);
  }
  ans;
}

cn_makeClassifier<-function # make random forest classifier
(trainScores, # a matrix of values
 ctt, # what cell type/tissues is this?
 stTrain,
 dLevel='description1'){
  resp<-as.vector(stTrain[,dLevel]);
  xi<-which(resp!=ctt);
  resp[xi]<-'other';  
  myClass<-randomForest(t(trainScores), as.factor(resp), ntree=2e3);
  myClass;
}

# Assessing


cn_classAssess<-function# make ROCs for each classifier
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
    tmpROC <- cn_eval(ct_scores[xname,],
                           stVal,
                           classLevels,
                           xname,threshs=seq(0,1, by=resolution));
    allROCs[[xname]]<-tmpROC;
    i<-1+i;
  }
  allROCs;
}

cn_findSensAt<-function # return the sens at max FPR)<5%
(rocRes,
 at=0.05){
  xi<-which(rocRes[,"FPR"]<=at);
  tmp<-rocRes[xi,];
  cat(nrow(tmp),"\n")
  as.vector(tmp[which.max(tmp[,"FPR"]),'TPR']);
}

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

cn_eval<-function# return a data frame of the number of TP, FN, FP, and TN, and pval cutoff
(vect, # named vector
 sampTab,
 dLevel, # description level)
 classification,
 threshs=seq(0,1,by=0.05) # pval cutoffs
){
  ans<-matrix(0,nrow=length(threshs), ncol=7);
  for(i in seq(length(threshs))){
    thresh<-threshs[i];
    ans[i,1:4]<-cn_clPerf(vect, sampTab, dLevel, classification, thresh);
  }
  ans[,5]<-threshs;
  colnames(ans)<-c("TP", "FN", "FP", "TN", "thresh","FPR", "TPR");
  TPR<-ans[,'TP']/(ans[,'TP']+ans[,'FN']);
  FPR<-ans[,'FP']/(ans[,'TN']+ans[,'FP']);
  ans[,'TPR']<-TPR;
  ans[,'FPR']<-FPR;
  ans;
}


cn_clPerf<-function # assumes rownames(sampTab) == sampTab identifier used as colname for vect
(vect,
 sampTab,
 dLevel,
 classification, # actual classification
 thresh){
  TP<-0;
  FN<-0;
  FP<-0;
  TN<-0;
  sampIDs<-names(vect);  
  classes<-as.vector(sampTab[sampIDs,dLevel]);
  
  actualPos<-as.vector(sampTab[sampTab[,dLevel]==classification,]$sample_id);#which(classes==classification));
  actualNeg<-setdiff(sampIDs, actualPos);
  
  calledPos<-names(which(vect>thresh));
  calledNeg<-names(which(vect<=thresh));
  
  TP <- length(intersect(actualPos, calledPos));
  FP <- length(intersect(actualNeg, calledPos));
  FN <- length(actualPos)-TP;
  TN <- length(actualNeg)-FP;
  c(TP, FN, FP, TN);  
}


