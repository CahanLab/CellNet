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
  for(ctt in ctts){
    stTmp<-sampTab[sampTab[,dLevel]==ctt,];
    stTrain<-rbind(stTrain, .samp_for_class(stTmp,prop));
  }
  stTrain;
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




