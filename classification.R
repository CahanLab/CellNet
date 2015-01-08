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




