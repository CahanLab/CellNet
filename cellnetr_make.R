# CellNet
# (C) Patrick Cahan 2012-2014
# functions to set up for sample classification

# To train a CellNet processor:
# cnProc<-CN3_make_processor(expTrain, stTrain, igGRNs, TCT_order, dLevel) 
# 
# To analyze a query data set:
# cnRes<-CN3_analyze(cnProc, expQuery, stQuery, dLevel)
# 
# To produce standard output (classification HM, donor, target, aberant GRN establishments):
# CN3_stdOut(cnProc, cnRes, prefix)



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