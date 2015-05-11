# CellNet
# (C) Patrick Cahan 2012-2014

# make a CellNet object

cn_make_processor<-function # train a CellNet object
(expTrain,
 stTrain,
 ctGRNs, # result of running cn_grnDoRock
 tctOrder=NULL,
 dLevel="description1",
 classWeight=TRUE, # weight GRN est by classification importance
 exprWeight=TRUE # weight GRN est by gene expression
 ){
  
  grnList<-ctGRNs$geneLists;
  
  # Make classifiers
  cat("Making classifiers (this can take awhile) ...\n");
#  classList<-make_classifiers(ctGRNs, expTrain, stTrain, dLevel);
  classList<-cn_makeRFs(expTrain, stTrain, grnList, dLevel=dLevel);
  cat("Done making classifiers :)\n");
  
  trainNorm<-cn_trainNorm(expTrain, stTrain, subNets=grnList, classList=classList, dLevel=dLevel, classWeight=classWeight, exprWeight=exprWeight);
  meanExpTrain<-as.matrix(GEP_makeMean(expTrain,as.vector(stTrain[,dLevel]),type='mean'));
  
  ans<-list(expTrain=expTrain,
            expMeanTrain=meanExpTrain,
            stTrain=stTrain,
            dLevelTrain=dLevel,
            ctGRNs=ctGRNs,
            grnList=grnList,
            classList=classList,
            trainingScores=trainNorm[['trainingScores']],
            normVals=trainNorm[['normVals']],
            raw_scores=trainNorm[['raw_scores']],
            minVals=trainNorm[['minVals']],
            tVals=trainNorm[['tVals']],
            classWeight=classWeight,
            exprWeight=exprWeight);
  class(ans)<-"cnProc";
  ans;
}



make_classifiers<-function# wrappper to cn_makeRFs
(ctGRNs, ### result of running cn_grnDoRock
 expDat, ### training data,
 sampTab, ### sample table
 dLevel ### description level to indicate cell types
 ){
	geneLists<-ctGRNs$ctGRNs$general$geneLists;
 	cn_makeRFs(expDat,sampTab,geneLists, dLevel=dLevel);
}

