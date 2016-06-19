# CellNet
# (C) Patrick Cahan 2012-2014

# make a CellNet object

cn_make_processorOnly<-function # train a CellNet object
(expTrain,
 stTrain,
 ctGRNs, # result of running cn_grnDoRock
 tctOrder=NULL,
 dLevel="description1",
 classWeight=TRUE, # weight GRN est by classification importance
 exprWeight=TRUE # weight GRN est by gene expression
 ){
  
  geneLists<-ctGRNs[['ctGRNs']][['geneLists']];
  commonTypes<-intersect(stTrain[,dLevel], names(geneLists));
  geneLists<-geneLists[commonTypes];

  # Make classifiers
  cat("Making classifiers (this can take awhile) ...\n");
#  classList<-make_classifiers(ctGRNs, expTrain, stTrain, dLevel);
  classList<-cn_makeRFs(expTrain, stTrain, geneLists, dLevel=dLevel);
  cat("Done making classifiers :)\n");
  
  trainNorm<-cn_trainNorm(expTrain, stTrain, subNets=geneLists, classList=classList, dLevel=dLevel, classWeight=classWeight, exprWeight=exprWeight);
   
  cat("done training norm\n")
  meanExpTrain<-as.matrix(GEP_makeMean(expTrain,as.vector(stTrain[,dLevel]),type='mean'));
  cat("done GEP make mean\n")


  ans<-list(expTrain=expTrain,
            expMeanTrain=meanExpTrain,
            stTrain=stTrain,
            dLevelTrain=dLevel,
            ctGRNs=ctGRNs,
            grnList=geneLists,
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


cn_make_processor<-function # train a CellNet object
(expTrain,
 stTrain,
 ctGRNs, # result of running cn_grnDoRock
 expVal, # Data for Classifier validation
 stVal, # Data for Classifier validation
 tctOrder=NULL,
 dLevel="description1",
 classWeight=TRUE, # weight GRN est by classification importance
 exprWeight=TRUE # weight GRN est by gene expression
 ){
  
  geneLists<-ctGRNs[['ctGRNs']][['geneLists']];
  commonTypes<-intersect(stTrain[,dLevel], names(geneLists));
  geneLists<-geneLists[commonTypes];

  # Make classifiers
  cat("Making classifiers (this can take awhile) ...\n");
#  classList<-make_classifiers(ctGRNs, expTrain, stTrain, dLevel);
  classList<-cn_makeRFs(expTrain, stTrain, geneLists, dLevel=dLevel);
  cat("Done making classifiers :)\n");
  
  trainNorm<-cn_trainNorm(expTrain, stTrain, subNets=geneLists, classList=classList, dLevel=dLevel, classWeight=classWeight, exprWeight=exprWeight);
   
  cat("done training norm\n")
  meanExpTrain<-as.matrix(GEP_makeMean(expTrain,as.vector(stTrain[,dLevel]),type='mean'));
  cat("done GEP make mean\n")

  #Make Object for Classifier Assessment
  ###ansVal<-cn_classify(classList, expVal, ctGRNs[['ctGRNs']][['geneLists']]) 
  ansVal<-cn_classify(classList, expVal, geneLists) 
  cat("done classifying validation data\n");

  assessed<-cn_classAssess(ansVal, stVal, classLevels="description1", resolution=0.01) 
  cat("done assessing CN\n")  


  ans<-list(expTrain=expTrain,
            expMeanTrain=meanExpTrain,
            stTrain=stTrain,
            dLevelTrain=dLevel,
            ctGRNs=ctGRNs,
            grnList=geneLists,
            classList=classList,
            trainingScores=trainNorm[['trainingScores']],
            normVals=trainNorm[['normVals']],
            raw_scores=trainNorm[['raw_scores']],
            minVals=trainNorm[['minVals']],
            tVals=trainNorm[['tVals']],
            classWeight=classWeight,
            exprWeight=exprWeight,
            ClassifierAssessment = assessed);
  class(ans)<-"cnProc";
  ans;
}

make_classifiers<-function# wrappper to cn_makeRFs
(ctGRNs, ### result of running cn_grnDoRock
 expDat, ### training data,
 sampTab, ### sample table
 dLevel ### description level to indicate cell types
 ){
	geneLists<-ctGRNs[['ctGRNs']][['grnLists']];
 	cn_makeRFs(expDat,sampTab,geneLists, dLevel=dLevel);
}



