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
  
  
  geneLists<-ctGRNs[['ctGRNs']][['geneLists']];

  commonTypes<-intersect(stTrain[,dLevel], names(geneLists));

  ctGRNs[['ctGRNs']][['geneLists']]<-ctGRNs[['ctGRNs']][['geneLists']][commonTypes];
  ctGRNs[['ctGRNs']][['graphLists']]<-ctGRNs[['ctGRNs']][['graphLists']][commonTypes];
  ctGRNs[['ctGRNs']][['tfTargets']]<-ctGRNs[['ctGRNs']][['tfTargets']][commonTypes];
  geneLists<-ctGRNs[['ctGRNs']][['geneLists']];
  
  # Make classifiers
  cat("Making classifiers (this can take awhile) ...\n");
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


