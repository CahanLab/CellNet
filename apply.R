# CellNet
# (C) Patrick Cahan 2012-2016

cn_apply<-function
### run CellNet on expQuery
(expQuery,
 ### expression matrix
 stQuery,
 ### sample table
 cnProc,
 ### cnProc object
 dLevelQuery="description1"
 ### stQuery column name
){
  

  # make sure stTrain has a sample_id col, and rownames are equivalent to it
  if(!any(colnames(stQuery)=='sample_id')){
    stQuery<-cbind(stQuery, sample_id=colnames(expQuery));
    rownames(stQuery)<-as.vector(stQuery[,"sample_id"])
  }


  # 08-07-14, set default to sample_name so that individual samples are plotted
  
  if(length(dLevelQuery)==0 ){
    dLevelQuery<-"sample_name";
  }
  else{
    xx<-which(colnames(stQuery)==dLevelQuery);
    if(length(xx)==0){
      dLevelQuery<-"sample_name";
    }
  }
  ctrlScores<-cnProc[['trainingScores']];   # subnet est scores in control  -- normalized, df
  normVals<-cnProc[['normVals']];              # average subnet est scores in control samples,  (list of ctt->subnet ave)
  minVals<-cnProc[['minVals']];             # min raw vals of grn establishments to shift by
  rawCtrlScores<-cnProc[['raw_scores']];
  tVals<-cnProc[['tVals']];
  subNets<-cnProc[['grnList']];
  classList<-cnProc[['classList']];
  classWeight<-cnProc[['classWeight']];
  exprWeight<-cnProc[['exprWeight']];
  
  # score the query data
  #cat("Scoring query data...\n")
  scoresQuery<-cn_score(expQuery, subNets, tVals, minVals=minVals, classList=classList, classWeight=classWeight, exprWeight=exprWeight);
  
  # normalize query scores
  #cat("normalizing grn scores\n");
  normScoresQuery<-cn_normalizeScores(normVals, scoresQuery, rownames(scoresQuery));
  
  # classify the query data
  #cat("Classifying query data ...\n");
  myClass<-cn_classify(classList,expQuery,subNets);
  colnames(myClass)<-as.vector(stQuery[,dLevelQuery])
  
  # it's really useful to have the averaged training data, too
  meanExpQuery<-as.matrix(GEP_makeMean(expQuery,as.vector(stQuery[,dLevelQuery]),type='mean'));
  
  ans<-list(expQuery=expQuery,
            stQuery=stQuery,
            dLevelQuery=dLevelQuery,            
            queryScores=scoresQuery,
            classRes=myClass,
            normScoresQuery=normScoresQuery,
            expMeanQuery=meanExpQuery);
  class(ans)<-"cnRes";
  ans;  
  # cnRes object
}



cn_outputRes<-function
### write out grn scores, classification scores, normalized data to csv, then tar and  zip
(cnRes,
 ### cnRes object
 tfScores,
 ### NIS
 prefix
 ### filename prefix
 ){
  
  # write out classRes
  classRes<-cnRes[['classRes']];
  fname1<-paste(prefix, "_classRes.csv", sep='');
  write.csv(classRes, file=fname1);
  
  # write out normalized GRN statuses
  grnScores<-cnRes[['normScoresQuery']];
  fname2<-paste(prefix, "_grnScores.csv", sep='');
  write.csv(grnScores, file=fname2);
  
  # write out normalized query expression data
  expDat<-cnRes[['expQuery']];
  fname3<-paste(prefix, "_expNorm.csv", sep='');
  write.csv(expDat, file=fname3);
  
  # write out one data file of TF scores
  tfStab<-.cn_makeTFtable(tfScores);
  fname4<-paste(prefix, "_NIS.csv", sep='');
  write.csv(tfStab, file=fname4, row.names=F);
  
  fnameOut<-paste(prefix,"_data.tar", sep='')
  cmd<-paste("tar -c ", fname1, " ",fname2," ", fname3," > ", fnameOut, sep='');
  system(cmd);
  cmd<-paste("gzip -f ", fnameOut, sep='');
  system(cmd);
  fnameOut;
  ### name of output file
}

.cn_makeTFtable<-function
### convert a tf nis list to a DF
(tfScores){
  allTFs<-data.frame();
  grnNames<-names(tfScores);
  for(grnName in grnNames){
    x<-tfScores[[grnName]];
    genes<-rownames(x);
    #rownames(x)<-'';
    x2<-data.frame(reg=genes, grn=rep(grnName, length(genes)));
    x2<-cbind(x2, x);
    allTFs<-rbind(allTFs, x2);
  }
  allTFs
}
