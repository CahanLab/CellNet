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


cn_normalizeScores<-function
### divide the query scores by the mean values in the training data.
(ctrlScores, 
 ### a list of subnet->mean value, all subnets
 queryScores, 
 ### a matrix, rownames = subnet names, all subnets
 subNets 
 ### a vector of subnets names to use
){
  
  ans<-matrix(0, nrow=length(subNets), ncol=ncol(queryScores));
  rownames(ans)<-subNets
  #subNets<-rownames(queryScores);
  for(subNet in subNets){
    cat(subNet,"\n")
    ans[subNet,]<- queryScores[subNet,] / ctrlScores[[subNet]];
  }
  colnames(ans)<-colnames(queryScores);
  ans;
  ### normalized expression matrix
}

cn_netScores<-function
### return the GRN establishment score for a given expression matrix
(expDat, 
 ### expression matrix
 genes, 
 ### vector of genes in the GRN
 tVals, 
 ### tVals list of ctt->list(means->genes, sds->genes)
 ctt, 
 ### ctt as basis of reference
 classList, 
 ### list of classifiers to extract weights
 classWeight=TRUE, 
 ### whether to weight by classifier or Not
 exprWeight=TRUE
 ### whether to weight by gene expression
){
  cat("genes:",length(genes),"\n");
  cat("classifier:", nrow(classList[[ctt]]$importance),"\n");
  aMat<-matrix(0, nrow=length(genes), ncol=ncol(expDat));
  rownames(aMat)<-genes;
  
  weights<-rep(1, length(genes));
  names(weights)<-genes;
  
  if(exprWeight){
    cat("expr weights\t");
    meanVect<-unlist(tVals[[ctt]][['mean']][genes]);
    weights<-(2**meanVect)/sum(2**meanVect);
  }
  
  
  if(classWeight){
    cat("class weights\n")
    classImp<-classList[[ctt]]$importance[genes,1];
    weights<-weights*classImp;
  }
  
  for(gene in genes){
    zzs<-cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]);
    aMat[gene,]<-zzs;
  }
  xscores<-apply(aMat, 2, weighted.mean, w=weights);
  xscores;
  # grn scores (normalized)
}

cn_rawScore<-function
### computes the raw score for a gene as xmax-abs(zscore). this makes better values higher.
(vect,
 ### a vector
 mmean,
 ### mean value in training data
 ssd
 ### standard deviation in training data
){
  zcs<-zscore(vect, mmean, ssd);
  xmax<-1000; # arbitrary, and corrected for later, but want some high enough that it should not be exceeded
  xmax-abs(zcs);
  # transformed (but not normalized) GRN score
}

cn_score<-function
### runs cn_netScores on each GRN
(expDat,
 ### expression matrix
 subList, 
 ### list of tct => genes
 tVals,
 ### a
 classList,
 ###
 minVals=NULL, 
 ### vector of subnet->minVal, # for shifting raw grn establishment values
 classWeight=TRUE,
 ### class weights?
 exprWeight=TRUE
 ### expression weights
){
  #nSubnets<-sum(sapply(subList, length));
  nSubnets<-length(subList);
  ans<-matrix(0, nrow=nSubnets, ncol=ncol(expDat));
  ctts<-names(subList);
  rnames<-vector();
  rIndex<-1;
  for(ctt in ctts){
    # cat(ctt,"\n");
    genes<-subList[[ctt]];
    #    snNames<-names(subnets);
    #    rnames<-append(rnames, snNames);    
    #    for(sName in snNames){
    ans[rIndex,]<-cn_netScores(expDat, genes, tVals=tVals, ctt=ctt,classList, classWeight=classWeight,exprWeight=exprWeight);
    rnames<-append(rnames, ctt);
    rIndex<-rIndex+1;
    #   }
  }
  rownames(ans)<-rnames;
  colnames(ans)<-colnames(expDat);
  if(!is.null(minVals)){
    minVals<-minVals[rownames(ans)];
    ans<-ans-minVals;
  }
  ans;
  ### GRN scores
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
