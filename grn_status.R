# CellNet
# (C) Patrick Cahan 2012-2014

# GRN status or establishment related functions. Setting up using training data, and calculating for query data.


cn_rawScore<-function
### computes the raw score for a gene as xmax-abs(zscore). this makes better values higher.
(vect, ### a vector of gene expression values for multiple samples
 mmean, ### mean value in training data
 ssd ### standard deviation in training data
){
  zcs<-zscore(vect, mmean, ssd);
  xmax<-1000; # arbitrary, and corrected for later, but want some high enough that it should not be exceeded
  xmax-abs(zcs);
  # transformed (but not normalized) GRN score
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
 classList=NULL, 
 ### list of classifiers to extract weights
 classWeight=FALSE, 
 ### whether to weight by classifier or Not
 exprWeight=TRUE
 ### whether to weight by gene expression
){

##  cat("genes:",length(genes),"\n");
##  cat("classifier:", nrow(classList[[ctt]]$importance),"\n");
  aMat<-matrix(0, nrow=length(genes), ncol=ncol(expDat));
  rownames(aMat)<-genes;
  
  weights<-rep(1, length(genes));
  names(weights)<-genes;
  
  if(exprWeight){
    meanVect<-unlist(tVals[[ctt]][['mean']][genes]);
    weights<-(2**meanVect)/sum(2**meanVect);
  }
    
  if(classWeight){
    classImp<-classList[[ctt]]$importance[genes,1];
    weights<-weights*classImp;
  }
  
  for(gene in genes){
    zzs<-cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]);
    aMat[gene,]<-zzs;
  }
  xscores<-apply(aMat, 2, weighted.mean, w=weights);
  xscores;
  # grn scores (not normalized)
}

cn_score<-function
### runs cn_netScores on each GRN
(expDat,
 ### expression matrix
 subList, 
 ### list of tct => genes
 tVals,
 ### a
 classList=NULL,
 ###
 minVals=NULL, 
 ### vector of subnet->minVal, # for shifting raw grn establishment values
 classWeight=FALSE,
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
    ### cat(subNet,"\n")
    ans[subNet,]<- queryScores[subNet,] / ctrlScores[[subNet]];
  }
  colnames(ans)<-colnames(queryScores);
  ans;
  ### normalized expression matrix
}

#
# functions to enable GRN status metric
#

cn_trainNorm<-function # Figure out normalization factors for GRNs, and norm training data
(expTrain,
 stTrain,
 subNets, # named list of genes, one list per CTT, tct=>gene vector
 classList = NULL, # list of classifiers
 dLevel = "description1",
 tVals=NULL, ### useful when debugging
 classWeight=FALSE, ### weight GRN status by importance of gene to classifier
 exprWeight=TRUE ### weight GRN status by expression level of gene?
){
  if(is.null(tVals)){
    tVals<-cn_make_tVals(expTrain, stTrain, dLevel);
  }
  
  ctts<-as.vector(unique(stTrain[,dLevel]));
  scoreList<-list();
  normList<-list(); # a list of ctt->subnet->mean value
  minVect<-vector(); # a list of ctt->subnet->min value, used to shift raw grn est scores
  
  cat("calculating GRN scores on training data ...\n");
  tmpScores<-cn_score(expTrain, subNets, tVals, classList, minVals=NULL, classWeight=classWeight, exprWeight=exprWeight); 
  
  minVect<-apply(tmpScores, 1, min);
  names(minVect)<-rownames(tmpScores);
  
  # shift the raw scores so that min=0;
  tmpScores<-tmpScores - minVect;
  cat("norm factors\n");
  for(ctt in ctts){
    # determine nomalization factors
    ##snets<-names(subNets[[ctt]]);
    snets<-ctt;
    scoreDF<-cn_extract_SN_DF(tmpScores, stTrain, dLevel, snets);
    scoreDF<-cn_reduceMatLarge(scoreDF, "score", "description", "subNet");
    xdf<-scoreDF[which(scoreDF$grp_name==ctt),];
    tmpSNS<-as.list(xdf$mean);
    names(tmpSNS)<-xdf$subNet;
    normList[names(tmpSNS)]<-tmpSNS;      
  }
  
  # normalize training scores
  nScores<-cn_normalizeScores(normList, tmpScores, rownames(tmpScores));
  scoreDF<-cn_extract_SN_DF(nScores, stTrain, dLevel);
  scoreDF<-cn_reduceMatLarge(scoreDF, "score", "description", "subNet");

  list(trainingScores=scoreDF,
       normVals=normList,
       raw_scores=tmpScores,
       minVals=minVect,
       tVals=tVals);
}



cn_make_tVals<-function### estimate gene expression dist in CTs
(expDat, ### training data 
 sampTab, ### training sample table 
 dLevel="description1", ### column to define CTs
 predictSD=FALSE ### whether to predict SD based on expression level
){
  
  if(predictSD){
    ans<-cn_make_tVals_predict(expDat, sampTab, dLevel);
  }
  else{
    # Note: returns a list of dName->gene->mean, sd, where 'dName' is a ctt or lineage 
    # make sure everything is lined up
    expDat<-expDat[,rownames(sampTab)];
    tVals<-list();
    dNames<-unique(as.vector(sampTab[,dLevel]));
    allGenes<-rownames(expDat);
    for(dName in dNames){
      #cat(dName,"\n");
      xx<-which(sampTab[,dLevel]==dName);
      sids<-rownames(sampTab[xx,]);
      xDat<-expDat[,sids];
      means<-apply(xDat, 1, mean);
      sds<-apply(xDat, 1, sd);
      tVals[[dName]][['mean']]<-as.list(means);
      tVals[[dName]][['sd']]<-as.list(sds);
    }
    ans<-tVals;
  }
  ans;
}


cn_make_tVals_predict<-function ### predicts SD based on mean expression
(expDat, ### training data 
 sampTab, ### training sample table 
 dLevel="description1" ### column to define CT
){
  # Note: returns a list of dName->gene->mean, sd, where 'dName' is a ctt or lineage 
  # make sure everything is lined up
  expDat<-expDat[,rownames(sampTab)];
  tVals<-list();
  dNames<-unique(as.vector(sampTab[,dLevel]));
  allGenes<-rownames(expDat);
  
  # make a model to predict SD given average expression level across all samples
  sdT<-apply(expDat, 1, sd);
  mT<-apply(expDat, 1, mean);
  myModel<-lm(sdT~mT);
  for(dName in dNames){
    xx<-which(sampTab[,dLevel]==dName);
    sids<-rownames(sampTab[xx,]);
    xDat<-expDat[,sids];
    means<-apply(xDat, 1, mean);
    sds<-predict(myModel, data.frame(mT=means));
    tVals[[dName]][['mean']]<-as.list(means);
    tVals[[dName]][['sd']]<-as.list(sds);
  }
  tVals;
}






