#
# Copyright (C) 2014 Patrick Cahan
#
cn_outputRes<-function
### write out grn scores, classification scores, normalized data to csv, then tar and  zip
(cnRes,
 ### cnRes object
 tfScores,
 ## tfscores
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
  
  # write out normalized query expression data
  expDat<-cnRes[['expQuery']];
  fname3<-paste(prefix, "_expNorm.csv", sep='');
  write.csv(expDat, file=fname3);
  
  
  # write out NIS
  nisTab<-.cn_makeTFtable(tfScores);
  fname4<-paste(prefix, "_NIS.csv", sep='');
  write.csv(nisTab, file=fname4, row.names=F);
  
  c(fname1, fname2, fname3, fname4);
  
 ### file names of pdfs 
}


CN3_apply_SN<-function### get subnet establishment levels; need to figure out ctt based on snName
(expQuery,
 ### expression matrix
 stQuery,
 ### sample table
 cnProc, 
 ### cnProc object
 sn_norm, 
 ### result of running CN3_trainNormSN
 dLevelQuery="description1"
 ### stQuery column name
){
  
  ctrlScores<-sn_norm[['trainingScores']];   # subnet est scores in control  -- normalized, df
  normVals<-sn_norm[['normVals']];              # average subnet est scores in control samples,  (list of ctt->subnet ave)
  minVals<-sn_norm[['minVals']];             # min raw vals of grn establishments to shift by
  rawCtrlScores<-sn_norm[['raw_scores']];
  tVals<-cnProc[['tVals']];
  subNets<-sn_norm[['grnList']];
  glLists<-.CN3_extractGL(subNets);
  classList<-cnProc[['classList']];
  classWeight<-sn_norm[['classWeight']];
  
  # score the query data
  cat("Scoring query data...\n")
  scoresQuery<-CN3_NS_subnet(expQuery, glLists, tVals, minVals=minVals, classList, classWeight=classWeight);
  
  # normalize query scores
  cat("normalizing grn scores\n");
  normScoresQuery<-CN3_normalizeScores(normVals, scoresQuery, rownames(scoresQuery));
  
  ans<-list(queryScores=scoresQuery,
            normScoresQuery=normScoresQuery,
            stQuery=stQuery,
            dLevelQuery=dLevelQuery);            
  ans;  
  
  ### Subnet establishment values
}

CN3_NS_subnet<-function#runs CN3_netScores on each GRN
(expDat,
 ### expression matrix
 subList, 
 ### list of sub_network->genes
 tVals,
 ### tvals
 classList,
 ### class list
 minVals=NULL, 
 ### vector of subnet->minVal, # for shifting raw grn establishment values
 classWeight=TRUE
 ### classWeight? default=TRUE
){
  #nSubnets<-sum(sapply(subList, length));
  nSubnets<-length(subList);
  ans<-matrix(0, nrow=nSubnets, ncol=ncol(expDat));
  snNames<-names(subList);
  rnames<-vector();
  rIndex<-1;
  for(snName in snNames){
    ctt<-.get_cttName(snName);
#    cat(ctt,"\n");
    genes<-subList[[snName]];
    cat(snName, " ",length(genes), "\n");
    #    snNames<-names(subnets);
    #    rnames<-append(rnames, snNames);    
    #    for(sName in snNames){
    ans[rIndex,]<-CN3_netScores(expDat, genes, tVals=tVals, ctt=ctt,classList, classWeight=classWeight);
    rnames<-append(rnames, snName);
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
  ### subnet-GRN values
}

.get_cttName<-function
### splits a sub-network name into a ctt
(snName
 ### subnet name
 ){
  x<-strsplit(snName, "_")[[1]][1];
  x;
}

cn_nis_all<-function
### NIS for all GRNs
(cnRes,
 ### cnRes object
 cnProc,
 ### CellNEt object used to produce cnRes
 ctt,
 ### target cell tissue type
 relaWeight=1
 ### whether to weight by overall expression such that TFs with higher expression in ctt are more important (1=do the weighting)
 ){
  snNames<-names(cnProc$ctGRNs$ctGRNs$graphLists);
  ans<-list()
  for(snName in snNames){
    cat("scoring ", snName,"\n")
    x<-cn_nis(cnRes, cnProc, snName, ctt,relaWeight);
    ans[[snName]]<-x;
#    x<-cbind(x, grn=rep(ctGRN, nrow(x)));
#    ans<-rbind(ans, x);
  }
  ans;
  # list of numeric matrix of TF scores
}


cn_nis<-function
### compute network influence score. see manuscript for details
(cnRes,
 ### cnRes object
 cnProc,
 ### CellNet object used to produce cnRes
 subnet,
 ### what subnet to evaluage
 ctt,
 ### what is the reference cell/tissue type
 relaWeight=1
 ### whether to weight by overall expression such that TFs with higher expression in ctt are more important (1=do the weighting)
){
  
  tfTargList<-cnProc[['ctGRNs']][['ctGRNs']][['tfTargets']];
  # return a DF of : tfs, nTargets, targetScore, tfScore, totalScore
  nTargets<-vector();
  targetScore<-vector();
  tfScore<-vector();
  totalScore<-vector();
  tfWeights<-vector();
  
  tfs<-names(tfTargList[[subnet]]);
  netGenes<-cnProc[['grnList']][[subnet]];
  netGenes<-intersect(netGenes, rownames(cnProc[['expTrain']]))
  
  expDat<-cnRes[['expQuery']];
  stQuery<-cnRes[['stQuery']];
  sids<-as.vector(stQuery$sample_id);
  
  ans<-matrix(0, nrow=length(tfs), ncol=nrow(stQuery));
  rownames(ans)<-tfs;
  colnames(ans)<-sids;
  
  tVals<-cnProc[['tVals']];
  
  # compute a matrix of zscores.
  zzzMat<-matrix(0, nrow=length(netGenes), ncol=nrow(stQuery));
  
  for(i in seq(length(sids))){
    sid<-sids[i];
    #cat("computing zscores ", sid,"\n");
    xvals<-as.vector(expDat[netGenes,sid]);
    names(xvals)<-netGenes;
    zzzMat[,i]<-cn_zscoreVect(netGenes, xvals, tVals, ctt);    
  }
  
  rownames(zzzMat)<-netGenes;
  colnames(zzzMat)<-rownames(stQuery);
  
  for(sid in sids){
    #cat("tf scoring ", sid,"\n");
    xvals<-as.vector(expDat[,sid]);
    names(xvals)<-rownames(expDat);
    
  
    # assign weights
    
    ### # meanVect<-unlist(tVals[[ctt]][['mean']][netGenes]);
    meanVect<-unlist(tVals[[subnet]][['mean']][netGenes]);
    weights<-(2**meanVect)/sum(2**meanVect);
    
    for(i in seq(length(tfs))){
      
      tf<-tfs[i];
      
      # zscore of TF relative to target C/T
##      tfScore[i]<-zscore(xvals[tf], tVals[[ctt]][['mean']][[tf]], tVals[[ctt]][['sd']][[tf]]);
      
      tfScore[i]<-zzzMat[tf,sid];
      
      targs<-tfTargList[[subnet]][[tf]];
      targs<-intersect(targs, rownames(cnProc[['expTrain']]));
      
      # Zscores of TF targets, relative to C/T
##      tmp<-cn_zscoreVect(targs, xvals, tVals, ctt );
      tmp<-zzzMat[targs,sid];
      targetScore[i]<-sum(tmp*weights[targs]);
      
      ## new one:
      totalScore[i]<-targetScore[i] + (length(targs)*tfScore[i]*weights[tf]);
      
      if(relaWeight!=1){ # don't weight by expression
        meanW<-mean(weights)
        totalScore[i]<- sum(tmp)*meanW + (length(targs)*tfScore[i])*meanW
      }
      nTargets[i]<-length(targs) ;
      tfWeights[i]<-weights[tf];
    }
    xxx<-data.frame(tf=tfs, tfScore=tfScore, targetScore=targetScore, nTargets=nTargets,tfWeight=tfWeights, totalScore=totalScore);
    xxx<-xxx[order(xxx$totalScore),]; # puts the worst ones at top when plotting
    xxx$tf<-factor(xxx$tf, as.vector(unique(xxx$tf)));
    ans[as.vector(xxx$tf),sid]<-as.vector(xxx$totalScore);
  }
  ans;
  # returns network influence score.
}

cn_subnet_dysregUp<-function
### find the subnets that are 'dysregulated' in the selected samples, higher in query than in select train ctts
(cnObj,
 ### cnRes obj
 cnProc, 
 ### CellNet object
 ctts, 
 ### cell type to which to compare to query samples
 grpName,
 ### which samples to test for dysregulation
 zThresh=2 
 ### z-score threshold for calling a GRN dysregulated.
){
  
  qScores<-cnObj[['queryScores']];
  sn_names<-rownames(qScores);
  
  ctrlScores<-cnProc[['raw_scores']];
  xx<-cn_extract_SN_DF(ctrlScores, cnProc[['stTrain']], cnProc[['dLevelTrain']],rnames=sn_names);
  xx3<-cn_reduceMatLarge(xx, "score", "description", "subNet");
  ctrlScores<-xx3;
  
  # convert into a data.frame
  aa<-cn_extract_SN_DF(qScores, cnObj[['stQuery']], cnObj[['dLevelQuery']], rnames=sn_names);
  aa3<-cn_reduceMatLarge(aa, "score", "description", "subNet");
  aa3<-cbind(aa3, src=rep('query', nrow(aa3)));
  
  ans<-list();
  for(i in seq(length(ctts))){
    ctt<-ctts[i];
    if(i==1){
      ans<-CN3_dr1(sn_names, ctrlScores, ctt, aa3, grpName,zThresh);
    }
    else{
      ans<-intersect(ans, CN3_dr1(sn_names, ctrlScores, ctt, aa3, grpName,zThresh));
    }
  }
  ans;
  # not implemented
}
