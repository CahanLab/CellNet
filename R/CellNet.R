# CellNet
# (C) Patrick Cahan 2012-2016


#
#' make a classifier, with a randomized class too
#'
#' rmake a classifier, with a randomized class too
#' @param expTrain training data
#' @param genes vector of genes to use as predictors
#' @param groups named vector of cells to groups or classes
#' @param nRand =50 num of randomized profiles to make
#' @param ntrees =2000 number of trees to build

#' @return RF
#' @export
#'
sc_makeClassifier<-function(
  expTrain,
  genes,
  groups,
  nRand=50,
  ntrees=2000){


  randDat<-randomize(expTrain, num=nRand)
  expTrain<-cbind(expTrain, randDat)

  allgenes<-rownames(expTrain)

  missingGenes<-setdiff(unique(genes), allgenes)
  cat("Number of mussing genes ", length(missingGenes),"\n")
  ggenes<-intersect(unique(genes), allgenes)
  randomForest(t(expTrain[ggenes,]), as.factor(c(groups, rep("rand", ncol(randDat)))), ntree=2000)

}

#
#' classify samples
#'
#' classify samples
#' @param rfObj result of running sc_makeClassifier
#' @param expQuery expQuery
#' @param numRand numRand

#' @return classRes matrix
#' @export
#'
rf_classPredict<-function(
  rfObj,
  expQuery,
  numRand=50){

    randDat<-randomize(expQuery, num=numRand)
    expQuery<-cbind(expQuery, randDat)

    preds<-rownames(rfObj$importance)
    xpreds<-t(predict(rfObj, t(expQuery[preds,]), type='prob'))
  colnames(xpreds)<-colnames(expQuery)
  xpreds
}


#' Apply CellNet to query data
#'
#' Classifies query data and computes GRN status
#'
#' @param expQuery matrix of expression values
#' @param stQuery data.frame expression meta data
#' @param cnProc CellNet object result of running cn_make_processor
#' @param develLevelQuery column name in stQuery to grou query samples for grn status plotting
#'
#' @examples
#' cnRes<-cn_apply(expQuery, stQuery, cnProc, "description1")
#'
#' @return cnRes object
#' @import igraph
#' @import ggplot2
#' @import randomForest
#' @import pheatmap
#' @import plyr
#' @import tidyr
#' @import parallel
#'
#' @export
cn_apply<-function
(expQuery,
 stQuery,
 cnProc,
 dLevelQuery="description1"
){
  
  # make sure stTrain has a sample_id col, and rownames are equivalent to it
  if(!any(colnames(stQuery)=='sample_id')){
    stQuery<-cbind(stQuery, sample_id=colnames(expQuery));
    rownames(stQuery)<-as.vector(stQuery[,"sample_id"])
  }

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
  scoresQuery<-cn_score(expQuery, subNets, tVals, minVals=minVals, classList=classList, classWeight=classWeight, exprWeight=exprWeight);

  # normalize query scores
  normScoresQuery<-cn_normalizeScores(normVals, scoresQuery, rownames(scoresQuery));

  # classify the query data
  myClass<-cn_classify(classList,expQuery,subNets);

  ###colnames(myClass)<-as.character(stQuery[,dLevelQuery])

  # it's really useful to have the averaged training data, too
  meanExpQuery<-as.matrix(GEP_makeMean(expQuery,as.character(stQuery[,dLevelQuery]),type='mean'));

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

#' Helper function to write a markdown line for the Readme
#'
#' This writes a string that includes the number of samples, the name and number of C/Ts, and S3 links for a given cnProc
#'
#'
#' @param cnProc cnProc
#' @param species species
#' @param date date made
#' @param url url 
#'
#' @return string
#'
#' @examples
#' myMDstring<-cn_writeReadme(cnProc, "mouse", date="Oct_25_2016")
#'
#' @export
cn_writeReadme<-function
(cnProc,
  species,
  date=NA,
  url="https://s3.amazonaws.com/CellNet/rna_seq/"
){
  url<-paste0(url, species,"/")
  if(is.na(date)){
    date<-utils_myDate()
  }
  ans<-paste0("| ", species, " | ", date, " | ")
  stTrain<-cnProc[['stTrain']]
  ctCounts<-table(stTrain[, cnProc[['dLevelTrain']]])
  for(i in 1:(length(ctCounts)-1)){
    ct<-names(ctCounts)[i]
    ccount<-ctCounts[[ct]]
    ans<-paste0(ans, ct , " (", ccount, "), ")
  }
  ct<-names(ctCounts)[i+1]
  ccount<-ctCounts[[ct]]
  ans<-paste0(ans, ct , " (", ccount, ") | ")
  # cnProc
  myLink<-paste0(url, "cnProc_RS_", species, "_",date,".rda")
  ans<-paste0(ans, "[cnProc](",myLink,") | ")

  # metadata
  myLink<-paste0(url, "sampTab_RS_", species, "_",date,".rda")
  ans<-paste0(ans, "[metadata](",myLink,") | ")

  # expression data
   myLink<-paste0(url, "expList_RS_", species, "_",date,".rda")
   ans<-paste0(ans, "[expression data](",myLink,") |")

  ans;
}


#' Write out classification scores, normalized data to csv
#'
#' and it compresses them, too
#'
#' @param cnRes result of running cn_apply
#' @param tfScores network inlfuence scores
#' @param prefix filename prefix
#'
#' @return filename
#'
#' @export
cn_outputRes<-function
(cnRes,
 tfScores,
 prefix
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
  tfStab<-cn_makeTFtable(tfScores);
  fname4<-paste(prefix, "_NIS.csv", sep='');
  write.csv(tfStab, file=fname4, row.names=F);
  
  fnameOut<-paste(prefix,"_data.tar", sep='')
  cmd<-paste("tar -c ", fname1, " ",fname2," ", fname3," > ", fnameOut, sep='');
  system(cmd);
  cmd<-paste("gzip -f ", fnameOut, sep='');
  system(cmd);
  fnameOut;
  ### 
}



#' make a CellNet object
#'
#' Does lots of stuff, but mostly you will care about the fact that it trains C/T classifiers
#' @param expTrain expression matrix
#' @param stTrain sample table
#' @param ctGRNs result of running cn_grnDoRock
#' @param dLevel sample table column name to operate on
#' @param classWeight weight GRN est by classification importance
#' @param exprWeight weight GRN est by gene expression
#' @param sidCol column in stTrain that contains unique id for each sample and is the colname in expTrain
#' @return cnProc
#'
#' @examples
#'
#'
#' @export
cn_make_processor<-function # train a CellNet object
(expTrain,
 stTrain,
 ctGRNs,  
 dLevel="description1",
 classWeight=TRUE,
 exprWeight=TRUE,
 sidCol="sample_id"
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
  ###classList<-cn_makeRFs(expTrain, stTrain, gListsSub, dLevel=dLevel);
  cat("Done making classifiers :)\n");
  
  trainNorm<-cn_trainNorm(expTrain, stTrain, subNets=geneLists, classList=classList, dLevel=dLevel, classWeight=classWeight, exprWeight=exprWeight, sidCol=sidCol);
  ###trainNorm<-cn_trainNorm(expTrain, stTrain, subNets=gListsSub, classList=classList, dLevel=dLevel, classWeight=classWeight, exprWeight=exprWeight);
   
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

#' make a CellNet object from an existing one with subset of genes and/or new expTrain
#'
#' Does lots of stuff, but mostly you will care about the fact that it trains C/T classifiers
#' @param cnProc old cnProc
#' @param newGenes optional gene subset to use 
#' @param expTrain optional new expTrain
#' @param ctGRNs result of running cn_grnDoRock
#' @param dLevel sample table column name to operate on
#' @param classWeight weight GRN est by classification importance
#' @param exprWeight weight GRN est by gene expression
#' @param sidCol column in stTrain that contains unique id for each sample and is the colname in expTrain
#' @return cnProc
#'
#' @examples
#'
#'
#' @export
cn_remake_processor<-function # train a CellNet object
(cnProc,
 newGenes=NA,
 expTrain=NA,
 sidCol="sample_id"
 ){
  
  if(is.na(newGenes) & is.na(expTrain)){
    ans<-cnProc;
  }
  else{
    # if a new expTrain is supplied
    if(is.na(expTrain)){
      expTrain<-cnProc$expTrain
    }
    # if newGenes is supplied, get right subset
    if(!is.na(newGenes)){

      oldGenes<-rownames(expTrain)
      newGenes<-intersect(newGenes, oldGenes)
      expTrain<-expTrain[newGenes,]
    }

    ctGRNs<-cnProc$ctGRNs 
    stTrain<-cnProc$stTrain
    dLevel<-cnProc$dLevelTrain
    classWeight<-cnProc$classWeight
    exprWeight<-cnProc$exprWeight

    geneLists<-ctGRNs[['ctGRNs']][['geneLists']]
  
    commonTypes<-intersect(stTrain[,dLevel], names(geneLists));

    ctGRNs[['ctGRNs']][['geneLists']]<-ctGRNs[['ctGRNs']][['geneLists']][commonTypes];
    ctGRNs[['ctGRNs']][['graphLists']]<-ctGRNs[['ctGRNs']][['graphLists']][commonTypes];
    ctGRNs[['ctGRNs']][['tfTargets']]<-ctGRNs[['ctGRNs']][['tfTargets']][commonTypes];
    geneLists<-ctGRNs[['ctGRNs']][['geneLists']];
  

    # Make classifiers
   cat("Making classifiers (this can take awhile) ...\n");
    classList<-cn_makeRFs(expTrain, stTrain, geneLists, dLevel=dLevel);
    ###classList<-cn_makeRFs(expTrain, stTrain, gListsSub, dLevel=dLevel);
    cat("Done making classifiers :)\n");
  
    # fix the spec genes so as to exclude genes in training but not in the q


    trainNorm<-cn_trainNorm(expTrain, stTrain, subNets=geneLists, classList=classList, dLevel=dLevel, classWeight=classWeight, exprWeight=exprWeight, sidCol=sidCol);
    ###trainNorm<-cn_trainNorm(expTrain, stTrain, subNets=gListsSub, classList=classList, dLevel=dLevel, classWeight=classWeight, exprWeight=exprWeight);
   
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
  }
  ans;
}




# convenience function to split make and assess classifiers

#' make classifiers and assess
#'
#' split the data by prop, make classifiers, assess them on held out data
#' @param stDat sample table
#' @param expDat expression data
#' @param grns result of running cn_grnDoRock
#' @param prop numeric 0<prop<1 what proportion of data to train classifiers on
#' @param dLevel sample table column to use as response value in classifier (e.g. cell type)
#' @param dLevelStudy column name to indicate experiment or study id
#' @param dLevelSID column to indicate sample 
#'
#' @return list if classifiers, classifier result as applied to held out data, ROCs
#'
#' @examples
#' assessmentResults<-cn_splitMakeAssess(stTrain, expTrain, ctGRNs, prop=.5, dLevel="description1")
#'
#' @export
cn_splitMakeAssess<-function
(stDat,
 expDat,
 grns,
 prop=0.5,
 dLevel="description1",
 dLevelStudy="exp_id",
 dLevelSID="sample_id"){

  cat("splitting data...\n");
  stList<-samp_for_class(stDat, prop=prop, dLevel=dLevel, dLevelStudy=dLevelStudy)
  stVal<-stList[['stVal']]; 
  stTrain<-stList[['stTrain']]; 


  expTrain<-expDat[,rownames(stTrain)] 
  expVal<-expDat[,rownames(stVal)] 

  cat("training size: ",ncol(expTrain),"\n");
  cat("making classifiers...\n");


  if(FALSE){  ### Was used when I was trying out a differnt method to ID CT GRNs
    gListsSub<-lapply(grns$ctGRNs$geneLists, sample, 100);
    for(ct in names(gListsSub)){
      gListsSub[[ct]]<- union(gListsSub[[ct]], names(ctGRNs[['ctGRNs']][['tfTargets']][[ct]]));
    }
  }

  ###testRFs<-cn_makeRFs(expTrain, stTrain, gListsSub, dLevel=dLevel)
  testRFs<-cn_makeRFs(expTrain, stTrain, grns$ctGRNs$geneLists, dLevel=dLevel)

  cat("assessing...\n")
  ansVal<-cn_classify(testRFs, expVal, grns$ctGRNs$geneLists)
  ###ansVal<-cn_classify(testRFs, expVal, gListsSub)
  assessed<-cn_classAssess(ansVal, stVal, classLevels=dLevel, dLevelSID=dLevelSID, resolution=0.01);
  list(classifiers=testRFs, classRes=ansVal, PRs=assessed);
}

#' classify data
#'
#' run binary classifiers on input expression matrix
#' @param classList result of running CN3_makeClassifier
#' @param expDat  expression data matrix
#' @param cttComms list of genes that were used to train each classifer cttComms[[className]]<-c(gene1, gene2, ...)
#'
#' @return classification matrix nrow=length(classList) ncol=ncol(expDat)
#'
cn_classify<-function
(classList, 
 expDat, 
 cttComms 
){
  ctts<-names(classList);
  ans<-matrix(0, nrow=length(ctts), ncol=ncol(expDat));
  rownames(ans)<-ctts;
  for(ctt in ctts){
    #cat(ctt,"\n")
    ##    xgenes<-cttComms[[ctt]][[1]];
    xgenes<-cttComms[[ctt]];
    ### 06-05-16
    xgenes<-intersect(xgenes, rownames(expDat));
    myProbs<-predict(classList[[ctt]], t(expDat[xgenes,]), type='prob');
    ans[ctt,]<-myProbs[,ctt];
  }
  colnames(ans)<-colnames(expDat);
  ans;
  # classification matrix
}


#' split data into train vs test
#'
#' Splits a sample table into 2 sample tables roughly by prop in which no samples with sampe exp_id are in both the test and train
#' @param sampTab sample table must have exp_id column
#' @param prop proportion of samples to use as training
#' @param dLevel column name for 
#' @param dLevelStudy column name to indicate experiment or study id
#'
#' @return list of stTrain stVal
samp_for_class<-function# return a sampTab for training a test classifier
(sampTab,
  prop=0.5,
  dLevel="description1",
  dLevelStudy="exp_id"){

  ctts<-unique(as.vector(sampTab[,dLevel]));
  stTrain<-data.frame();
  stVal<-data.frame();
  for(ctt in ctts){
    stTmp<-sampTab[sampTab[,dLevel]==ctt,];
    stTrain<-rbind(stTrain, subSamp_for_class(stTmp,prop, dLevelStudy));

    idsval<-setdiff( rownames(stTmp), rownames(stTrain) );
    stVal<-rbind(stVal, stTmp[idsval,]);
  }
  list(stTrain=stTrain, stVal=stVal);
}

#' select from sample table
#'
#' that's it. helper function
#' @param sampTab sample table
#' @param prop fraction of samples to select
#' @param dLevelStudy column name to indicate experiment or study id
#'
#' @return stTrain
subSamp_for_class<-function
(sampTab,
  prop=0.5,
  dLevelStudy="exp_id"){


  if(is.null(dLevelStudy)){
    ccount<-floor(nrow(sampTab)*prop)
    stTrain<-sampTab[sample(rownames(sampTab), ccount),]
  }
  else{
    expIDcounts<-sort(table(sampTab[,dLevelStudy]));
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
      stTrain<-rbind(stTrain, sampTab[sampTab[,dLevelStudy]==expID,]);
    }
  }
  stTrain;
}

#' make classifiers
#'
#' Make one Random Forest classification per unique dLevel
#' @param expTrain expression matrix
#' @param stTrain sample table
#' @param list of dLevel=>vector of genes
#' @param dLevel stTrain column name to split on
#'
#' @return list of classifiers
cn_makeRFs<-function#
(expTrain,
 stTrain,
 geneLists, 
 dLevel='description1'){
  ans<-list();
  cnames<-names(geneLists);
  cnames<-intersect(cnames, unique(as.vector(stTrain[,dLevel])));
  for(cname in cnames){
    cat(cname,"\n");
    xgenes<-intersect(rownames(expTrain), geneLists[[cname]]);  
    cat(cname ,":", length(xgenes),"\n");
    ans[[cname]]<-cn_makeClassifier(expTrain[xgenes,],cname,stTrain, dLevel=dLevel);
  }
  ans;
}

#' make a single RF classifier
#'
#' uses R's randomForest package
#' @param trainScores numeric matrix of predictors
#' @param ctt what cell type/tissues is this?
#' @param stTrain sample table
#'
#' @return randomForest classifier
cn_makeClassifier<-function 
(trainScores, 
 ctt,
 stTrain,
 dLevel='description1'){
  resp<-as.vector(stTrain[,dLevel]);
  xi<-which(resp!=ctt);
  resp[xi]<-'other';  
  myClass<-randomForest(t(trainScores), as.factor(resp), ntree=2e3);
  myClass;
}







