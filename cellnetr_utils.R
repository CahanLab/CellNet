cn_geoFetch<-function# Fetch GEO RAW data from GEO and decompress
(sampTab,
 ### data.frame of sample_id, sample_name, description[1-x], exp_id, file_name
 maxTries=10,
 ### number of times to try to fetch a GEO archive
 decompress=FALSE
){
  accessions<-as.vector(sampTab$sample_id);
  keep<-vector();
  
  # keep track of downloaded files
  fnames<-vector();
  xi<-1;
  tryCount<-1;
  while(xi<=length(accessions)){
    fname<-accessions[xi];
    # sometimes the lines are blank
    if(nchar(fname)>0){
      cat(fname,"\n");
      myTry<-try(getGEOSuppFiles(fname, makeDirectory=FALSE));
      
      if(class(myTry) == "try-error" && tryCount<=maxTries){
        cat("tryCount:",tryCount,"\n");
        tryCount<-(tryCount+1);
        Sys.sleep(5);
        cat(fname,"\n");
        next();
      }
      fnames<-append(fnames, paste(fname, "_RAW.tar", sep=''));
      if(decompress){
        cmd<-paste("gzip -df *.gz");
        system(cmd, ignore.stdout=TRUE, wait=TRUE);
      }
    }
    xi<-xi+1;
    tryCount<-1;
  }
  fnames;
  ### data.frame, same columns as sampTab, but with corrected file_names    
}


utils_loadObject<-function
### loads an R object when you don't know the name
(fname
 ### file name
){
  x<-load(fname);
  get(x);
}

utils_stripwhite<-function
### strip whitespace from a string
(string
 #### string
 ){
  gsub("^\\s+|\\s+$", "", string)
}

utils_myDate<-function
### print date
()
{
  format(Sys.time(), "%b_%d_%Y");
}

GEP_makeMean<-function
### return a dataframe of mean or median-ed data based on given groupings
(exp,
 ### exp df
 groupings,
 ### vector of groupings
 type='mean'
 ### mean or median
){
  
  ##<<note colnames become the column name of the first sample in each group from the original data
  ans<-data.frame();
  grps<-unique(groupings);
  if(type=='mean'){
    for(grp in grps){
      gi<-which(groupings==grp);
      if(length(gi)==1){
        cat("\nhere")
        if(nrow(ans)==0){
          ans<-data.frame(exp[,gi]);
        }else{
          ans<-cbind(ans, exp[,gi]);
        }
      }
      else{
        xxx<-apply(exp[,gi],1,mean);
        if(nrow(ans)==0){
          ans<-data.frame(xxx);
        }
        else{
          ans<-cbind(ans, xxx);
        }
      }
    }
  }
  else{
    for(grp in grps){
      gi<-which(groupings==grp);
      xxx<-apply(exp[,gi],1,median);
      if(nrow(ans)==0){
        ans<-data.frame(xxx);
      }
      else{
        ans<-cbind(ans, xxx);
      }
    }
  }
  
  colnames(ans)<-grps;
  ans;
  ### data.frame of mean or median-ed data based on given groupings
}


zscore<-function
### compute zscore
(x,
 ### numeric vector
 meanVal, 
 ### mean of distribution to compute zscore of x against 
 sdVal
 ### standard deviation of distribution to compute zscore of x agains
 ){ 
  (x-meanVal)/sdVal;
  ### zscore
}

cn_zscoreVect<-function
### Compute the mean zscore of given genes in each sample
(genes,
 ### genes
 xvals,
 ### named vector
 tVals,
 ### tvals
 ctt
 ### ctt
 ){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
  ### zscore vector
}

cn_reduceMatLarge<-function
### reduce a data.frame
(datFrame,
 ### df
 valCol="score",
 ### value to merge
 cName="description",
 ### column to reduce on 
 iterOver="subNet"
 ### iterate over
 ){
  
  iterOvers<-unique(as.vector(datFrame[,iterOver]));
  ans<-data.frame();
  for(io in iterOvers){
    #  cat(io,"\n");
    xi<-which(datFrame[,iterOver]==io);
    dfTmp<-datFrame[xi,];
    x<- utils_reduceMat(dfTmp,valCol=valCol,cName=cName);
    x<-cbind(x, subNet=rep(io, nrow(x)));
    ans<-rbind(ans, x);
  }
  ans;
  ### ans
}

cn_extract_SN_DF<-function
### returns a DF of: sample_id, description, ctt, subnet_name, score
(scores,
 ### a matrix of subNet scores
 sampTab,
 ### sample table
 dLevel,
 ### column name of sampTab to group values by
 rnames=NULL
 ### rownames to extract
){
  
  if(is.null(rnames)){
    rnames<-rownames(scores);
    #cat("GOT NULL\n");
  }
  tss<-scores[rnames,];
  if(length(rnames)==1){
    tss<-t(as.matrix(scores[rnames,]));
    rownames(tss)<-rnames;
  #  cat(dim(tss),"\n")
  }
  nSamples<-ncol(tss);
  stTmp<-sampTab[colnames(tss),]; ####
  snNames<-rownames(tss);
  num_subnets<-length(snNames);    
  snNames<-unlist(lapply(snNames, rep, times=nSamples));
  sample_ids<-rep(as.vector(stTmp$sample_id), num_subnets);
  descriptions<-rep(as.vector(stTmp[,dLevel]), num_subnets);
    # myCtts<-rep(ctt, length(snNames));
  scores<-as.vector(t(tss));
  data.frame(sample_id=sample_ids,
             description=descriptions,
             #         ctt=myCtts,
             subNet = snNames, 
             score=scores);
  ### data.frame
}

utils_reduceMat<-function
### reduce the data.matrix values by averaging and getting st dvs
(datFrame,
 ### data frame to reduce
 valCol,
 ### value column
 cName='ann'
 ### cname to reduce on e.g. ann
){
    
  mids<-unique(as.vector(datFrame[,cName]));
  means<-rep(0, length(mids));
  sds<-rep(1, length(mids));
  indexI<-rep(0, length(mids)); # index to extract other columns from original data.frame
  
  for(i in seq(length(mids))){
    mid<-mids[i]
    xi<-which(datFrame[,cName]==mid);
    tmpDat<-datFrame[xi,];
    means[i]<-mean(tmpDat[,valCol]);
    sds[i]<-sd(tmpDat[,valCol]);
  }  
  data.frame(grp_name=mids, mean=means, stdev=sds);
  ### df of grp_name, mean, sd
}

utils_stderr<-function
### calculate standard error
(x){
  sqrt(var(x)/length(x));
  ### stderr
}

utils_concord<-function 
### concordance is a useful measure of set overalp 
(vect1,
 ### vector 1
 vect2
 ### vector 2
 ){
  length(intersect(vect1, vect2))/length(union(vect1, vect2));
  ### length(intersection) / length(union)
}

cellnetr_log<-function
#### make a CellNet logging string
(str
 ### str to log
 ){
  lstring1<-paste("CellNet [",format(Sys.time(), "%Y_%b_%d_%H:%M:%S"),"]: ", sep='');
  paste(lstring1, str, "\n", sep='');
  ### output str
}

cn_remergeGRNs<-function### reduce redundancy of the lineage and general GRNS
(ctGRNsgen,### from grnStuff [['ctGRNs']]
 ctGRNgl){
   
  ctts<-names(ctGRNgl[['graphs']]);
  newGRs<-list();
  newGLs<-list();
  for(ct in ctts){
    cat(ct, "\t");
    a1<-ctGRNsGen[['graphs']][[ct]];
    cat(length(a1),"\t");
    a2<-ctGRNgl[['graphs']][[ct]];
    cat(length(a2), "\n");
    a1v2<-cn_graphMerge(a1);
    a2v2<-cn_graphMerge(a2);
    mergedGRN<-cn_graphMerge(list(l1=a1v2, l2=a2v2));
    newGRs[[ct]]<-mergedGRN;
    newGLs[[ct]]<-unique(V(mergedGRN)$name);
  }
  list(graphs=newGRs,geneLists=newGLs);
}

cn_graphMerge<-function### merge graphs
(
  igsList ### list og iGraph graphs
){
  nnames<-names(igsList);
  newGraph<-igsList[[nnames[1]]];
  for(i in 2:length(nnames)){
    newGraph<-graph.union(newGraph, igsList[[nnames[i]]], byname=TRUE);
    # select the maximum zscore for weighting an edge
    wNames<-list.edge.attributes(newGraph);
    x1<-get.edge.attribute(newGraph, wNames[1]);
    weights<-matrix(0, ncol=length(wNames), nrow=length(x1));
    for(j in seq(length(wNames))){
      weights[,j]<- get.edge.attribute(newGraph, wNames[j]);
    }
    gWeights<-apply(weights, 1, max, na.rm=T);
    E(newGraph)$weight<-gWeights;
  }
  newGraph;
}




