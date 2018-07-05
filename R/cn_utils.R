# CellNet
# (C) Patrick Cahan 2012-2016

# commonly used or misc functions


#
#' randomize data matrix 
#'
#' randomize data matrix 
#' @param expDat expDat
#' @param num number of profiles to return
#'
#' @return exp matrix random
#' @export
#'
randomize<-function(
 expDat,
 num=50){


  randDat<-t(apply(expDat, 1, sample))  
  randDat<-apply(randDat, 2, sample)

  randDat<-randDat[,sample(1:ncol(randDat), num)]
  colnames(randDat)<-paste0(rep("rand_", num), 1:num)
  rownames(randDat)<-rownames(expDat)
  randDat
}



#' 1-PCC distance
#'
#' 1-PCC distance
#' @param x numeric matrix
#' 
#' @return distance matrix  
#'
#' @examples
#' xdist<-utils_myDist(t(expDat))
#' plot(hclust(xdist, 'ave'), hang=-1)
#'
#' @export
utils_myDist<-function
(x
){
  as.dist(1-cor(t(x)));
}

#' loads an R object when you don't know the name
#'
#' loads an R object when you don't know the name
#' @param fname file
#'
#' @return variable
#'
#' @export
utils_loadObject<-function
(fname
 ### file name
){
  x<-load(fname);
  get(x);
}

#' strip whitespace from a string
#'
#' strip whitespace from a string
#' @param string string
#'
#' @return new string
#'
#' @export
utils_stripwhite<-function
### 
(string
 #### string
 ){
  gsub("^\\s+|\\s+$", "", string)
}

#' print date
#'
#' print date
#' @return string
#'
#' @export
utils_myDate<-function
### 
()
{
  format(Sys.time(), "%b_%d_%Y");
}

#' count comments in a file
#'
#' count comments in a file, assumes #
#' @param fname file
#'
#' @return numeric
#'
#' @export
utils_count_comments<-function
(fname){
  cmd<-paste("grep \"#\" ", fname,sep='');
  x<-system(cmd, intern=T);
  length(x);
}

#' tars and compresses (gzip) dir
#'
#' tars and compresses (gzip) dir
#' @param dir directory
#'
#' @return file name
#'
#' @export
utils_tarcompress_dir<-function 
(dir){
  fname<-paste(dir,".tar.gz", sep='');
  cmd<-paste("tar -cvzf ", fname," ",dir, sep='');
  system(cmd);
  fname;
}

#' reduces full path to filename
#'
#' reduces full path to filename
#' @param string
#'
#' @return something
#'
#' @export
utils_strip_fname<-function #
(str){
  a<-strsplit(str, "/")[[1]];
  a[length(a)];
}

#' make a logging string
#'
#' make a logging string
#' @param str to log
#' @param str to prepend if wanted
#'
#' @return output str
#'
#' @export
utils_log<-function
(str, 
  prepend="" 
 ){
  lstring1<-paste(prepend, "[",format(Sys.time(), "%Y_%b_%d_%H:%M:%S"),"]\t", sep='');
  paste(lstring1, str, "\n", sep='');
}

#' gzip and untar
#'
#' gzip and untar
#' @param fname file
#'
#' @return nothing
#' @export
utils_unpack<-function#
(fname){
  cmd<-paste("gzip -d ", fname, sep='');
  system(cmd);
  fname<-strsplit(fname, ".gz")[[1]][1];
  cmd<-paste("tar -xvf ", fname, sep='');
  system(cmd);
}

#' row average (or median) based on groups
#'
#' row average (or median) based on groups
#' @param exp expression df
#' @param groupings groupings
#' @param type mean or media
#'
#' @return return a dataframe of mean or median-ed data based on given groupings.  colnames become the column name of the first sample in each group from the original data
#'
#' @export
GEP_makeMean<-function
(exp,
 groupings,
 type='mean'
){
  
  
  ans<-data.frame();
  grps<-unique(groupings);
  if(type=='mean'){
    for(grp in grps){
      gi<-which(groupings==grp);
      if(length(gi)==1){
        
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


#' returns a DF of: sample_id, description, ctt, subnet_name, score
#'
#' returns a DF of: sample_id, description, ctt, subnet_name, score
#' @param scores a matrix of subNet scores
#' @param sampTab sample table
#' @param dLevel column name of sampTab to group values by
#' @param rnames rownames to extract
#' @param sidCol sample identifier column name
#'
#' @return returns a DF of: sample_id, description, ctt, subnet_name, score
cn_extract_SN_DF<-function
(scores,
 sampTab,
 dLevel,
 rnames=NULL,
 sidCol="sample_id"
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
  sample_ids<-rep(as.vector(stTmp[,sidCol]), num_subnets);
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



list_intersect<-function### applies 'intersect' to members of a list
(aList){
  ii<-2;
  a<-aList[[1]];
  while(ii <= length(aList)){
    a<-intersect(a, aList[[ii]]);
    ii<-ii+1;
  }
  a;
}

# s3 utils
s3_get_par<-function# get files in parallel
(bucket,
  path,
  fnames){
  tfname<-'tmpfile.txt';
  write.table(fnames,tfname, col.names=FALSE, row.names=FALSE, quote=FALSE);
  cmd<-paste("cat ", tfname, " | parallel aws s3 cp s3://", bucket, "/", path, "/{} ./", sep='');
  system(cmd);
}

#' get file from S3 using CLI
#'
#' get file from S3 using CLI, assumes AWS credentials are set
#' @param dir path
#' @param target filename
#' @param bucket bucket
#'
#' @return nothing
#'
#' @export
s3_get<-function
(dir,
 target, # name of target file
 bucket="pcahanrnaseq"
){
  fpath<-paste(bucket, "/",dir,sep='');
  cmd<-paste("aws s3 cp s3://", bucket, "/", dir, "/", target, " ./", sep='');
  system(cmd);  
}

#' put file from S3 using CLI
#'
#' put file from S3 using CLI, assumes AWS credentials are set
#' @param dir path
#' @param target filename
#' @param bucket bucket
#'
#' @return nothing
#'
#' @export
s3_put<-function
(dir,
 target, # name of target file
 bucket="pcahanrnaseq"
){
  fpath<-paste(bucket, "/",dir,sep='');
  cmd<-paste("aws s3 cp ",target," s3://", bucket, "/", dir, "/", sep='');
  cat(cmd,"\n");
  system(cmd);  
}


#' List directories in a bucket
#'
#' List directories in a bucket
#' @param bucket bucket
#' @param path path
#'
#' @return vector of directories
#'
#'
#' @export
s3_listDir<-function
(bucket,
 path
){
#  fpath<-paste(bucket, "/",path,sep='');
  cmd<-paste("aws s3 ls s3://", bucket, "/", path, "/", sep='');
  cat(cmd,"\n");
  aa<-system(cmd, intern=T)
  bb<-unlist(lapply(aa, utils_stripwhite))
  cc<-unlist(strsplit(bb, " "))
  dd<-cc[grep("/", cc)]
  unlist(strsplit(dd, "/"))
}


s3_file_append<-function
#### add a line to a file on S3, always start with date
(dir, ### directory on bucket
 target, ### filename 
 str, ### str to log
 bucket="pcahanrnaseq"
 ){
  # fetch file from S3
  s3_get(dir, target, bucket);
  # add a line
  cat(utils_log(str), file=target, append=TRUE);
  # send file back
  s3_put(dir, target, bucket);
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

#' make Inf and -Inf values sensible
#'
#' make Inf and -Inf values sensible
#' @param zMat zMat
#'
#' @return corrected zMat
#'
cn_correctZmat<-function
(zmat){
  myfuncInf<-function(vect){
    xi<-which(vect=='Inf')
    if(any(xi)){
      mymax<-max(vect[-xi])
      vect[xi]<-mymax
    }
    vect
  }
  zmat<-apply(zmat,2, myfuncInf)
  zmat[is.na(zmat)]<-0
  zmat
}

########################################
##
## UNUSED
##
########################################


