# CellNet
# (C) Patrick  2012-2016


#' set up directories on ephemeral drives that will be needed to store large files
#'
#' Creates a dir /media/ephemeral0/analysis, changes ownership ti ec2-user, sets cwd to this dir, sets up dir in /media/ephemeral1 for reference data
#' 
#' @param wrDir working directory
#'
#' @return nothing
#' @export
cn_setup<-function
(local=FALSE,
  wrDir="/media/ephemeral0/analysis")
{
  if(local) {
    system("mkdir CellNetLocal")
    ## Working Directory
    setwd("CellNetLocal")
    ## Directory for Indices and gene_trans file
    system("mkdir ref")
    ## tmp directory is for parallel
    system("mkdir ./tmp")
    Sys.setenv(TMPDIR="./tmp")
  }
  else {

   cmd<-paste0("sudo chown ec2-user /media/ephemeral0")
  system(cmd)
  cmd<-paste0("sudo mkdir ", wrDir)
  system(cmd)
  cmd<-paste0("sudo chown ec2-user ",wrDir)
  system(cmd)
  setwd(wrDir)

   cmd<-paste0("mkdir /media/ephemeral0/tmp")
  system(cmd) 
  Sys.setenv(TMPDIR="/media/ephemeral0/tmp")

  cmd<-paste0("mkdir /media/ephemeral0/analysis/tmp")
  system(cmd) 

 
    cmd<-paste0("sudo mount /dev/xvdc /media/ephemeral1")
    system(cmd)

    cmd<-paste0("sudo mkdir /media/ephemeral1/dat")
   system(cmd)
   cmd<-paste0("sudo mkdir /media/ephemeral1/dat/seq")
   system(cmd)
   cmd<-paste0("sudo mkdir /media/ephemeral1/dat/ref")
   system(cmd)

    cmd<-paste0("sudo chown ec2-user /media/ephemeral1/dat")
    system(cmd)

    cmd<-paste0("sudo chown ec2-user /media/ephemeral1/dat/seq")
   system(cmd)

    cmd<-paste0("sudo chown ec2-user /media/ephemeral1/dat/ref")
    system(cmd)
  }  
}

#' wrapper to fetch files needed to run salmon in case of errors
#'
#' fetch files needed to run salmon
#' @param destination where to put the indices
#' @param species mouse or human
#' @param bucket where to get them
#' @param dirw hat path on s3 to get them
#' @param salmonVersion version of salmon, if > 0.6 gets latest indices
#'
#' @return nothing
#' @export
fetchIndexHandler <- function
(destination = "/media/ephemeral1/dat/ref",
 species = "mouse",
 bucket='cellnet-rnaseq',
 dir="ref",
 iFile="salmon.index.mouse.050316.tgz") {
  tryCatch({
    fetch_salmon_indices(destination,species,bucket,dir, iFile=iFile)
  }, error = function(e) {
    print("Doesn't Always work the first time. Trying again!")
    setwd("..")
    fetch_salmon_indices(destination,species,bucket,dir, iFile=iFile)
  }
  )
}




#' fetch files needed to run salmon
#'
#' fetch files needed to run salmon
#' @param destination where to out the indices
#' @param species mouse or human
#' @param bucket where to get them
#' @param dir what path on s3 to get them
#'
#' @return nothing
#' @export
fetch_salmon_indices = function (destination = "/media/ephemeral1/dat/ref", species = "mouse", 
          bucket = "cellnet-rnaseq", dir = "ref", iFile = NA) 
{
  curdir <- system("pwd", intern = T)
  setwd(destination)
  if (species == "mouse") {
    fname2 = "salmon.index.mouse.050316.tgz"
    fname3 <- "geneToTrans_Mus_musculus.GRCm38.80.exo_Jun_02_2015.R"
  }
  else {
    fname2 = "salmon.index.human.050316.tgz"
    fname3 <- "geneToTrans_Homo_sapiens.GRCh38.80.exo_Jul_04_2015.R"
  }
  if(!is.na(iFile)) {
    fname2 = iFile
  }
  cat("fetching and unpacking stuff needed for Salmon ...\n")
  pref <- paste0("https://s3.amazonaws.com/", bucket, "/", 
                 dir, "/")
  download.file(paste0(pref, fname2), destfile = fname2)
  cmd <- paste("tar zxvf ", fname2, sep = "")
  system(cmd)
  download.file(paste0(pref, fname3), destfile = fname3)
  #download.file(paste0(pref, fname4), destfile = fname4)
  #cmd <- paste("gzip -d ", fname4, sep = "")
  system("rm *.tgz")
  setwd(curdir)
}





#' Derive gene expression estimates compatible with CellNet
#'
#' Derive gene expression estimates compatible with CellNet
#' @param sampTab sample table with fname column point to fastq files
#' @param total number of reads to normalzie to
#' @param fnameCol column name of fastq files
#' @param finalLength length of reads after trimming
#' @param species mouse or human
#' @param bucket where to get them
#' @param what path on s3 to get them
#'
#' @return nothing
#' @export
cn_salmon<-function
(sampTab,
 total=1e5,
 fnameCol='fname',
 finalLength=40,
 delOrig=FALSE,
 salmonIndex="MM_GRCh38.SalmonIndex.030816",
 geneTabfname="geneToTrans_Mus_musculus.GRCm38.80.exo_Jun_02_2015.R",
 refDir="/media/ephemeral1/dat/ref/",
 salmonPath="~/rnaseq/SalmonBeta-0.6.1_DebianSqueeze/bin"
  ){

  # make sure we are in the right place
 # setwd("/media/ephemeral0/analysis/")

  cat("determining read length\n")
  sampTab<-cbind(sampTab, readLength=unlist(lapply(as.vector(sampTab[,fnameCol]), fastq_readLength)))

  cat("Trimming reads\n")
  stTmp<-fastq_trim(sampTab, finalLength=finalLength, outDir="./")

  # remove original fastqs
  if(delOrig){
    delete_par(as.vector(stTmp$fname))
    system("rm *txt")
  }
 
  # Salmon
  cat("Salmon\n")
  resNames<-salmon_par(stTmp, salmonIndex=paste0(refDir,"/",salmonIndex), salmonPath=salmonPath)
  cleanup()

  sampTab<-cbind(stTmp, salmonDir=resNames)

  # transcript-level estimates
  transList<-salmon_load_tranEst(sampTab)

  # gene level estimates
  expGeneList<-gene_expr_sum(transList,geneTabfname=paste0(refDir,"/",geneTabfname))
    
  # normalize data
  expGeneList[['normalized']]<-trans_rnaseq(expGeneList[['counts']], total=total)
  expGeneList
}

############
# RNA-SEQ
############

# functions to pre-process RNA-seq data for use in CellNet
# requires salmon, trim_galore, gnu parallel to be installed

#' weighted subtraction from mapped reades
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param vector of total mapped reads per gene/transcript
#' @param total post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
downSampleW<-function
(vector,
total=1e5){ 

  totalSignal<-sum(vector)
  wAve<-vector/totalSignal
  resid<-sum(vector)-total #num to subtract from sample
  residW<-wAve*resid # amount to substract from each gene
  ans<-vector-residW
  ans[which(ans<0)]<-0
  ans
}

#' weighted subtraction from mapped reades and log applied to all
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
#' @export
trans_rnaseq<-function
(expRaw,
 total
 ){
    expCountDnW<-apply(expRaw, 2, downSampleW, total)
    log(1+expCountDnW)
  }

#' Determines read length of fastq file
#'
#' Uses awk to determine read length of first read in fastq file. Assumes all reads in file are same length.
#' @param fastq filename
#'
#' @return numeric read length
#'
#' @export
fastq_readLength<-function
(fastq
  ){
  cmd<-paste("head -n 4 ", fastq, " | awk \'{if(NR%4==2) print length($1)}\'", sep='');
  ans<-system(cmd, intern=TRUE);
  as.numeric(ans);
}

#' Trim reads to specified length
#'
#' Use cutadapt to trim reads from both ends to finalLength
#' @param sampTab sample table
#' @param fnameCol sampTab column with name of fastq files to trim
#' @param finalLength final length of each read post trimming
#' @param outDir where to store the trimmed fastq files

#' @return sampTab with filenames of trimmed fastqs appended
#'
#' @export
fastq_trim<-function### trim reads 
(sampTab, 
  fnameCol="fname",
###  cname="sra_id", not needed but still called in Clouseq processCustom
  finalLength=40,
  outDir="./")
{
 
  fnames<-paste(as.vector(sampTab[,fnameCol]), collapse=" ");
  lengths<-as.vector(sampTab$readLength);
  trimLefts<-ceiling((lengths-finalLength)/2);
  trimRights<-floor((lengths-finalLength)/2);
  trimRights<- -1 * trimRights;

  nnames<-paste(unlist(strsplit(fnames, ".fastq")), "_trimmed.fq", sep='')
  nnames<-utils_stripwhite(as.vector(nnames))
  ##cbind(sampTab, trimNames=sra_addFnames(sampTab,suffix="trimmed.fq", cname=cname)) ;
  sampTab<-cbind(sampTab, trimNames=nnames);
  tmpDF<-data.frame(fnamesIn=as.vector(sampTab[,fnameCol]), fnamesOut=nnames, lefts=trimLefts, rights=trimRights);
  tfname<-"tmpLengths.txt";
  write.table(tmpDF, file=tfname, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t");

  thecall<-"cutadapt -m 30 -u " ###20 -u -20 -o testTrimmed.fastq test1.fastq
  ###cmd<-paste("parallel --colsep \"\\t\" \"",thecall," {3} -u {4} -o ", outDir,"{2} {1}\" :::: ",tfname,sep='');
  cmd<-paste("parallel --tmpdir ./tmp/ --colsep \"\\t\" \"",thecall," {3} -u {4} -o ", outDir,"{2} {1}\" :::: ",tfname,sep='');

###  cat(cmd,"\n");
  system(cmd);
  sampTab;
}

####################################################################################
# SALMON
####################################################################################

#' quantify transcript levels using 'pseudo-alignments'
#'
#' Uses salmon
#' @param sampTab sample table
#' @param sname colname to id sample 
#' @param salmonIndex salmon index directory
#' @param cname trimNames
#' @param libraryType library type (unstranded, etc)
#' @param numThreads number of threads to use
#' @param destPrefix where should the resulting files be placed 
#' @param njobs use shell parallel too, and allot njobs
#' @param salmonPath path to salmon binary 
#'
#' @return vector of output names
salmon_par<-function
(sampTab,
  sname='sra_id',
  salmonIndex="/media/ephemeral1/dat/ref/MM_GRCh38.SalmonIndex.030816",
  cname="trimNames",
  libraryType="U",
  numThreads=5,
  destPrefix="salmonRes",
  njobs=4,
  salmonPath="~/rnaseq/SalmonBeta-0.6.1_DebianSqueeze/bin"
  ){
  
  # write a file that has the following columns:
  # {1} fastq input file name
  # {2} directory output name

  fnames<-as.vector(sampTab[,cname]);
  sras<-as.vector(sampTab[,sname]);
  oNames<-paste(destPrefix, "_", sras, sep='');
  tfname<-"tmpSalmon.txt";
  xxx<-data.frame(fastq=fnames, outnames=oNames);
  write.table(xxx, file=tfname, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t");

  # use parallel to execute salmon on multiple input files simultaneously
  cmd<-paste("parallel --tmpdir ./tmp/ --jobs ",njobs," --colsep \"\\t\" \"",salmonPath, "/salmon quant -p ", numThreads," -i ",salmonIndex, " -l \"",libraryType,"\" -r {1} -o {2}\" :::: ",tfname, sep='');
  system(cmd);
  oNames;
}

#' load Salmon-based transcript estimates
#'
#' load Salmon-based transcript estimates from quant.sf files
#' @param sampTab sample table 
#' @param cname used as convetion to name output files 
#' @param prefix text prefix to output files
#'
#' @return list of exp matricies of both NumReads and TPM
salmon_load_tranEst<-function
(sampTab, 
 cname="sra_id",
 prefix="salmonRes"
){

  # descend into each directory and load table, then merge
  tmpList<-list();
  for(i in seq(nrow(sampTab))){
    xname<-as.vector(sampTab[i,cname]);
    fname<-paste("./",prefix,"_",xname,"/quant.sf", sep='');
    cat(fname,"\n");
    #skips<-(utils_count_comments(fname)-1);

    # xname<-strsplit(as.vector(sampTab[i,"r1"]), "_R1")[[1]][1];
    tmpList[[xname]]<-read.csv(fname, sep="\t", header=1);
  }
  x<-tmpList[[xname]];
  TPM<-matrix(0, nrow=nrow(x), ncol=length(tmpList));
  colnames(TPM)<-names(tmpList);
  rownames(TPM)<-tmpList[[xname]][,1];
  NumReads<-TPM;

  for(xname in names(tmpList)){
    TPM[,xname]<-tmpList[[xname]][,'TPM'];
    NumReads[,xname]<-tmpList[[xname]][,'NumReads'];
  }
  list(TPM=TPM, NumReads=NumReads);
}

#' get files needed to run Salmon and Hisat2 QC pipeline
#'
#' get files needed to run Salmon and Hisat2 QC pipeline from S#
#' @param destination where to store them
#' @param species species. e.g. 'mouse'
#' @param bucket name of source s3 bucket
#' @param dir path on S3 bucket 
#'
#' @return nothing
trans_fetch_index<-function # get files needed to run Salmon and Hisat2 QC pipeline
(destination="/media/ephemeral1/dat/ref",
  species='mouse',
  bucket='cellnet-rnaseq',
  dir="ref"){

  curdir<-system('pwd', intern=T);
  setwd(destination);

  # fetch 4 files:
  # 1. hisat2Indices_050216.tar.gz
  fname1<-"hisat2Indices_050216.tar.gz";
  cat("fetching hisat2 indices\n");
  s3_get(dir, fname1, bucket);
  cat("unpackagin indices\n");
  utils_unpack(fname1);

  # depending on species
  
  # if mouse:
  #   2. salmon.index.mouse.050316.tgz
  #   3. geneToTrans_Mus_musculus.GRCm38.80.exo_Jun_02_2015.R 
  #   4. Mus_musculus.GRCm38.83.gtf.gz
  if(species=='mouse'){
    fname2<-"salmon.index.mouse.050316.tgz";
    fname3<-"geneToTrans_Mus_musculus.GRCm38.80.exo_Jun_02_2015.R";
    fname4<-"Mus_musculus.GRCm38.83.gtf.gz";
  }

  # if human:
  #   2. salmon.index.human.050316.tgz
  #   3. geneToTrans_Homo_sapiens.GRCh38.80.exo_Jul_04_2015.R
  #   4. Homo_sapiens.GRCh38.83.gtf.gz
  else{
    fname2<-"salmon.index.human.050316.tgz";
    fname3<-"geneToTrans_Homo_sapiens.GRCh38.80.exo_Jul_04_2015.R";
    fname4<-"Homo_sapiens.GRCh38.83.gtf.gz";
  }

  cat("fecthing and unpacking other stuff...\n");
  s3_get(dir, fname2, bucket);
  cmd<-paste("tar zxvf ", fname2, sep='');
  system(cmd);

  s3_get(dir, fname3, bucket);
  s3_get(dir, fname4, bucket);
  cmd<-paste("gzip -d ", fname4, sep='');
  system(cmd);

  setwd(curdir);
}

#' Sum transcription expression estimates to gene-level expression measures
#'
#' Sums all trnascript level expression estimates to a single gene level estimate. Needs a table that maps transcript IDs to gene IDs.
#' @param expDatList result of running salmon_load_tranEst
#' @param numCores num of cores to use for parallel 
#' @param geneTabfname gene <-> transcript mapping, df that needs cols: gene_id, transcript_id
#' @param nameCol gene ann table column name to average over
#'
#' @return list(TPM=ansTPM, counts=ansCounts);
gene_expr_sum<-function
(expDatList,
 numCores=10,
 geneTabfname="/media/ephemeral1/dat/ref/geneToTrans_Mus_musculus.GRCm38.80.exo_Jun_02_2015.R",
 nameCol="gene_name"
){
  
  matchFunc<-function(val, vect){
    which(vect==val);
  }
  
  expTPM<-expDatList[['TPM']];
  save_colnames = colnames(expTPM)
  expCounts<-expDatList[['NumReads']];
  if(!is.matrix(expTPM)){
    expTPM<-as.matrix(expTPM);
  }
  if(!is.matrix(expCounts)){
    expCounts<-as.matrix(expCounts);
  }
  
  geneTab<-utils_loadObject(geneTabfname);
  
  # keep only common variables
  sameProbes<-intersect(rownames(expTPM), rownames(geneTab));
  expTPM<-as.matrix(expTPM[sameProbes,]);
  expCounts<-as.matrix(expCounts[sameProbes,]);
  geneTab<-geneTab[sameProbes,];
  
  # all genes
  allgenes<-as.vector(geneTab[,nameCol]);
  
  # unique genes
  eids<-unique(allgenes)
  
  # make a cluster
  aClust<-parallel::makeCluster(numCores, type='SOCK')
  
  # a list of indices into all genes
  geneIndexList<-parLapply(aClust, eids, matchFunc, vect=allgenes);
  names(geneIndexList)<-eids
  uSymbols<-vector(length=length(eids));
  
  stopCluster(aClust)
  
  
  ansTPM<-matrix(0, nrow=length(eids), ncol=ncol(expTPM));
  ansCounts<-ansTPM;
  
  for(i in seq(length(geneIndexList))){
    eid<-eids[i]
    #cat(".");
    xi <-  geneIndexList[[i]];
    
    ## desProbes <- as.vector(geneTab[xi,]$probe_id);
    desProbes <- as.character(geneTab[xi,]$transcript_id);
    if(length(xi)>1){
      ansTPM[i,]<- apply(as.matrix(expTPM[desProbes,]), 2, sum);
      ansCounts[i,]<-apply(as.matrix(expCounts[desProbes,]), 2, sum);
    }
    else{
      ansTPM[i,]<-as.matrix(expTPM[desProbes,]);
      ansCounts[i,]<-as.matrix(expCounts[desProbes,]);
    }
    uSymbols[i]<-as.vector(geneTab[ xi[1] ,nameCol]);
  }
  rownames(ansTPM)<-uSymbols;
  rownames(ansCounts)<-uSymbols;
  #colnames(ansTPM)<-colnames(expTPM);
  colnames(ansTPM) = save_colnames
  #colnames(ansCounts)<-colnames(expCounts);
  colnames(ansCounts) = save_colnames
  
  list(TPM=ansTPM, counts=ansCounts);
}

if(FALSE){
#' Sum transcription expression estimates to gene-level expression measures
#'
#' Sums all trnascript level expression estimates to a single gene level estimate. Needs a table that maps transcript IDs to gene IDs.
#' @param expDatList result of running salmon_load_tranEst
#' @param numCores num of cores to use for parallel 
#' @param geneTabfname gene <-> transcript mapping, df that needs cols: gene_id, transcript_id
#' @param nameCol gene ann table column name to average over
#'
#' @return list(TPM=ansTPM, counts=ansCounts);
gene_expr_sum<-function# returns the gene summed matrix
(expDatList,
  numCores=10,
  geneTabfname="/media/ephemeral1/dat/ref/geneToTrans_Mus_musculus.GRCm38.80.exo_Jun_02_2015.R",
  nameCol="gene_name"
){

  matchFunc<-function(val, vect){
    which(vect==val);
  }

  expTPM<-expDatList[['TPM']];
  expCounts<-expDatList[['NumReads']];

  if(!is.matrix(expTPM)){
    expTPM<-as.matrix(expTPM);
  }
  if(!is.matrix(expCounts)){
    expCounts<-as.matrix(expCounts);
  }
  
  geneTab<-utils_loadObject(geneTabfname);

  # keep only common variables
  sameProbes<-intersect(rownames(expTPM), rownames(geneTab));
  expTPM<-expTPM[sameProbes,];
  expCounts<-expCounts[sameProbes,];
  geneTab<-geneTab[sameProbes,];
  
  # all genes
  allgenes<-as.vector(geneTab[,nameCol]);

  # unique genes
  eids<-unique(allgenes)

  # make a cluster
  aClust<-parallel::makeCluster(numCores, type='SOCK')

  # a list of indices into all genes
  geneIndexList<-parLapply(aClust, eids, matchFunc, vect=allgenes);
  names(geneIndexList)<-eids
  uSymbols<-vector(length=length(eids));

  stopCluster(aClust)


  ansTPM<-matrix(0, nrow=length(eids), ncol=ncol(expTPM));
  ansCounts<-ansTPM;

  for(i in seq(length(geneIndexList))){
    eid<-eids[i]
    #cat(".");
    xi <-  geneIndexList[[i]];
    
    ## desProbes <- as.vector(geneTab[xi,]$probe_id);
    desProbes <- as.character(geneTab[xi,]$transcript_id);
    if(length(xi)>1){
      ansTPM[i,]<- apply(expTPM[desProbes,], 2, sum);
      ansCounts[i,]<-apply(expCounts[desProbes,], 2, sum);
    }
    else{
      ansTPM[i,]<-expTPM[desProbes,];
      ansCounts[i,]<-expCounts[desProbes,];
    }
    uSymbols[i]<-as.vector(geneTab[ xi[1] ,nameCol]);
  }
  rownames(ansTPM)<-uSymbols;
  rownames(ansCounts)<-uSymbols;
  colnames(ansTPM)<-colnames(expTPM);
  colnames(ansCounts)<-colnames(expCounts);

  list(TPM=ansTPM, counts=ansCounts);
}}

#' re-order the sampTab, by dLevel, given dLevel names
#'
#' re-order the sampTab, by dLevel, given dLevel names
#' @param sampTab sample table
#' @param dLevel col to re-order on
#' @param nnames lelves to reorder
#'
#' @return new sample table
expr_reorder<-function# 
(sampTab,
 dLevel,
 nnames){
  
  newTab<-data.frame();
  for(nname in nnames){
    newTab<-rbind(newTab,sampTab[which(sampTab[,dLevel]==nname),]);
  }  
  newTab;    
}

#' properly read a csv file as a sampTab
#'
#' properly read a csv file as a sampTab
#' @param fname csv file
#'
#' @return data frame, should have sample_id, sample_name, description1, description2, exp_id, file_name columns
expr_readSampTab<-function
(fname
){
  sampTab<-read.csv(fname, strip.white=T, as.is=T, blank.lines.skip=T);
  # remove duplicate sample_ids
  counts<-table(sampTab$sample_id);
  if(any(counts>1)){
    xi<-names(counts)[which(counts>1)];
    ix<-names(counts)[which(counts==1)];
    sampTab<-sampTab[match(c(ix,xi), sampTab$sample_id),];
  }
  rownames(sampTab)<-sampTab$sample_id;
  sampTab;
}

#' replace the files_names column with the actual files names, and remove samples not found
#'
#' replace the files_names column with the actual files names, and remove samples not found
#' @param sampTab data.frame of sample_id, sample_name, description[1-x], exp_id, file_name
#'
#' @return data.frame, same columns as sampTab, but with corrected file_names
geo_fixNames<-function
(sampTab
){  
  ##<<note Assumes that CEL files have already been uncompressed in current directory
  # remove blank lines
  sids<-as.vector(sampTab$sample_id);
  xi<-unlist(lapply(sids, nchar));
  sampTab<-sampTab[which(xi>0),];
  if(any(colnames(sampTab)=='file_name')){
    fnames<-as.vector(sampTab$file_name);
  }
  else{
    sids<-as.vector(sampTab$sample_id);
    fnames<-sids;
    for(i in seq(length(sids))){
      fnames[i]<-paste(sids[i], ".CEL", sep='');
    }
  }
  
  cat("Matching up provided with actual file names...\n");
  # get actual filenames (all of them)
  cmd<-paste("ls *.cel");
  aFnames<-system(cmd, intern=T);
  cmd<-paste("ls *.cel*");
  aFnames<-append(aFnames,system(cmd, intern=T));
  cmd<-paste("ls *.CEL");
  aFnames<-append(aFnames,system(cmd, intern=T));
  cmd<-paste("ls *.CEL*");
  aFnames<-append(aFnames,system(cmd, intern=T));
  aFnames<-unique(aFnames);
  
  jCount<-0;
  st2<-data.frame();
  nnames<-vector();
  for(i in seq(length(fnames))){
    fname<-fnames[i];
    fnamePref<-strsplit(fnames[i], ".CEL")[[1]][1];
    ai<-grep(fnamePref, aFnames, ignore.case=T);
    if( length(ai)>0 ){
      st2<-rbind(st2, sampTab[i,]);
      jCount<-jCount+1;
      nnames<-append(nnames, aFnames[ ai[1] ]);
    }
  }
  if(any(colnames(st2)=='file_name')){
    st2$file_name<-nnames
  }
  else{
    st2<-cbind(st2, file_name=nnames);
  }
  rownames(st2)<-st2$sample_id
  st2;
  ### 
}

#' load and normalize each group separately, then return expProp (proportional expression)
#'
#' load and normalize each group separately, then return expProp (proportional expression)
#' @param sampTab sample table as returned by geo_fixNames
#' @param platform platform name. Options are ...
#' @param inc number of samples to process simultaneously
#'
#' @return raw expression expressed as a proportion of overall signal per sample
#'
#' @export
Norm_cleanPropRaw<-function
(sampTab,
 platform,
 inc=1e3
 ){
  
  expRaw<-Norm_preNormRaw(sampTab, inc=inc);
  expRaw<-exprs(expRaw);
  geneTab<-Norm_loadAnnotation(platform);
  rownames(geneTab)<-as.vector(geneTab$probe_id)
  
  ggs<-intersect(rownames(expRaw), as.vector(geneTab$probe_id));
  expRaw2<-expRaw[ggs,];
  gTab<-geneTab[ggs,];
  
  expClean<-Norm_cleanExp(expRaw2, gTab, "symbol");
  trans_dNorm(expClean);
}

#' Run justRMA on samples, nor norm, just BG
#'
#' Run justRMA on samples, nor norm, just BG
#' @param sampTab sample table 
#' @param inc number of CEL files to load at one time
#'
#' @return cleaned eSet obj
Norm_preNormRaw<-function
(sampTab,
 inc=1e3
){

  strs<-seq(from=1, to=nrow(sampTab), by=inc);
  strs<-append(strs, nrow(sampTab));
  stps<-strs[2:length(strs)];
  strs<-strs[1:(length(strs)-1)]
  xList<-list();
  for(i in seq(length(strs))){
    cat(i,"\n");
    stTmp<-sampTab[strs[i]:stps[i],];
    expTmp<-justRMA(filenames=stTmp[,'file_name'], normalize=FALSE);
    colnames(exprs(expTmp))<-as.vector(stTmp[,"sample_id"]);
    xList[[i]]<-expTmp;
  }
  
  expAll<-matrix(0, nrow=nrow(exprs(xList[[1]])), ncol=nrow(sampTab));
  for(i in seq(length(strs))){
    cat(".");
    expAll[,strs[i]:stps[i]]<-exprs(xList[[i]])
  }
  colnames(expAll)<-as.vector(sampTab[,"sample_id"]);
  rownames(expAll)<-as.vector(rownames(exprs(expTmp)));
  exprs(expTmp)<-expAll;
  expTmp;
}

#' summariz expression when there are multiple measures per gene in a matrix
#'
#' returns the gene averaged matrix
#' @param expDat exp matrix
#' @param geneTab gene annotation
#' @param nameCol how genes are uniquely IDed in geneTab
#'
#' @return expression matrix
Norm_cleanExp<-function
(expDat,
 geneTab,
 nameCol="symbol"
){

  ### re-work this to run in parallel
  if(!is.matrix(expDat)){
    expDat<-as.matrix(expDat);
  }
  
  # Make sure the ann table and expDat have same rows:
  rownames(geneTab)<-as.character(geneTab$probe_id);
  sameProbes<-intersect(rownames(expDat), as.character(geneTab$probe_id));
  expDat<-expDat[sameProbes,];
  geneTab<-geneTab[sameProbes,];
  
  eids<-unique(as.vector(geneTab[,nameCol]));
  uSymbols<-vector(length=length(eids));
  ans<-matrix(nrow=length(eids), ncol=ncol(expDat));
  for(i in seq(length(eids))){
    eid<-eids[i];
    xi <-  which( geneTab[,nameCol]==eid );
    desProbes <- as.character(geneTab[xi,]$probe_id);
    if(length(xi)>1){
      ans[i,]<- apply(expDat[desProbes,], 2, mean);
    }
    else{
      ans[i,]<-expDat[desProbes,];
    }
    uSymbols[i]<-as.vector(geneTab[ xi[1] ,nameCol]);
  }
  rownames(ans)<-uSymbols;
  colnames(ans)<-colnames(expDat);
  ans;
}

#' load the specified gene annotation table
#'
#' load the specified gene annotation table
#' @param pName
#'
#' @return table of probe_id, entrez id and gene symbols
Norm_loadAnnotation<-function
(pName
){
  
  # Need to add code to specifically download and install proper library versions.
  geneTab<-'';
  
  if(pName=="mogene10stv1"){
    pName<-"mogene10sttranscriptcluster";
  }
  if(pName=="hugene10stv1"){
    pName<-"hugene10sttranscriptcluster";
  }
  libName<-paste(pName, ".db",sep='')
  x<-require(libName, character.only=TRUE);
  
  if(FALSE){
  if(!x){
    source("http://bioconductor.org/biocLite.R")
    biocLite(libName);
    x<-require(libName, character.only=TRUE);
  }
  }
  if(!x){
    cat(".loadAnnotation\tUnable to install ",libName,"\n")
    return;
  }
  else{
    entrezCmd<-paste(pName, "ENTREZID", sep='');
    idType<-'entrezgeneid';
    probes<-eval(parse(text=entrezCmd));
    probeTable<-links(probes);
    
    symbolsCmd<-paste(pName,"SYMBOL", sep='');
    symbols<-eval(parse(text=symbolsCmd));
    symbols<-links(symbols);
    
    geneTab<-merge(probeTable, symbols);
  }
  geneTab;
  
}


#########################################
# transformations and normalizations
#########################################

#' express as a fraction of the column total
#'
#' express as a fraction of the column total
#' @param expDat expression matrix
#' @param xFact scale by this value
#'
#' @return transformed data matrix
trans_dNorm<-function 
(expDat,
 xFact=1e5
){
  ans<-matrix(0, nrow=nrow(expDat), ncol=ncol(expDat));
  for(i in seq(ncol(expDat))){
    ans[,i]<-expDat[,i]/sum(expDat[,i]);    
  }
  ans<-ans*xFact;
  colnames(ans)<-colnames(expDat);
  rownames(ans)<-rownames(expDat);
  ans;
}

#' shift the expression of each sample so that the min=minVal
#'
#' shift the expression of each sample so that the min=minVal
#' @param expDat expression matrix
#' @param minVal
#'
#' @return transformed expression matrix
trans_eShift<-function #
(expDat,
 minVal=0){
  
  myFunc<-function(vect, myMin){
    vect-min(vect)+myMin;
  }
  
  ans<-apply(expDat, 2, myFunc, myMin=minVal);
  rownames(ans)<-rownames(expDat);
  colnames(ans)<-colnames(expDat);
  ans;
}

#' quantile normalize expression matrix
#'
#' equalize the distributions of columns
#' @param expDat
#'
#' @return quantile normalized matrix
Norm_quantNorm<-function#
(expDat){
  require(preprocessCore);
  if(is.data.frame(expDat)){
    expDat<-as.matrix(expDat);
  }
  ans<-normalize.quantiles(expDat);
  colnames(ans)<-colnames(expDat);
  rownames(ans)<-rownames(expDat);
  ans;
}

##################################
#
# functions to 
# ... fetch fastqs (if needed)
# ... process data locally
#uu
##################################

#' Fetch fastq files from S3
#'
#' It assumes that the files are publicly accessible. If not then you should use s3_get after you have run aws configure and enter your aws credentials
#' @param bucket name of bucket
#' @param path path to fastqs
#' @param sampTab sample table with a column fname listing the fastq files
#' @param fnameCol name of column listing fastq files
#'
#' @return sammple table with de-compressed file names replacing old file names
#' @export
cn_s3_fetchFastq<-function
(bucket,
 path,
 sampTab,
 fnameCol="fname",
 compressed=NA)
{
  pref<-paste0("https://s3.amazonaws.com/", bucket,"/",path,"/");
  fnames<-as.vector(sampTab[,fnameCol])
  nfnames<-vector()
  for(fname in fnames){
    destName<-fname
    target<-paste0(pref, fname)
    download.file(target, destfile=destName)
    if(!is.na(compressed)){
      if(compressed=="gz"){
        cmd<-paste0("gzip -d ", destName)
        nfile<-strsplit(destName, ".gz")[[1]][1]
     }
      if(compressed=="bz2"){
        cmd<-paste0("bzip2 -d ", destName)
        nfile<-strsplit(destName, ".bz2")[[1]][1]
      }
      system(cmd)
    }
    else{
      nfile<-destName
    }
    nfnames<-append(nfnames, nfile)
  }
  sampTab[,fnameCol]<-nfnames
  sampTab
}



########################################################
#
# SRA
#
########################################################

#' processes a SRA study
#'
#' This function is a wrapper to qcAndSalmon. It also summarizes transcript into gene expression estimates
#' @param sampTab sample table
#' @param inc increment size
#' @param studyID study id
#' @param bucket S3 bucket where to put the results
#' @param startDir S3 bucket path to store results
#' @param finalLength final length of reads 
#' @param subProp proportion of reads to sample from for QC
#' @param target species/genome for expression estimates
#' @param gtfFile genomic feature of target genome for QC
#' @param geneTabfname transcript <-> gene data frame
#' @param salmonIndex path to salmon index for expression estimation
#'
#' @return sample table with QC measures appened. side effect is storage of expression estimates + sample table to S3
#'
#' @examples
#' processSRA(stX, inc=5)
#'
#' @export
processSRA<-function
(sampTab,
  inc=5,
  studyID = NA,
  bucket="cellnet-rnaseq",
  startDir='singlecell',
  finalLength=40,
  subProp=0.001,
  target="mouse",
  gtfFile="/media/ephemeral1/dat/ref/Mus_musculus.GRCm38.83.gtf",
  geneTabfname="/media/ephemeral1/dat/ref/geneToTrans_Mus_musculus.GRCm38.80.exo_Jun_02_2015.R",
  salmonIndex="/media/ephemeral1/dat/ref/MM_GRCh38.SalmonIndex.030816"
){

  if(is.na(studyID)) {
        studyID<-as.vector(sampTab$study_id)[1];
  } else {
    studyID = as.vector(studyID)
  }
  print(studyID[[1]])

  stList<-sra_breakTable(sampTab, inc=inc);
  ansList<-list();
  for(i in seq(length(stList))){
    stTmp<-stList[[i]];
    ansList[[i]]<-qcAndSalmon(stTmp, studyID, bucket=bucket, finalLength=finalLength, subProp=subProp, target=target, gtfFile=gtfFile, salmonIndex=salmonIndex)
  }

  newSampTab<-ansList[[1]];
  if(length(ansList)>1){
    for(i in 2:length(ansList)){
      newSampTab<-rbind(newSampTab, ansList[[i]]);
    }
  }

  # save and send sample table
  bPath<-paste(startDir,"/",studyID,"/results",sep='');
  fname<-paste("sampTab_", studyID, "_QC.R", sep='');
  save(newSampTab, file=fname);
  s3_put(bPath, fname, bucket);

  # transcript-level estimates
  transList<-salmon_load_tranEst(newSampTab);
  fname<-paste("expTransList_", studyID, ".R", sep='');
  save(transList, file=fname);
  s3_put(bPath, fname, bucket);

  # gene level estimates
  expGeneList<-gene_expr_sum(transList,geneTabfname=geneTabfname);
  fname<-paste("expGeneList_", studyID, ".R", sep='');
  save(expGeneList, file=fname);
  s3_put(bPath, fname, bucket);
  newSampTab;
}

#' Runs the qc and salmon pipeline used by processSRA
#'
#' Fetch raw sra files named in sampTab and convert to fastq (using SRA toolkit), trim reads to standard length, run QC by aligning sub-sampled fastq files aginst several genomes (human, mouse, zebrafish, fly, ecoli, yeast), and counting hits to genomic features of the target genome. Finish by estimating gene expression levels with Salmon. 
#' @param sampTab sample table
#' @param studyID study id
#' @param bucket S3 bucket where to put the results
#' @param finalLength final length of reads 
#' @param subProp proportion of reads to sample from for QC
#' @param target species/genome for expression estimates
#' @param gtfFile genomic feature of target genome for QC
#' @param salmonIndex path to salmon index for expression estimation
#'
#' @return sample table with QC measures appened
qcAndSalmon<-function
(sampTab,
  studyID,
  bucket="cellnet-rnaseq",
  finalLength=40,
  subProp=0.001,
  target="mouse",
  gtfFile="/media/ephemeral1/dat/ref/Mus_musculus.GRCm38.83.gtf",
  salmonIndex="/media/ephemeral1/dat/ref/MM_GRCh38.SalmonIndex.030816"
  ){

  cat("Fetching and converting SRAs ___________________________________\n");
  stTmp<-fetchAndConvert(sampTab);
  cat("Done ___________________________________________________________\n");

  study<-studyID
  bPath<-paste("trainingdata/", study,"/fastq",sep='');
  cat("Packing raw data and sending to S3 _____________________________\n");
  bundleAndSend(bucket, bPath, as.vector(stTmp$fname));
  cat("Done ___________________________________________________________\n");

  cat("Trimming reads, sending to S3 __________________________________\n");
  stTmp<-fastq_trim(stTmp, finalLength=finalLength, outDir="./");
  bPath<-paste("trainingdata/", study,"/trimmed",sep='');
  bundleAndSend(bucket, bPath, as.vector(stTmp$trimNames));
  cat("Done ___________________________________________________________\n");

  # remove original fastqs
  delete_par(as.vector(stTmp$fname));
  system("rm *txt");

  #QC
  stTmp<-qc_analysis(stTmp,
    subProp=subProp,
    target=target,
    gtfFile=gtfFile,
    cname="trimNames");

  # Salmon
  resNames<-salmon_par(stTmp, salmonIndex=salmonIndex);
  cleanup()

  cbind(stTmp, salmonDir=resNames);

}

####################################################################################
# SRA functions
####################################################################################

#' remove txt and trimmed.fq after running salmon
#'
#' remove txt and trimmed.fq after running salmon
#'
#' @return nothing
cleanup<-function(){
  cmd<-paste("rm *txt");
  system(cmd);
  cmd<-"rm *trimmed.fq";
  system(cmd);
}

#' fetch sras and converts them to fastqs, can take a REALLY long time and is unreliable
#'
#' fetch sras and converts them to fastqs, can take a REALLY long time and is unreliable. SRA toolkit must be installed and configured. Make sure the store the SRAs on a drive with lots of space. Must be consistent with tmpPath.
#' @param sampTab sample table
#' @param tmpPath where sras are stored
#'
#' @return sample table with read length appended
fetchAndConvert<-function# does just that
(sampTab,
 tmpPath="/media/ephemeral0/dat"){
  sraids<-as.vector(sampTab$sra_id);
  # puts sra files inot /media/ephemeral0/dat/sra
  sra_fetch(sraids);
  # converts the sra files to fastqs and puts them in the current directory, and deletes the sra files
  sra_parse(sraids, outdir=tmpPath);

  # I am only dealing with single end reads, or one end of paired end reads
  sra_delete("./", "_2.fastq")

  # update the sample table
  sampTab<-cbind(sampTab, fname=sra_addFnames(sampTab));
  sampTab<-cbind(sampTab, readLength=unlist(lapply(sampTab$fname, fastq_readLength)));
  sampTab;
}

#' fetch sras to path defined when configuring sra tools.
#'
#' fetch sras to path defined when configuring sra tools. does it in parallel. for now limited to files <=600GB.
#' @param sra_ids sra accessions to fetch
#'
#' @return nothing
sra_fetch<-function# 
(sra_ids){
  cmd<-paste("parallel prefetch -X 600G -v ::: ", paste(sra_ids, collapse=" "))
  system(cmd);
}

#' unpacks the SRA files, splits them, and deletes SRAs
#'
#' unpacks the SRA files, splits them, and deletes SRAs
#' @param sra_ids sra ids
#' @param delete whether or not to delete the originals after fastq-dump
#' @param outdir where to store resulting fastw files
#' @param indir where are the sras stored by sra-toolkit
#' @return nothing
sra_parse<-function# 
(sra_ids,
 delete=TRUE,
 outdir="/media/ephemeral0/dat",
 indir="/media/ephemeral1/dat/seq/sra"){

  ### assumes that sra filename is: sraid.sra  
  curdir<-system('pwd', intern=T);
  setwd(indir);
  fnames<-paste(sra_ids, ".sra", sep='');
  fnames<-paste(fnames, collapse=" ");
  cmd<-paste("parallel fastq-dump -Q 33 -M 25 --skip-technical --split-3 --O ",outdir," ::: ",fnames, sep='');
  cat("parsing files...\n")
  system(cmd);
  
  if(delete){
    cmd<-paste("rm *sra");
    cat("deleting sras...\n")
    system(cmd)
  }
  setwd(curdir);
}

#' delete files in parallel
#'
#' used in qcAndSalmon
#' @param fnames files to delete
#'
#' @return nothing
delete_par<-function#
(fnames){
  cmd<-paste("parallel rm ::: ", paste(fnames, collapse=" "));
  cat(cmd,"\n");
  system(cmd);
}

#' Deletes files with specified suffix in specified directory
#'
#' Deletes files with specified suffix in specified directory. Used in fetchAndConvert
#' @param xdir path to files
#' @param suff what types of files to rm
#'
#' @return nothing
sra_delete<-function
(xdir="/media/ephemeral0/dat/seq/tmpOut",
  suff="gz"){
  cmd<-paste("rm ",xdir,"/*",suff,sep='');
  system(cmd);
}

#' break a sample table into _inc_ size pieces
#'
#' break a sample table into _inc_ size pieces. used in processSRA.
#' @param sampTab sample table
#' @param inc increment size
#'
#' @return list of sample tables
sra_breakTable<-function# 
(sampTab, 
 inc=6
 ){
  ans<-list();
  
  if(nrow(sampTab)<inc){
    ans[[1]]<-sampTab;
  }
  else{
    xi<-1;
    str<-1;
    stp<-inc;
    while(str<nrow(sampTab)){
      
      # I don't want there to be any 1 row tables, which can only happen on last step
      # if it's a multi row table to begin with
      tmpStp<-stp+inc;
      if(tmpStp>nrow(sampTab)){
        tmpStp<-nrow(sampTab)
      }
      tmpStr<-stp+1;
      if( tmpStp==tmpStr ){
        stp<-nrow(sampTab);
        
      }
      ans[[xi]]<-sampTab[str:stp,];
      xi<-xi+1;
      tmpStp<-stp+inc;
      if(tmpStp>nrow(sampTab)){
        tmpStp<-nrow(sampTab)
      }
      str<-stp+1;
      stp<-tmpStp;
    }
  }
  ans;
}

#######################################################
#
# functions for running the QC pipeline
#
#######################################################

#' Run qc pipeline
#'
#' subsamples reads, aligns to each genome, and counts the reads mapping to features in the target genome
#' @param sampTab sample table
#' @param sname column of default prefix
#' @param subProp fraction of reads to sample
#' @param target species/genome
#' @param gtfFile gtf path/filename
#' @param cname column name containing input file names
#'
#' @return sample tablein which qc measures have been appended
qc_analysis<-function #wrapper to do it all and clean up
(sampTab,
  sname="sra_id",
  subProp=0.001,
  target="mouse",
  gtfFile="/media/ephemeral1/dat/ref/Mus_musculus.GRCm38.83.gtf",
  cname="trimNames"
  ){

  cat("sub-sampling\n")
  # sub-sample reads
  subsamp_fastqs(as.vector(sampTab[,cname]), prop=subProp)

  # align
  indexTab<-hisat_findIndices()
  cat("aligning\n");
  newTab<-hisat_QC_aligns(sampTab, indexTab, sname=sname, target=target, cname=cname);
  cat("htseq\n")
  newTab<-feature_count_par(newTab, gtfFile, sname=sname, target=target);

  # clean up
  cmd<-paste("rm subset_*");
  system(cmd);
  system("rm *.SAM");
  system("rm *.bam");
  system("rm *.htseq.ann");
  newTab;
}

#' randomly samples reads from single end fastqs
#'
#' randomly samples reads from single end fastqs uses python script in parallel
#' @param fnames files to subsample
#' @param path to python script subsetOne.py
#' @param proporation of reads to select
#'
#' @return nothing
subsamp_fastqs<-function# 
(fnames,
  codeDir="~/code/",
  prop=0.01){

  tmpfname<-"tmpFile.txt";
  write.table(fnames, file=tmpfname, quote=FALSE, col.names=FALSE, row.names=FALSE);
  cmd<-paste("cat ", tmpfname, " | parallel --tmpdir ./tmp/ ",codeDir,"subsetOne.py ", prop, " {} subset_{}",sep='');
  system(cmd);

}

#' finds hisat indecies
#'
#' finds hisat indecies
#' @param dir base path to indecies
#'
#' @return df of nickname -> path
hisat_findIndices<-function
(dir="/media/ephemeral1/dat/ref/hisat2Indices"){
  cmd<-paste("ls", dir);
  nnames<-system(cmd, intern=TRUE);
  paths<-paste(dir, "/", nnames, sep='');
  inames<-vector();

  curdir<-system('pwd', intern=T);

  for(i in seq(length(paths))){
    setwd(paths[i]);
    cmd<-paste("ls *.ht2",sep='');
    x<-system(cmd, intern=T);
    paths[i]<-paste(paths[i],"/",strsplit(x[1], ".1.ht2")[[1]][1], sep='');
    setwd(curdir);
  }
  setwd(curdir);
  data.frame(paths=paths, nickname=nnames)
}

#' wrapper to hisat_QC_align
#'
#' wrapper to hisat_QC_align
#' @param sampTab sample table
#' @param indexDF index data frame result of running hisat_findIndices
#' @param sname col name
#' @param target species/genome
#' @param cname colname indicating files to align
#' @param hisatDir path to hisat 
#'
#' @return sampTab with hisat QC stats added AND features counts vs target genome
hisat_QC_aligns<-function### returns 
(sampTab,
  indexDF,
  sname="sra_id",
  target='mouse',
  cname="trimNames",# 
  hisatDir="~/rnaseq/hisat2-2.0.3-beta"){

  tnames<-as.vector(sampTab[,cname]);
  alirates<-matrix(0,nrow=length(tnames), ncol=nrow(indexDF));

  for(i in seq(length(tnames))){
    fastq<-paste("subset_", tnames[i], sep=''); 
    cat(fastq,"\n");
    tmpX<-hisat_QC_align(fastq, indexDF, as.vector(sampTab[i,sname]), target);
    alirates[i,]<-tmpX[,"ali.rate"]; #build in code to extract feature count and delete unneeded files

  }
  colnames(alirates)<-as.vector(tmpX[,'genome']);
  sampTab<-cbind(sampTab, alirates);
}

#' performs qc on fastq file
#'
#' aligns fastq to genomes indexed and pointed to by indexDF, counts features aligning to target
#' @param fastq file name
#' @param indexDF index data frame result of running hisat_findIndices
#' @param sname col name
#' @param target species/genome
#' @param hisatDir path to hisat 
#'
#' @return sampTab with hisat QC stats added AND features counts vs target genome
hisat_QC_align<-function
(fastq,
  indexDF,
  sname,
  target='mouse',
  hisatDir="~/rnaseq/hisat2-2.0.3-beta"){
  tfname<-paste(sname,"_tmpfile.txt", sep='');

  thecall<-paste(hisatDir, "/hisat2", sep='')

  outnames<-paste(sname, "_", as.vector(indexDF[,2]), sep='');
  samNames<-paste(outnames, ".SAM", sep='');
  summaries<-paste(outnames, ".summary", sep='');
  indexDF<-cbind(indexDF, sams=samNames);
  indexDF<-cbind(indexDF, summ=summaries);

  write.table(indexDF, file=tfname, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t");

  cmd<-paste("parallel --colsep \"\\t\" \"",thecall," -x {1} -U ",fastq," -S {3} 2> {4}\" :::: ",tfname,sep='');
  system(cmd);


  system(paste("rm ",tfname,sep=''));

  alis<-vector();
  nnames<-as.vector(indexDF$nickname);
  for(nname in nnames){

    cmd<-paste("grep overall ",sname,"_",nname,".summary | cut -f1 -d\" \"", sep='');
    tmp<-system(cmd, intern=TRUE);
    tmp<-strsplit(tmp, "%")[[1]][1];
    alis<-append(alis, tmp);
    if(target!=nname){
      cmd<-paste("rm ", sname,"_",nname,".SAM", sep='');
      system(cmd);
    }
    cmd<-paste("rm ",sname,"_",nname,".summary", sep='');
    system(cmd);
  }
  data.frame(genome=nnames, ali.rate=as.numeric(alis));
}

#' count the number of reads that overlap with genomic features
#'
#' count the number of reads that overlap with genomic features uses samtools/htseq-count
#' @param sampTab sample table
#' @param gtfpath gtf of genomic features
#' @param sname column to use to name output files
#' @param target species, e.g. 'mouse'
#'
#' @return sample table with hit rates appended
feature_count_par<-function
(sampTab,
 gtfpath,
 sname="sra_id",
 target='mouse'){

  sras<-as.vector(sampTab[,sname]);
  prefix<-paste(sras, "_", target,sep='');
  
  samFiles  <- paste(prefix, ".SAM", sep='');
  bamFiles  <- paste(prefix, ".bam", sep='');
  sbamFiles <- paste(prefix,".sorted.bam", sep='');
  countFiles<- paste(prefix,".htseq.ann", sep='');

  tfname<-paste("tmpCount.txt");
  x<-data.frame(sam=samFiles, bam=bamFiles, sbam=sbamFiles, out=countFiles)
  write.table(x, tfname, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  # convert to BAM
  cmd<-paste("parallel --colsep \"\\t\" \"samtools view -bS {1} > {2}\" :::: ",tfname,sep='');
  system(cmd);
  ###cat(cmd, "\n");

  # sort them
  cmd<-paste("parallel --colsep \"\\t\" \"samtools sort -n {2} -o {3}\" :::: ",tfname,sep='');
  system(cmd);
  ###cat(cmd, "\n");

  # count reads that overlap features in GTF
  cmd<-paste("parallel --colsep \"\\t\" \"htseq-count -r name -f bam -s no -i gene_biotype {3} ",gtfpath, " > {4}\" :::: ",tfname,sep='');
  ###cat(cmd,"\n");
  system(cmd);

  ansList<-list();
  for(i in seq(length(countFiles))){
    oname<-countFiles[i];
    resDF<-read.csv(oname, sep="\t", as.is=T, strip.white=TRUE);
    ansList[[i]]<-htseqCount(resDF);
  }
  cat(nrow(ansList[[1]]),"\n");
  hitrates<-matrix(0,nrow=length(sras), ncol=nrow(ansList[[1]]));
  for(i in seq(length(sras))){
    hitrates[i,]<-as.numeric(ansList[[i]][,"prop"]);
  }
  colnames(hitrates)<-ansList[[1]][,"feature"]
  #rownames(hitrates)<-sras;
  cbind(sampTab, hitrates)

}

#' compiles read count hits
#'
#' compiles read count hits
#' @param df data.frame that results from htseq-count
#' @param promNames vector of features names of interest
#'
#' @return data.frame of feature, proportion of reads mapping to this feature
htseqCount<-function###
(df,
  promNames=c("rRNA","Mt_rRNA","protein_coding","__no_feature"," __ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")){

# write out the total number of reads/piared reads
# and the percentage falling into each category in promNames

  props<-vector();
  total<-sum(df[,2]);
  ansVect<-vector();
  pnames<-vector();
  for(i in seq(length(promNames))){
    promName<-promNames[i];
    promName<-utils_stripwhite(promName);
    tmp<-as.vector(df[df[,1]==promName,2]);
    ansVect<-append(ansVect, tmp/total);
    pnames<-append(pnames, promName)
  }

  # remove annoying "__" from some feature names
  pnames<-gsub("__", "", pnames);
  data.frame(feature=pnames, prop=ansVect);
}


