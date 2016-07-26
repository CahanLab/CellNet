# (C) Patrick Cahan 2014-2016
# functions to pre-process RNA-seq data for use in CellNet
# requires salmon, trim_galore, gnu parallel to be installed

downSampleW<-function(vector, total=1e5){ # function to simulate expression profile of  _total_ mapped reads

   totalSignal<-sum(vector);
   wAve<-vector/totalSignal;
     resid<-sum(vector)-total; ### num to subtract from sample
     residW<-wAve*resid; #amount to substract from each gene
     ans<-vector-residW;
     ans[which(ans<0)]<-0;
     ans;
}

fastq_readLength<-function# returns read length assumes all reads in file are same length
(fastq### filename
  ){
  cmd<-paste("head -n 4 ", fastq, " | awk \'{if(NR%4==2) print length($1)}\'", sep='');
  ans<-system(cmd, intern=TRUE);
  as.numeric(ans);
}

fastq_trim<-function### trim reads 
(sampTab, 
  fnameCol="fname1",
  cname="sra_id",
  finalLength=40,
  outDir="./",
  galorePath="~/rnaseq")
{

  fnames<-paste(as.vector(sampTab[,fnameCol]), collapse=" ");
  lengths<-as.vector(sampTab$readLength);
  trimLefts<-ceiling((lengths-finalLength)/2);
  trimRights<-floor((lengths-finalLength)/2);
  tmpDF<-data.frame(fnames=as.vector(sampTab[,fnameCol]), lefts=trimLefts, rights=trimRights);
  tfname<-"tmpLengths.txt";
  write.table(tmpDF, file=tfname, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t");
  thecall<-paste(galorePath, "/trim_galore --no_report_file --clip_R1 ",sep='');
  cmd<-paste("parallel --colsep \"\\t\" \"",thecall," {2} --three_prime_clip_R1 {3} --length 20 -o ", outDir," {1}\" :::: ",tfname,sep='');
  system(cmd);

  nnames<-paste(unlist(strsplit(fnames, ".fastq")), "_trimmed.fq", sep='')
  nnames<-utils_stripwhite(as.vector(nnames))
  cbind(sampTab, trimNames=nnames);
}



####################################################################################
# SALMON
####################################################################################


salmon_par<-function# quantify transcript levels using 'pseudo-alignments'
(sampTab,     ### sample table
  sname='sra_id',
 salmonIndex="/media/ephemeral1/dat/ref/MM_GRCh38.SalmonIndex.030816", ### index directory
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
  cmd<-paste("parallel --jobs ",njobs," --colsep \"\\t\" \"",salmonPath, "/salmon quant -p ", numThreads," -i ",salmonIndex, " -l \"",libraryType,"\" -r {1} -o {2}\" :::: ",tfname, sep='');
  system(cmd);
  oNames;
}


salmon_load_tranEst<-function# load Salmon-based transcript estimates
(sampTab, 
 cname="sra_id",
 prefix="salmonRes"
){
  
  # retuns a list of exp matricies of both NumReads and TPM

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

gene_expr_sum<-function# returns the gene summed matrix
(expDatList,
  numCores=10,
 ### list of exp matrix
 geneTabfname="/media/ephemeral1/dat/ref/geneToTrans_Mus_musculus.GRCm38.80.exo_Jun_02_2015.R",
 ### gene annotation, needs cols: gene_id, transcript_id
 nameCol="gene_name"
 ### gene ann table column name to average over
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
  aClust<-makeCluster(numCores);

  # a list of indices into all genes
  geneIndexList<-parLapply(aClust, eids, matchFunc, vect=allgenes);
  names(geneIndexList)<-eids
  uSymbols<-vector(length=length(eids));

  stopCluster(aClust)

  # which of these need to be summed? -- not implemented yet
  ##llengths<-unlist(lapply(geneIndexList, length));
  ##multiTransGenes<-eids[which(llengths)>1];
  ##TransGenes<-setdiff(eids, multiTransGenes);

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

  
  ### transcript summed matrix
    list(TPM=ansTPM, counts=ansCounts);
}
