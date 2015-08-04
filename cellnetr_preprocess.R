

expr_readSampTab<-function
### properly read a csv file as a sampTab
(fname
 ### csv file
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
  ### data frame, should have sample_id, sample_name, description1, description2, exp_id, file_name columns
}


geo_fixNames<-function
### replace the files_names column with the actual files names, and remove samples not found
(sampTab
 ### data.frame of sample_id, sample_name, description[1-x], exp_id, file_name
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
  ### data.frame, same columns as sampTab, but with corrected file_names
}

Norm_cleanPropRaw<-function
### load and normalize each group separately, then return expProp (proportional expression)
(sampTab,
 ### sample table as returned by geo_fixNames
 platform,
 ### platform name. Options are ...
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
  ### raw expression expressed as a proportion of overall signal per sample
}

Norm_preNormRaw<-function
### Run justRMA on samples, nor norm, just BG
(sampTab,
 ### sample table
 inc=1e3
 ### number of CEL files to load at one time
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
  ###   cleaned eSet obj
}

Norm_cleanExp<-function
### returns the gene averaged matrix
(expDat,
 ### exp matrix
 geneTab,
 ### gene annotation
 nameCol="symbol"
 ### gene ann table column name to average over
){
  if(!is.matrix(expDat)){
    expDat<-as.matrix(expDat);
  }
  
  # Make sure the ann table and expDat have same rows:
  ## altered 05-03-13 
  ## rownames(geneTab)<-as.vector(geneTab$probe_id);
  ## sameProbes<-intersect(rownames(expDat), as.vector(geneTab$probe_id));
  rownames(geneTab)<-as.character(geneTab$probe_id);
  
  sameProbes<-intersect(rownames(expDat), as.character(geneTab$probe_id));
  
  expDat<-expDat[sameProbes,];
  ## cat(length(sameProbes),"\n");
  geneTab<-geneTab[sameProbes,];
  
  ##    cat(length(sameProbes),"\n");
  
  eids<-unique(as.vector(geneTab[,nameCol]));
  uSymbols<-vector(length=length(eids));
  ##  cat(length(eids),"\n");
  ans<-matrix(nrow=length(eids), ncol=ncol(expDat));
  for(i in seq(length(eids))){
    eid<-eids[i];
    #cat(".");
    xi <-  which( geneTab[,nameCol]==eid );
    
    ## desProbes <- as.vector(geneTab[xi,]$probe_id);
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
  ### gene averaged matrix
}


Norm_loadAnnotation<-function
### load the specified gene annotation table
(pName
 ### platform name
){
  
  # Need to add code to specifically download and install proper library versions.
  geneTab<-'';
  
  if(pName=="mogene10stv1"){
    pName<-"mogene10sttranscriptcluster";
  }
  if(pName=="hugene10stv1"){
    pName<-"hugene10sttranscriptcluster";
  }
  libName<-paste(pName, ".db",sep='');
  ## cat("Loading ", libName,"\n");
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
  ### table of probe_id, entrez id and gene symbols
}

trans_dNorm<-function 
### express as a fraction of the column total
(expDat,
 ### expression matrix
 xFact=1e5
 ### scale by this value
){
  ans<-matrix(0, nrow=nrow(expDat), ncol=ncol(expDat));
  for(i in seq(ncol(expDat))){
    ans[,i]<-expDat[,i]/sum(expDat[,i]);    
  }
  ans<-ans*xFact;
  colnames(ans)<-colnames(expDat);
  rownames(ans)<-rownames(expDat);
  ans;
  ### transformed data matrix
}


trans_eShift<-function # shift the expression of each sample so that the min=minVal
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

find_tfs<-function# find transcript factors
(species='Hs' # species abbreviation
  ){

  cat("Loading gene annotations ...\n")
  require(GO.db);

  if(species=='Hs'){
    require(org.Hs.eg.db);
    egSymbols<-as.list(org.Hs.egSYMBOL);
    goegs<-as.list(org.Hs.egGO2ALLEGS);
  }
  else{
    require(org.Mm.eg.db);
    egSymbols<-as.list(org.Mm.egSYMBOL);
    goegs<-as.list(org.Mm.egGO2ALLEGS);
  }

  goterms<-as.list(GOTERM);
  goids<-names(goegs);
  onts<-lapply(goids, Ontology);
  bps<-onts[onts=='BP'];
  goids<-names(unlist(bps));

  cat("matching gene symbols and annotations")
  gobpList<-list();
  for(goid in goids){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gobpList[[goterm]]<-genes;
  }

  ### newHsTRs<-gobpList[['regulation of transcription, DNA-dependent']];
  regNames<-names(gobpList)[grep("regulation of transcription", names(gobpList))];
  trs<- unique(unlist(gobpList[regNames]));
  cat("Regulation of transcription: ", length(trs),"\n");

  mfs<-onts[onts=='MF'];
  goidsMF<-names(unlist(mfs));

  gomfList<-list();
  for(goid in goidsMF){
    egs <- goegs[[ goid ]];
    goterm<-Term(goterms[[goid]]);
    genes<-sort(unique(as.vector(unlist( egSymbols[egs] ))));
    gomfList[[goterm]]<-genes;
  }
  dbs<-gomfList[['DNA binding']];
  cat("DNA binding: ", length(dbs),"\n");
  sort(intersect(trs, dbs));
}


