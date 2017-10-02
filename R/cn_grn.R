# CellNet
# (C) Patrick Cahan 2012-2016
# GRN reconstruction functions


#' make a CLR-like matrix for assessment purposes
#'
#' computes gene-gene correlation, then clr-like zscores
#'
#' @param samptTab sample table
#' @param expDat properly normalized expression matrix
#' @param tfs vector of transcription factors
#' @param grnSampSize min number of samples per ct or group 
#' @param normDat whether to quantile normalize the data prior to computing correlations
#' @param dLevel column name in sample table to group samples by
#'
#' @return zscore matrix
#' 
#' @export
#'
cn_clr<-function
(sampTab,
 expDat,
 species="Mm",
 tfs=NA,
 grnSampSize=40, #min number of samples per ct or group
 normDat=FALSE,
 dLevel='description1')
{
  targetGenes<-rownames(expDat)
  if(is.na(tfs)){
    cat("Defining transcriptional regulators...\n");
    tfs<-find_tfs(species)
  }
  tfs<-intersect(targetGenes, tfs)

  stGRN<-sample_profiles_grn(sampTab, minNum=grnSampSize)
  cat("Number of samples per CT: ",mean(table(stGRN[,dLevel])),"\n")
  expGRN<-expDat[,rownames(stGRN)]
  if(normDat){
    cat("Normalizing expression data...\n")
    expGRN<-Norm_quantNorm(expGRN)
  }

  cat("Calculating correlation...\n");
  corrs<-grn_corr_round(expGRN)
  cat("Calculating context dependent zscores...\n");
  grn_zscores(corrs, tfs)
}



#' Reconstruct CT-specific GRNs
#'
#' runs getRawGRN, findSpecGenes, and specGRNs
#'
#' @param samptTab sample table
#' @param expDat properly normalized expression matrix
#' @param tfs vector of transcription factors
#' @param grnSampSize min number of samples per ct or group 
#' @param dLevel column name in sample table to group samples by
#' @param normDat whether to quantile normalize the data prior to computing correlations
#' @param corrs matrix of gene-gene pearson correlations
#' @param zscores CLR-defined context dependent zscores
#' @param cval template matching threshold for overall CT specific expression
#' @param cvalGK template matching threshold for developmentally shared CT specific expression
#' @param zthresh GRN zscore threshold at which to remove weak edges
#' @param holmSpec pvalue threshold for template matching
#'
#' @return list of overallGRN (desc), specGenes (desc), ctGRNs (desc)
#'
#' @examples
#' system.time(grnAll<-cn_make_grn(stAll, expAll, species='Mm', zThresh=6) )
#' 
#' @export
#'
cn_make_grn<-function
(sampTab,
 expDat,
 species="Mm",
 tfs=NA,
 ###snName, ### network name prefix
 grnSampSize=0, #min number of samples per ct or group
 normDat=FALSE,
 corrs=NA, ### pearson correlation values
 zscores=NA, ### zscores
 cval=0.5,
 cvalGK=0.75,
 dLevel='description1',
 dLevelGK="description6",
 zThresh=4,
 holmSpec=1e-6,
 prune=FALSE)
{
  targetGenes<-rownames(expDat)
  if(is.na(tfs)){
    cat("Defining transcriptional regulators...\n");
    tfs<-find_tfs(species)
  }
  tfs<-intersect(targetGenes, tfs)

  if(grnSampSize==0){
    grnSampSize<-min(table(sampTab[,dLevel]))
  }

  stGRN<-sample_profiles_grn(sampTab, minNum=grnSampSize)
  cat("Number of samples per CT: ",mean(table(stGRN[,dLevel])),"\n")
  expGRN<-expDat[,rownames(stGRN)]
  if(normDat){
    cat("Normalizing expression data...\n")
    expGRN<-Norm_quantNorm(expGRN)
  }

  if(is.na(corrs)){
    cat("Calculating correlation...\n");
    corrs<-grn_corr_round(expGRN)
  }
  if(is.na(zscores)){
    cat("Calculating context dependent zscores...\n");
    zscores<-grn_zscores(corrs, tfs);
  }

  ### cat("n target genes ... ", length(targetGenes),"\n");
  ### cat("n tfs ... ",length(tfs), "\n");

  grnall<-cn_getRawGRN(zscores, corrs, targetGenes, zThresh=zThresh)#, snName=snName);

cat("specGenes all\n")
  specGenes<-cn_specGenesAll(expGRN, stGRN, holm=holmSpec, cval=cval, cvalGK=cvalGK, dLevel=dLevel, dLevelGK=dLevelGK, prune=prune)
cat(" done specGenes all\n")
cat(" specGenes\n")
  ctGRNs<-cn_specGRNs(grnall, specGenes);
  cat("done specGenes\n")
  ###list(overallGRN=grnall, specTFs=specTFs,ctGRNs=ctGRNs);  
  list(overallGRN=grnall, specGenes=specGenes,ctGRNs=ctGRNs, grnSamples=rownames(stGRN));  
}


#' get raw GRN from zscores, and corr
#'
#' get raw GRN from zscores, and corr
#' @param zscores zscores matrix
#' @param corrs correlation matrix
#' @param targetGenes target genes
#' @param zThresh zscore threshold
#'
#' @return list of grnTable and corresponding graph
cn_getRawGRN<-function# get raw GRN, communities from zscores, and corr
(zscores, ### zscores matrix
 corrs, #### correlation matrix
 targetGenes,#### target genes
 zThresh=4### zscore threshold
### snName="raw" ### subnetwork prefix
 ){
 
  # make a grn table
  cat("Making GRN table...\n")
  grn<-cn_extractRegsDF(zscores, corrs, targetGenes, zThresh);
  cat("Done making GRN table...\n")
  colnames(grn)[1:2]<-c("TG", "TF");
  
  # make an iGraph object and find communities
  cat("Make iGraph...\n")
  igTmp<-ig_tabToIgraph(grn, directed=FALSE, weights=TRUE);
   cat("done with iGraph...\n")
  list(grnTable=grn, graph=igTmp);
}

#' finds general and context dependent specifc genes
#'
#' finds general and context dependent specifc genes
#' @param expDat expression matrix
#' @param sampTab sample table
#' @param holm pvalue threshold for template matching
#' @param cval template matching threshold for overall CT specific expression
#' @param cvalGK template matching threshold for developmentally shared CT specific expression
#' @param dLevel "description1",
#' @param dLevelGK "description2"
#'
#' @return list of $matcher${cell_type}->{germ_layer}$context$general${cell_type}->gene vector etc
cn_specGenesAll<-function
(expDat, 
 sampTab,
 holm=1e-10,
 cval=0.5,
 cvalGK=0.75,
 dLevel="description1",
 dLevelGK="description2",
 prune=FALSE){
  matcher<-list();
  general<-cn_findSpecGenes(expDat, sampTab, holm=holm, cval=cval, dLevel=dLevel,prune=prune);
  ctXs<-list()# one per germlayer
  if(!is.null(dLevelGK)){
    
    germLayers<-unique(as.vector(sampTab[,dLevelGK]));
    for(germlayer in germLayers){
      stTmp<-sampTab[sampTab[,dLevelGK]==germlayer,];
      expTmp<-expDat[,rownames(stTmp)];
      xxx<-cn_findSpecGenes(expTmp, stTmp, holm=holm, cval=cvalGK,dLevel=dLevel, prune=prune);
      cts<-names(xxx);
      for(ct in cts){
        matcher[[ct]]<-germlayer;
        # remove general ct-specific genes from this set
        a<-general[[ct]];
        b<-xxx[[ct]];
        ba<-setdiff(b, a);
        both<-union(a,b);
        xxx[[ct]]<-ba;
      }
      ctXs[[germlayer]]<-xxx;
    }
  }
  ctXs[['general']]<-general;
  list(context=ctXs, matcher=matcher);
}

if(FALSE){ #10-06-16 had implemented in order to find CT GRNs by first finding CT-enriched TFs, then just finding all of these TFs
cn_specGRNs<-function### extract sub-networks made up of CT genes; don't bother finding communities
(rawGRNs, ### result of running cn_getRawGRN
 specTFs ### result of running cn_specGenesAll
 ){

  ## ct GRN is graph seeded by spec TFS, and all targets thereof (spec or not)
  # should return a list of gene lists and igraphs
  geneLists<-list();
  graphLists<-list();

  groupNames<-names(specTFs[['context']][['general']]);  
  big_graph<-rawGRNs[['graph']];

  matcher<-specTFs$matcher;
  
  allgenes<-V(big_graph)$name;
  tfTargets<-list();

  for(ct in groupNames){
    gll<-matcher[[ct]];
    myTFs<-union(specTFs[['context']][['general']][[ct]], specTFs[['context']][[gll]][[ct]]);
    cat(ct," ",gll, " n tfs: ",length(myTFs), "\n");
   # tfTargetsTmp<-cn_get_targets_of(aGraph, myTFs)
   tfTargetsTmp<-cn_get_targets_of(big_graph, myTFs)
    mygenes<-union(myTFs, unlist(tfTargetsTmp))
    cat("mygenes: ", length(mygenes),"\n");
    geneLists[[ct]]<-intersect(allgenes, mygenes);
    graphLists[[ct]]<-induced.subgraph(big_graph, geneLists[[ct]]);
    tfTargets[[ct]]<-tfTargetsTmp
  }

#  tfTargets<-cn_MakeTLs(graphLists);
   
  list(geneLists=geneLists, graphLists=graphLists, tfTargets=tfTargets);

}
}

#' extract sub-networks made up of CT genes; 
#'
#' extract sub-networks made up of CT genes; 
#' @param rawGRNsresult of running cn_getRawGRN
#' @param specGenes result of running cn_specGenesAll
#'
#' @return list(geneLists=geneLists, graphLists=graphLists, tfTargets=tfTargets)
cn_specGRNs<-function### don't bother finding communities
(rawGRNs, ### result of running cn_getRawGRN
 specGenes ### result of running cn_specGenesAll
 ){

  # should return a list of gene lists and igraphs
  geneLists<-list();
  graphLists<-list();

  groupNames<-names(specGenes[['context']][['general']]);  

  big_graph<-rawGRNs[['graph']];

  matcher<-specGenes$matcher;
  
  allgenes<-V(big_graph)$name;

  for(ct in groupNames){
    cat(ct,"\n")
    if(!is.null(names(matcher))){
      gll<-matcher[[ct]];
      cat(ct," ",gll,"\n");
      mygenes<-union(specGenes[['context']][['general']][[ct]], specGenes[['context']][[gll]][[ct]]);
    }
    else{
      mygenes<-specGenes[['context']][['general']][[ct]]
    }
    
    geneLists[[ct]]<-intersect(allgenes, mygenes);
    graphLists[[ct]]<-induced.subgraph(big_graph, geneLists[[ct]]);
  }

  tfTargets<-cn_MakeTLs(graphLists);
   
  list(geneLists=geneLists, graphLists=graphLists, tfTargets=tfTargets);
}

#' get targets of tFs
#'
#' get targets of tFs
#' @param graphList a list of networks represented as iGraphs
#'
#' @return list of tf=>targets
cn_MakeTLs<-function(graphList){
  tfTargs<-list();
  nnames<-names(graphList);
  for(nname in nnames){
    tfTargs[[nname]]<-cn_get_targets(graphList[[nname]]);
  }
  tfTargs;
}

#' get targets of a tf
#'
#' get targets of a tf
#' @param aGraph an iGraph
#'
#' @return target list
cn_get_targets<-function(aGraph){
  targList<-list();
  regs<-V(aGraph)$label[V(aGraph)$type=='Regulator'];
  if(length(regs)>0){ 
    for(reg in regs){
       targList[[reg]]<-unique(sort(V(aGraph)$label[neighbors(aGraph, reg)]));
    }
  }
  targList;
}

#' get targets of tfs
#'
#' get targets of tfs
#' @param aGraph iGraph network
#' @param tfs vector of TF names
#'
#' @return list of tf-> targets
cn_get_targets_of<-function(aGraph, tfs){
  targList<-list();
  regs<-V(aGraph)$label[V(aGraph)$type=='Regulator'];
  regs<-intersect(regs, tfs);
  if(length(regs)>0){ 
    for(reg in regs){
       targList[[reg]]<-unique(sort(V(aGraph)$label[neighbors(aGraph, reg)]));
    }
  }
  targList;
}

#' find genes that are preferentially expressed in specified samples
#'
#' find genes that are preferentially expressed in specified samples
#' @param expDat expression matrix
#' @param sampTab sample table
#' @param holm sig threshold
#' @param cval R thresh
#' @param dLevel annotation level to group on
#' @param prune boolean limit to genes exclusively detected as CT in one CT
#'
#' @return a list of something
#'
#' @export
#'
cn_findSpecGenes<-function# 
(expDat, ### expression matrix
 sampTab, ### sample tableß
 holm=1e-50, ### sig threshold
 cval=0.5, ### R thresh
 dLevel="description1", #### annotation level to group on
 prune=FALSE ### limit to genes exclusively detected as CT in one CT
 ){
  
  myPatternG<-cn_sampR_to_pattern(as.vector(sampTab[,dLevel]));
  specificSets<-apply(myPatternG, 1, cn_testPattern, expDat=expDat);

  # adaptively extract the best genes per lineage
  cvalT<-vector();
  ctGenes<-list();
  ctNames<-unique(as.vector(sampTab[,dLevel]));
  for(ctName in ctNames){
    x<-specificSets[[ctName]];
    tmp<-rownames(x[which(x$cval>cval),]);
    tmp2<-rownames(x[which(x$holm<holm),]);
    tmp<-intersect(tmp, tmp2)
    ctGenes[[ctName]]<-tmp;
###    cvalT<-append(cvalT, cval);
  }
  
  if(prune){
  # now limit to genes exclusive to each list
   specGenes<-list();
   for(ctName in ctNames){
     others<-setdiff(ctNames, ctName);
     x<-setdiff( ctGenes[[ctName]], unlist(ctGenes[others]));
     specGenes[[ctName]]<-x;
   }
   ans<-specGenes;
  }
  else{
   ans<-ctGenes;
  }
  ans;
}

#' make induced subgraphs from gene lists and iGraph object
#'
#' make induced subgraphs from gene lists and iGraph object
#' @param bigIG super-graph
#' @param geneLists gene lists
#'
#' @return list of graphs
cn_makeSGs<-function# 
(bigIG,### super-graph
 geneLists ### gene lists
){
  
  allGenes<-V(bigIG)$name;
  subGraphs<-list();
  nnames<-names(geneLists); 
  for(nname in nnames){
    genes<-intersect(allGenes, geneLists[[nname]]);    
    subGraphs[[nname]] <-  induced.subgraph(bigIG, genes);
  }
  subGraphs;
}

#' convert a table to an igraph
#'
#' convert a table to an igraph. This adds an nEnts vertex attribute to count the number of entities in the sub-net
#' 
#' @param grnTab table of TF, TF, maybe zscores, maybe correlations
#' @param simplify false
#' @param directed FALSE,
#' @param weights TRUE
#'
#' @return iGraph object
ig_tabToIgraph<-function#
(grnTab,
 simplify=FALSE,
 directed=FALSE,
 weights=TRUE
){
  
  tmpAns<-as.matrix(grnTab[,c("TF", "TG")]);
  regs<-as.vector(unique(grnTab[,"TF"]));
  ###cat("Length TFs:", length(regs), "\n");
  targs<-setdiff( as.vector(grnTab[,"TG"]), regs);
  
###  cat("Length TGs:", length(targs), "\n");
  myRegs<-rep("Regulator", length=length(regs));
  myTargs<-rep("Target", length=length(targs));
  
  types<-c(myRegs, myTargs);
  verticies<-data.frame(name=c(regs,targs), label=c(regs,targs), type=types);
  
  ### iG<-graph.data.frame(tmpAns,directed=directed,v=verticies);
  iG<-igraph::graph_from_data_frame(tmpAns,directed=directed,v=verticies);
  
  if(weights){
    #E(iG)$weight<-grnTab$weight;    
    E(iG)$weight<-grnTab$zscore;    
  }
  
  if(simplify){
    iG<-simplify(iG);
  }
  V(iG)$nEnts<-1;
  iG;
}

#' return a pattern for use in cn_testPattern (template matching)
#'
#' return a pattern for use in cn_testPattern (template matching)
#' @param sampR vector
#'
#' @return ans
cn_sampR_to_pattern<-function# 
(sampR){
  d_ids<-unique(as.vector(sampR));
  nnnc<-length(sampR);
  ans<-matrix(nrow=length(d_ids), ncol=nnnc);
  for(i in seq(length(d_ids))){
    x<-rep(0,nnnc);
    x[which(sampR==d_ids[i])]<-1;
    ans[i,]<-x;
  }
  colnames(ans)<-as.vector(sampR);
  rownames(ans)<-d_ids;
  ans;
}

#' template matching
#'
#' test correlation between idealized expression pattern and target gene
#' @param pattern vector of pattern
#' @param expDat expression matrix
#'
#' @return data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
#'
cn_testPattern<-function(pattern, expDat){
  pval<-vector();
  cval<-vector();
  geneids<-rownames(expDat);
  llfit<-ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE);
  xxx<-matrix( unlist(llfit$coef), ncol=8,byrow=TRUE);
  ccorr<-xxx[,6];
  cval<- sqrt(as.numeric(llfit$summary[,2])) * sign(ccorr);
  pval<-as.numeric(xxx[,8]);
  
  #qval<-qvalue(pval)$qval;
  holm<-p.adjust(pval, method='holm');
  #data.frame(row.names=geneids, pval=pval, cval=cval, qval=qval, holm=holm);
  data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
}

#' extracts the TRs, zscores, and corr values passing thresh
#'
#' extracts the TRs, zscores, and corr values passing thresh
#' @param zscores, # zscore matrix, non-TFs already removed from columns
#' @param corrMatrix, # correlation matrix
#' @param genes, # vector of target genes
#' @param threshold # zscore threshold
#'
#' @return data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations);
cn_extractRegsDF<-function# 
(zscores,
 corrMatrix,
 genes,
 threshold
){
    
  targets<-vector();
  regulators=vector();
  zscoresX<-vector();
  correlations<-vector();
  
  targets<-rep('', 1e6);
  regulators<-rep('', 1e6);
  zscoresX<-rep(0, 1e6);
  correlations<-rep(0, 1e6);
  
  str<-1;
  stp<-1;
  for(target in genes){
    x<-zscores[target,];
    regs<-names(which(x>threshold));
    if(length(regs)>0){
      zzs<-x[regs];
      corrs<-corrMatrix[target,regs];
      ncount<-length(regs);
      stp<-str+ncount-1;
      targets[str:stp]<-rep(target, ncount);
      #    targets<-append(targets,rep(target, ncount));
      regulators[str:stp]<-regs;
      #regulators<-append(regulators, regs);
      #    zscoresX<-append(zscoresX, zzs);
      zscoresX[str:stp]<-zzs;
      correlations[str:stp]<-corrs;
      str<-stp+1;
    }
    #    correlations<-append(correlations, corrs);
  }
  targets<-targets[1:stp];
  regulators<-regulators[1:stp];
  zscoresX<-zscoresX[1:stp];
  correlations<-correlations[1:stp];
  
  
  data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations);
}

#' sample equivalent numbers of profiles per cell type
#'
#' sample equivalent numbers of profiles per cell type
#' @param sampTab sample table
#' @param minNum min number of samples to get per CT
#' @param dLevel grouping
#'
#' @return subset of sample table
sample_profiles_grn<-function# 
(sampTab,### sample table
  minNum=NULL, # min number of samples to get per CT
  dLevel="description1" ### grouping
  ){

  nperCT<-table(sampTab[,dLevel]);   
  if(is.null(minNum)){
    minNum<-min(nperCT);
  }

  nsamps<-vector();
  ctts<-names(nperCT);
  for(ctt in ctts){
    stTmp<-sampTab[sampTab[,dLevel]==ctt,];
    cat(ctt,":",nrow(stTmp),"\n");
    ids<-sample( rownames(stTmp),minNum);
    nsamps<-append(nsamps, ids);
  }
  sampTab[nsamps,];
}

#' gene-gene correlations, and round
#'
#' gene-gene correlations, and round
#' @param expDat expression matrix
#'
#' @return correlation matrix
grn_corr_round<-function # 
(expDat
  ){
  corrX<-cor(t(expDat));
  round(corrX, 3);
}

#' compute CLR-like zscores
#'
#' compute CLR-like zscores
#' @param corrVals correlation matrix
#' @param tfs vector of  transcriptional regualtor names
#'
#' @return zscore matrix
grn_zscores<-function 
(corrVals,
 tfs
  ){
  zscs<-mat_zscores(corrVals);
  gc();
  zscs[,tfs];
}

#' compute context dependent zscores
#'
#' slightly modidied from JJ Faith et al 2007
#' @param corrMat correlation matrix
#'
#' @return matrix of clr zscores
#'
mat_zscores<-function# computes sqrt(zscore_row + zscore_col) .. 
(corrMat ### correlation matrix
){
  corrMat<-abs(corrMat);
  zscs_2<-round(scale(corrMat), 3);
  rm(corrMat);
  gc()
  zscs_2 + t(zscs_2);
}

#' make subnets from a GRN
#'
#' make subnets from a GRN
#' @param aGraph igraph
#' @param geneLists named list of genes
#'
#' @return list(subnets=list(graphs=graphs, geneLists=geneLists), general=list(graphs=graphs, geneLists=geneLists));
#'
cn_extractSubNets<-function# 
(aGraph, # igraph
 geneLists # named list of genes
){
    graphs<-list();
    nnames<-names(geneLists);
    genesInGraph<-V(aGraph)$name;
    for(nname in nnames){
      graphs[[nname]]<-induced.subgraph(aGraph, intersect(geneLists[[nname]], genesInGraph));
    }
    list(subnets=list(graphs=graphs, geneLists=geneLists), general=list(graphs=graphs, geneLists=geneLists));
}


#' find transcript factors
#'
#' find transcript factors
#' @param species defaul is 'Hs', can also be 'Mm;
#'
#' @return vector fo TF names
#' @export
#' @importFrom AnnotationDbi as.list
#'
find_tfs<-function# 
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


#' computes the raw score for a gene as xmax-abs(zscore).
#'
#' better values are higher
#' @param vect a vector of gene expression values for multiple samples
#' @param mmean mean value in training data
#' @param ssd standard deviation in training data
#'
#' @return transformed (but not normalized) GRN score
#'
cn_rawScore<-function
(vect,
 mmean,
 ssd,
 xmax=1e3
){
  zcs<-zscore(vect, mmean, ssd);
  ### xmax<-1000; # arbitrary, and corrected for later, but want some high enough that it should not be exceeded 
  xmax-abs(zcs);
}


#' min diff 
#'
#' computes mean gene A in CT 1 - mean gene A in CT 2, where CT2 has the non CT1 max value. does this for genes
#' @param tVals tVals 
#' @param genes vector of gene names
#' @param ct ct to compare to
#' @return vector of differences
minDif<-function
(tVals,
  genes,
  ct){
  octs<-setdiff(names(tVals), ct)
  qq<-lapply(tVals[octs], "[[", "mean")
  ##tVals[[ct]][["mean"]][[gene]]#-max(unlist(lapply(qq, "[[", gene)))
  tmpMat<-matrix(unlist(lapply(qq, "[", genes)), nrow=length(genes))
  rownames(tmpMat)<-genes
  maxes<-apply(tmpMat, 1, max)
  unlist(tVals[[ct]][["mean"]][genes])-maxes
}


#' GRN status
#'
#' Calculates the status of all GRNs in query samples as compared to training data for
#' @param expDat query expression matrix
#' @param subList of ct => genes
#' @param tVals tvals
#' @param classList classList
#' @param minVals minVals
#' @param classWeight class weight
#' @param exprWeight  expression weight
#' @return grn scores (not normalized)
cn_netScores<-function
### return the GRN establishment score for a given expression matrix
(expDat, 
 genes, 
 tVals, 
 ctt, 
 classList=NULL, 
 classWeight=FALSE, 
 exprWeight=TRUE,
 xmax=1e3
){
  cat(ctt,"\n")
  aMat<-matrix(0, nrow=length(genes), ncol=ncol(expDat));
  rownames(aMat)<-genes;
  
  weights<-rep(1, length(genes));
  names(weights)<-genes;
  
  #otherCTs<-setdiff(names(tVals), ct)
  
  cat(dim(aMat),"\n")
  if(exprWeight){
    meanVect<-unlist(tVals[[ctt]][['mean']][genes]);
    weights<-(2**meanVect)/sum(2**meanVect);
    if(FALSE){
      bDifs<-minDif(tVals, genes, ctt)
    # set to zero neg vals
      bDifs[bDifs<0]<-0
      weights<-bDifs/sum(bDifs)
    }
  }
    
  if(classWeight){
    classImp<-classList[[ctt]]$importance[genes,1];
    ### 04-19-17
    ###classImp<-classImp/sum(classImp)
    weights<-weights*classImp;
  }
  
  for(gene in genes){
   ### cat("***",gene,"\n")
    ###zzs<-as.matrix(cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]])[1,])


    zzs<-cn_rawScore(expDat[gene,], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]], xmax=xmax)
    aMat[gene,]<-zzs;

  }
  xscores<-apply(aMat, 2, weighted.mean, w=weights);
  xscores;
}

#' GRN status
#'
#' Calculates the GRN status in query samples as compared to training data
#' @param expDat query expression matrix
#' @param subList of ct => genes
#' @param tVals list of ctt->list(means->genes, sds->genes)
#' @param classList class list
#' @param minVals  min vals
#' @param classWeight classweght
#' @param exprWeight  expression weight
#' @return GRN scores 
#'
cn_score<-function
(expDat,
 subList, 
 tVals,
 classList=NULL,
 minVals=NULL, 
 classWeight=FALSE,
 exprWeight=TRUE,
 xmax=1e3
){
  #nSubnets<-sum(sapply(subList, length));
  nSubnets<-length(subList);
  ans<-matrix(0, nrow=nSubnets, ncol=ncol(expDat));
  ctts<-names(subList);
  rnames<-vector();
  rIndex<-1;
  for(ctt in ctts){
     cat(ctt,"\n");
    genes<-subList[[ctt]];
    # 06-06-16 -- added to allow for use of GRNs defined elsewhere
    genes<-intersect(genes, rownames(expDat));
    #    snNames<-names(subnets);
    #    rnames<-append(rnames, snNames);    
    #    for(sName in snNames){
    ans[rIndex,]<-cn_netScores(expDat, genes, tVals=tVals, ctt=ctt,classList, classWeight=classWeight,exprWeight=exprWeight, xmax=xmax);
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
}

#' Normalize grn status as compared to training data
#'
#' Divide the query scores by the mean values in the training data.
#' @param ctrlScores a list of subnet->mean value, all subnets
#' @param queryScores a matrix, rownames = subnet names, all subnets
#' @param subNets a vector of subnets names to use
#'
#' @return normalized grn status matrix
#'
cn_normalizeScores<-function
(ctrlScores, 
 queryScores, 
 subNets 
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
}

#
# functions to enable GRN status metric
#

#' Figure out normalization factors for GRNs, and norm training data
#'
#' Exactly that.
#' @param expTrain expression matrix
#' @param stTrain  sample table
#' @param subNets named list of genes, one list per CTT, tct=>gene vector
#' @param classList list of classifiers
#' @param dLevel column name to group on
#' @param tVals seful when debugging
#' @param classWeight weight GRN status by importance of gene to classifier
#' @param exprWeight weight GRN status by expression level of gene?
#' @param sidCol sample id colname
#' 
#' @return list of trainingScores, normVals, raw_scores, minVals, tVals=tVals
#' @export
cn_trainNorm<-function # 
(expTrain,
 stTrain,
 subNets,
 classList = NULL, 
 dLevel = "description1",
 tVals=NULL,
 classWeight=FALSE,
 exprWeight=TRUE,
 sidCol='sample_id',
 xmax=1e3,
 predSD=FALSE
){

  if(is.null(tVals)){
    tVals<-cn_make_tVals(expTrain, stTrain, dLevel, predictSD=predSD)
  }

  ctts<-as.vector(unique(stTrain[,dLevel]));
  scoreList<-list();
  normList<-list(); # a list of ctt->subnet->mean value
  minVect<-vector(); # a list of ctt->subnet->min value, used to shift raw grn est scores
  
  cat("calculating GRN scores on training data ...\n");
  tmpScores<-cn_score(expTrain, subNets, tVals, classList, minVals=NULL, classWeight=classWeight, exprWeight=exprWeight, xmax=xmax)


  minVect<-apply(tmpScores, 1, min);
  names(minVect)<-rownames(tmpScores);
  
  # shift the raw scores so that min=0;
  tmpScores<-tmpScores - minVect;
  cat("norm factors\n");
  for(ctt in ctts){
    # determine nomalization factors
    ##snets<-names(subNets[[ctt]]);
    snets<-ctt;

    scoreDF<-cn_extract_SN_DF(tmpScores, stTrain, dLevel, snets, sidCol=sidCol);
    scoreDF<-cn_reduceMatLarge(scoreDF, "score", "description", "subNet");
    xdf<-scoreDF[which(scoreDF$grp_name==ctt),];
    tmpSNS<-as.list(xdf$mean);
    names(tmpSNS)<-xdf$subNet;
    normList[names(tmpSNS)]<-tmpSNS;      
  }
  
  # normalize training scores
  nScores<-cn_normalizeScores(normList, tmpScores, rownames(tmpScores));

  scoreDF<-cn_extract_SN_DF(nScores, stTrain, dLevel, sidCol=sidCol);

  scoreDF<-cn_reduceMatLarge(scoreDF, "score", "description", "subNet");

  list(trainingScores=scoreDF,
       normVals=normList,
       raw_scores=tmpScores,
       minVals=minVect,
       tVals=tVals);
}

#' Estimate gene expression dist in CTs
#'
#' Calculate mean and SD 
#' @param expDat training data 
#' @param sampTab, ### training sample table 
#' @param dLevel="description1", ### column to define CTs
#' @param predictSD=FALSE ### whether to predict SD based on expression level
#'
#' @return tVals list of ct->mean->named vector of average gene expression, ->sd->named vector of gene standard deviation
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

#' cn_make_tVals_predict
#'
#' predicts SD based on mean expression
#' @param expDat training data 
#' @param sampTab training sample table 
#' @param dLevel="description1" column to define CT
#'
#' @return tVals list of ct->mean->named vector of average gene expression, ->sd->named vector of gene standard deviation
cn_make_tVals_predict<-function ### 
(expDat,
 sampTab,
 dLevel="description1" 
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

#' convert a tf nis list to a DF
#'
#' convert a tf nis list to a DF
#' @param tfScores tf scores
cn_makeTFtable<-function
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


# end grn_status.R
########################################################################################


########################################################################################
# start cellnetr_cnres_ops.R

#' Network Influence Score for all GRNs
#'
#' Runs cn_nis on all GRNs
#' @param cnRes object result of running cn_apply
#' @param cnProc object result of running cn_make_processor
#' @param ctt string indicating the CT to compare against
#' @param relaWeight whether to weight by overall expression such that TFs with higher expression in ctt are more important (1=do the weighting)
#'
#' @return list of numeric matrix of TF scores
#'
#' @export
cn_nis_all<-function
(cnRes,
 cnProc,
 ctt,
 relaWeight=1
 ){
  snNames<-names(cnProc$ctGRNs$ctGRNs$graphLists);
  ans<-list()
  for(snName in snNames){
    cat("scoring ", snName,"\n")
    x<-cn_nis(cnRes, cnProc, snName, ctt,relaWeight);
    ans[[snName]]<-x;
  }
  ans;
}

#' network influence score
#'
#' Computes network influence score (NIS). See paper for details.
#' @param cnRes object result of running cn_apply
#' @param cnProc object result of running cn_make_processor
#' @param subnet name of subnet to evaluage
#' @param ctt string indicating the CT to compare against
#'
#' @return numeric matrix where rows are TFs in CT GRN and columns are query samples
#'
#' @export
cn_nis<-function
(cnRes,
 cnProc,
 subnet,
 ctt,
 relaWeight=1
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

trimmean<-function
(vect){
  mean(quantile(vect, c(.25, .5, .5,.75)))
}

#' transcription factor score
#'
#' Computes TF score (TFS). See forthcoming paper for details.
#' @param cnRes object result of running cn_apply
#' @param cnProc object result of running cn_make_processor
#' @param tfname transcription factor
#' @param ctt string indicating the CT to compare against
#'
#' @return numeric matrix where rows are TFs and columns are query samples, values are zscores
#'
#' @export
cn_tfScore<-function
(cnProc,
 expDat,
 ctt,
 tfnames=NULL
){
  
  if(is.null(tfnames)){
    tfnames<-unique(unlist(lapply(cnProc$ctGRNs$ctGRNs$tfTargets, names)))
  }
  tVals<-cnProc[['tVals']][[ctt]]
  # compute a matrix of zscores


  ans<-matrix(0, nrow=length(tfnames), ncol=ncol(expDat))

  bGraph<-cnProc$ctGRNs$overallGRN$graph
  for(i in seq(length(tfnames))){

    tfname<-tfnames[i]
    cat(tfname,"\n")
    tgenes<-unique(V(bGraph)$label[neighbors(bGraph, tfname)])
    
    if(FALSE){
      zzzMat<-matrix(0, nrow=length(tgenes), ncol=ncol(expDat))  
      for(j in seq(length(tgenes))){
       gene<-tgenes[j]
       zzzMat[j,]<-zscore(expDat[gene,], tVals[['mean']][[gene]], tVals[['sd']][[gene]])
      }
      zzzMat<-cn_correctZmat(zzzMat)
    }
    zzzMat<-cn_zmat(expDat[tgenes,], tVals)
    ans[i,]<-apply(abs(zzzMat), 2, median)
  }
  
  tfzs<-cn_zmat(expDat[tfnames,], tVals)
  # scale this by max(expr_ctt, expr_query)
  geneScale<-cn_exprScale(expDat[tfnames,], tVals[['mean']])
  ans<-tfzs * ans * geneScale

  rownames(ans)<-tfnames
  colnames(ans)<-colnames(expDat)
  ans
}

cn_exprScale<-function
(expDat,
 tVals)
 {
  genes<-rownames(expDat)
  aMat<-matrix(0, nrow=length(genes), ncol=ncol(expDat))
  for(i in seq(ncol(expDat))){  
    aMat[,i]<-pmax(unlist(tVals[genes]), expDat[genes,i])
  }
  rownames(aMat)<-genes
  colnames(aMat)<-colnames(expDat)
  aMat
}


cn_zmat<-function
(expDat,
 tVals){
  tgenes<-rownames(expDat)
  zzzMat<-matrix(0, nrow=length(tgenes), ncol=ncol(expDat))  
  for(j in seq(length(tgenes))){
    gene<-tgenes[j]
    
    zzzMat[j,]<-zscore(expDat[gene,], tVals[['mean']][[gene]], tVals[['sd']][[gene]])
  }
  zzzMat<-cn_correctZmat(zzzMat)
  colnames(zzzMat)<-colnames(expDat)
  rownames(zzzMat)<-tgenes
  zzzMat
}


cn_nis_bd<-function
(cnRes,
 cnProc,
 subnet,# network
 ctt, # target ct
 relaWeight=1
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
  
 ### ans<-matrix(0, nrow=length(tfs), ncol=nrow(stQuery));
  ans<-list()
 ### rownames(ans)<-tfs;
 ### colnames(ans)<-sids;
  
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
  zzzMat<-cn_correctZmat(zzzMat)
   rownames(zzzMat)<-netGenes;
  colnames(zzzMat)<-rownames(stQuery);
 
  

###    meanVect<-unlist(tVals[[subnet]][['mean']][netGenes]);
    meanVect<-unlist(tVals[[ctt]][['mean']][netGenes]);
    weights<-(2**meanVect)/sum(2**meanVect);

  for(sid in sids){
    #cat("tf scoring ", sid,"\n");
    xvals<-as.vector(expDat[,sid]);
    names(xvals)<-rownames(expDat);
    
  
    # assign weights
    
    ### # meanVect<-unlist(tVals[[ctt]][['mean']][netGenes]);
    ####meanVect<-unlist(tVals[[subnet]][['mean']][netGenes]);
    ####weights<-(2**meanVect)/sum(2**meanVect);
    
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
      ###if(targetScore[i]=='Inf'){
      ###  ans<-list(tmp=tmp, sid=sid, weights=weights[targs])
      ###  break
      ###}
      
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
    ###ans[as.vector(xxx$tf),sid]<-as.vector(xxx$totalScore);
   ans[[sid]]<-xxx
  }
 ### ans;
 ans

}

###########################################################################
#
# UNUSED
#
###########################################################################

if(FALSE){
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
}


if(FALSE){
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
}

if(FALSE){
  .get_cttName<-function
### splits a sub-network name into a ctt
(snName
 ### subnet name
 ){
  x<-strsplit(snName, "_")[[1]][1];
  x;
}
}

if(FALSE){


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

}


if(FALSE){

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
}
if(FALSE){
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
}

if(FALSE){
  cn_findGene<-function # find what subents a gene is in
(sns,
 gene){
  ans<-vector();
  snNames<-names(sns);
  for(snName in snNames){
    if( any(sns[[snName]]==gene) ){
      ans<-append(ans, snName);
    }
  }
  ans;
}
}

if(FALSE){
cn_getGeneWeights<-function# extract the gene weighting for each gene in a c/t GRN
(cnProc,
 ctt){

  tVals<-cnProc[['tVals']];
  classList<-cnProc[['classList']];
  genes<-cnProc[['grnList']][[ctt]];
  weights<-rep(1, length(genes));
  names(weights)<-genes;
  if(cnProc[['exprWeight']]){
    cat("expr weights\t");
    meanVect<-unlist(tVals[[ctt]][['mean']][genes]);
    weights<-(2**meanVect)/sum(2**meanVect);
  }
  if(cnProc[['classWeight']]){
    cat("class weights\n")
    classImp<-classList[[ctt]]$importance[genes,1];
    weights<-weights*classImp;
  }

  names(weights)<-genes;
  weights;
}

}

if(FALSE){
  cn_commToNames<-function # return a named list, wrapper to celnet_commToNames
(commObj,
 prefix
){
  ans<-list();
  comms<-communities(commObj);
  for(i in seq(length(comms))){
    nname<-paste(prefix,"_sn_",i,sep='');
    ans[[nname]]<-commObj$names[comms[[i]]];
  }
  ans;
}
}

# end cellnetr_cnres_ops.R
##############################################################################









