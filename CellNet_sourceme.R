# CellNet 
# (C) Patrick Cahan 2012-2016

library(igraph)
###library(gplots)
library(ggplot2)
library(randomForest)
library(preprocessCore)
library(parallel)
library(GO.db)


utils_sourceRs<-function# source all .R files in given directory
(dirname ### directory containing R files for sourcing
  ){
  cmd<-paste("ls ",dirname,"*.R", sep='');
  rfiles<-system(cmd, intern=TRUE);
  for(rfile in rfiles){
    source(rfile);
  }
}




