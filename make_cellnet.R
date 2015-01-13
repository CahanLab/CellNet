# CellNet
# (C) Patrick Cahan 2012-2014

# make a CellNet object

make_classifiers<-function# wrappper to cn_makeRFs
(ctGRNs, ### result of running cn_grnDoRock
 expDat, ### training data,
 sampTab, ### sample table
 dLevel ### description level to indicate cell types
 ){
	geneLists<-ctGRNs$ctGRNs$general$geneLists;
 	cn_makeRFs(expDat,sampTab,geneLists, dLevel=dLevel);
}

