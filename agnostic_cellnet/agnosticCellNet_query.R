# ---
# title: "agnosticCellNet Query Instructions"
# author: "Emily Lo"
# Last updated: "10/15/2020"
# ---
# Based on cancerCellNet README

### Prerequisite Installations ----------------------------------------------
# Only need to install each package on your computer ONCE, then comment out this section
#install packages from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# 
# BiocManager::install("AnnotationDbi")
# BiocManager::install("GO.db")
# BiocManager::install("org.Hs.eg.db")
# 
# install.packages("pheatmap")
# install.packages("RColorBrewer")
# install.packages("randomForest")
# install.packages("ggplot2")
# install.packages("igraph")
# install.packages("stringr")
# library(devtools)
# install_github("pcahan1/cancerCellNet@v0.1.1", ref="master")
library(cancerCellNet)
library(CellNet)
library(ggplot2)


library(RColorBrewer)
library(pheatmap)
source("plotting.R")

### HUMAN Query sample classification ---------------------------------------

# Re-train the random forest classifier using only intersecting genes (`iGenes`)

setwd("AgnosticCellNet/Scripts") # CHANGE FOR YOUR OWN DIRECTORY STRUCTURE
expTrain <- utils_loadObject("../Hs_All_Train/Hs_expTrain_Jun-20-2017.rda")
stTrain <- utils_loadObject("../Hs_All_Train/Hs_stTrain_Jun-20-2017.rda")
broadClass <- utils_loadObject("../Hs_All_Train/Hs_TP_broadClassifier100_Apr_22_2020.rda")

setwd("~/Documents/JohnsHopkins/CahanLab/AgnosticCellNet/TP_train/Hs/neuron")
study_name <- "SRP032798" # REPLACE WITH YOUR STUDY NAME

# For author-preprocessed:
queryExpDat <- utils_loadObject(paste0(study_name, "/", study_name,"_expMat_counts.rda"))
querySampTab <- utils_loadObject(paste0(study_name, "/", study_name,"_sampTab.rda"))

# OR

# For CloudSeq-preprocessed:
queryExpDat <- utils_loadObject(paste0(study_name, "/expGeneList_", study_name,".R"))
queryExpDat <- queryExpDat$counts
querySampTab <- utils_loadObject(paste0(study_name, "/sampTab_", study_name,"_QC.R"))

# Normalize expression matrix
queryExpDat <- log(1+queryExpDat)
#queryExpDat <- apply(queryExpDat, 2, downSampleW, 1e5) # No longer downsampling when using TP classifier

# Find intersecting genes between query sample and training samples
iGenes <- Reduce(intersect, list(rownames(queryExpDat), rownames(expTrain)))
# OR when further subsetting:
#iGenes <- Reduce(intersect, list(rownames(queryExpDat), iGenes))

save(iGenes, file = paste0("neuron_14724-iGenes_", study_name, ".rda"))

# Split expTrain into training and validation
set.seed(39) # Set seed for reproducibility
stList <- splitCommon_proportion(sampTab = stTrain, proportion = 0.66, dLevel = "description1")
stTrainSubset <- stList$trainingSet
expTrainSubset <- expTrain[,stTrainSubset$sra_id]

# Train the random forest classifier
# Will take 3-5 minutes
system.time(broad_return <- broadClass_train(stTrain = stTrainSubset, 
                                expTrain = expTrainSubset[iGenes, ], 
                                colName_cat = "description1", 
                                colName_samp = "sra_id", 
                                nRand = 70,
                                nTopGenes = 100, 
                                nTopGenePairs = 100, 
                                nTrees = 2000, 
                                stratify=TRUE, 
                                sampsize=25, 
                                quickPairs=TRUE))

#save(broad_return, file=paste0("broadClassifier-100_10214-iGenes.rda"))
save(broad_return, file=paste0(study_name, "/broadClassifier100_", study_name, "_14724-iGenes.rda"))

# Validate the random forest classifier
stValSubset <- stList$validationSet
stValSubsetOrdered <- stValSubset[order(stValSubset$description1), ] #order by classification name
expValSubset <- expTrain[iGenes,rownames(stValSubsetOrdered)]
cnProc_broad <- broad_return$cnProc #select the cnProc from the broadclass training earlier 

classMatrix_broad <- broadClass_predict(cnProc_broad, expValSubset, nrand = 60)

# Rename "rand" to "Rand" for the sake of visualization
reorder <- rownames(classMatrix_broad)[c(1:12,14:15,13)]
classMatrix_broad <- classMatrix_broad[reorder,]
rownames(classMatrix_broad)[15] <- "Rand"
# Add random samples to validation heatmap
stValRand_broad <- addRandToSampTab(classMatrix_broad, stValSubsetOrdered, "description1", "sra_id")
# Rename "rand" to "Rand" in stVal
rand_ind <- which(stValRand_broad$description1 == "rand")
stValRand_broad$description1[rand_ind] <- "Rand"

grps <- as.vector(stValRand_broad$description1)
names(grps)<-rownames(stValRand_broad)

# Plot
Sys.setlocale("LC_COLLATE", "C") # So that sort() is case sensitive
pdf(file=paste0(study_name,"/Hs_TP_train100_", study_name,"_iGenes_validation_hm.pdf"), height=6, width=10)
#pdf(file="Hs-TP_train100_19004-iGenes_validation_hm.pdf", height=6, width=10)
ccn_hmClass(classMatrix_broad, grps=grps, fontsize_row=10)
dev.off()

# Classifier Assessment
assessmentDat <- ccn_classAssess(classMatrix_broad, stValRand_broad, "description1","sra_id")
pdf(file=paste0(study_name,"/Hs_TP_train100_", study_name,"_assessment_PR.pdf"), height=6, width=10)
#pdf(file="Hs-TP_train100_19004-iGenes_assessment_PR.pdf", height=6, width=10)
plot_class_PRs(assessmentDat)
dev.off()


#############################################################################
### Analyzing query studies
# Query the trained random forest classifier

#load("broadClassifier100_10215-iGenes.rda")
#load("broadClassifier100_intestine_colon_5909-iGenes.rda")

cnProc_broad <- broad_return$cnProc
classMatrixQuery <- broadClass_predict(cnProc = cnProc_broad, expDat = queryExpDat, nrand = 3)

#grp_names <- c(as.character(querySampTab$description1), rep("random", 3))
grp_names <- c(as.character(querySampTab$description3), rep("random", 3))
#grp_names <- c(as.character(querySampTab$study_id), rep("random", 3))
names(grp_names) <- c(as.character(querySampTab$sra_id), "rand_1", "rand_2", "rand_3")
#names(grp_names) <- c(as.character(querySampTab$sample_name), "rand_1", "rand_2", "rand_3")

# Re-order classMatrixQuery to match order of rows in querySampTab
classMatrixQuery <- classMatrixQuery[,names(grp_names)]

pdf_width <- ceiling(ncol(queryExpDat+3)/2) + 4
#pdf_width <- floor(ncol(queryExpDat)/26)
pdf(file=paste0(study_name,"/",study_name,"_TP_classHm.pdf"), height=4, width=pdf_width)
#pdf(file="liver-engineered-ref_TP_classHm.pdf", height=4, width=pdf_width)
acn_queryClassHm(classMatrixQuery, main = paste0("Classification Heatmap, ", study_name), 
            grps = grp_names, 
            fontsize_row=10, fontsize_col = 10, isBig = FALSE)
dev.off()



###########################################################
### Rank-based GRN Status and Network Influence Score (NIS) 
library(igraph)

# Subsetting original GRN and original TrainNormParam using iGenes

# # Load in the constructed GRN
# grnAll <- utils_loadObject("../Hs_grnAll_curatedTFs_Apr-22-2020.rda")
# # Load in the normalization parameters 
# trainNormParam <- utils_loadObject("../Hs_trainingNormalization_Apr-22-2020.rda")
# 
# 
# 
# #### SUBSET grnAll AND trainNormParam BASED ON iGenes ####
# # Subset grnTable based on iGenes
# allTargets <- grnAll$overallGRN$grnTable$TG
# newGRNTable <- grnAll$overallGRN$grnTable[which(allTargets %in% iGenes),]
# newTFsAll <- newGRNTable$TF
# newGRNTable <- newGRNTable[which(newTFsAll %in% iGenes),]
# grnAll$overallGRN$grnTable <- newGRNTable
# 
# # Subset overallGRN graph based on iGenes
# vertex_names <- V(grnAll$overallGRN$graph)$name
# graph_iGenes <- which(vertex_names %in% iGenes)
# newGraph <- induced_subgraph(graph=grnAll$overallGRN$graph, vids=graph_iGenes, impl="copy_and_delete")
# grnAll$overallGRN$graph <- newGraph
# 
# # Subset specGenes based on iGenes and tissue type
# tissueTypes <- names(grnAll$specGenes$context$general)
# newGeneral <- grnAll$specGenes$context$general
# for (tissue in tissueTypes) {
#    tissueSpecGenes <- newGeneral[[tissue]]
#    tissueSpecGenes <- tissueSpecGenes[which(names(tissueSpecGenes) %in% iGenes)]
#    newGeneral[[tissue]] <- tissueSpecGenes
# }
# grnAll$specGenes$context$general <- newGeneral
# 
# # Subset ctGRN geneLists, graphLists, and tfTargets  based on iGenes and tissue type
# grnAll$ctGRNs$geneLists <- newGeneral
# 
# newGraphLists <- grnAll$ctGRNs$graphLists
# for (tissue in tissueTypes) {
#    tissueGRN <- newGraphLists[[tissue]]
#    iVertices <- vertex_attr(tissueGRN, name="name")
#    iVertices <- iVertices[which(iVertices %in% iGenes)]
#    tissueGRN <- induced_subgraph(graph=tissueGRN, vids=iVertices, impl="copy_and_delete")
#    newGraphLists[[tissue]] <- tissueGRN
# }
# grnAll$ctGRNs$graphLists <- newGraphLists
# 
# newTFTargets <- grnAll$ctGRNs$tfTargets
# for (tissue in tissueTypes) {
#    tissueTFTargets <- newTFTargets[[tissue]]
#    tissueTFTargets <- tissueTFTargets[which(names(tissueTFTargets) %in% iGenes)]
#    for (TF in names(tissueTFTargets)) {
#       newTargets <- tissueTFTargets[[TF]]
#       newTargets <- newTargets[which(newTargets %in% iGenes)]
#       tissueTFTargets[[TF]] <- newTargets
#    }
#    newTFTargets[[tissue]] <- tissueTFTargets
# }
# grnAll$ctGRNs$tfTargets <- newTFTargets
# 
# save(grnAll, file="grnAll_neuron_studies_14724-iGenes.rda")
# 
# # Subset trainNormParam
# newTVals <- trainNormParam$tVals
# for (tissue in tissueTypes) {
#    newIndices <- which(names(newTVals[[tissue]][["mean"]]) %in% iGenes)
#    newTVals[[tissue]][["mean"]] <- newTVals[[tissue]][["mean"]][newIndices]
#    newTVals[[tissue]][["sd"]] <- newTVals[[tissue]][["sd"]][newIndices]
# }
# trainNormParam$tVals <- newTVals
# 
# save(trainNormParam, file="trainNormParam_neuron_studies_14724-iGenes.rda")

#### END OF iGenes SUBSETTING ####


##################################
load("../webApp/grnAll_liver-studies_10215-iGenes.rda")
load("../webApp/trainNormParam_liver-studies_10215-iGenes.rda")

#### Begin GRN analysis ####

queryExpDat_ranked <- logRank(queryExpDat, base = 0)

system.time(GRN_statusQuery <- ccn_queryGRNstatus(expQuery = queryExpDat_ranked, grn_return = grnAll, 
                                      trainNorm = trainNormParam, classifier_return = broad_return, prune = TRUE))

save(GRN_statusQuery, file=paste0(study_name, "/GRN_status_query_", study_name, ".rda"))
#save(GRN_statusQuery, file=paste0("liver-engineered-ref-all_GRN-status.rda"))

# Plot
cell_types <- rownames(GRN_statusQuery)
GRN_statusQuery <- GRN_statusQuery[,querySampTab$sra_id]
#GRN_statusQuery <- GRN_statusQuery[,querySampTab$sample_name]

pdf_width <- ceiling(ncol(queryExpDat)/3) + 1
pdf(file=paste0(study_name,"/GRN_status_", study_name, ".pdf"), height=8, width=pdf_width)
#pdf(file=paste0("liver-engineered-ref-all_GRN-status.pdf"), height=8, width=pdf_width)
# plot_list <- list()
# i <- 1
for(type in cell_types) {
   plot_df <-  data.frame("SampleNames" = paste(colnames(GRN_statusQuery), 
                                                #querySampTab$description3),
                                                querySampTab$description1),
                  "GRN_Status" = as.vector(GRN_statusQuery[type, ]))
   plot_df$SampleNames <- factor(plot_df$SampleNames, levels=plot_df$SampleNames)
   type_plot <- ggplot(plot_df) + geom_bar(stat="identity", data = plot_df, 
                              aes(x=SampleNames, y=GRN_Status), width = 0.7) +
                  ggtitle(paste0(type, " Network GRN Status")) + 
                  xlab("Samples") + ylab("GRN Status") + theme_bw() +
                  theme(text = element_text(size=10), 
                        legend.position="none",
                        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
                  geom_hline(yintercept=1, linetype="dashed", color = "steelblue")
   print(type_plot)
   # plot_list[[i]] <- type_plot
   # i <- i+1
}
#multiplot(plotlist=plot_list, cols=1)
dev.off()


# Calculate and plot TF scores

target_cell_type <- "liver" # CHANGE FOR SPECIFIC STUDY

system.time(TF_scores <- ccn_tfScores(expQuery = queryExpDat_ranked, grnAll = grnAll, trainNorm = trainNormParam,
                          classifier_return = broad_return, subnetName = target_cell_type,
                          exprWeight = FALSE, normTFscore = TRUE))
save(TF_scores, file=paste0(study_name, "/TF_scores_", target_cell_type, "_", study_name, ".rda"))

pdf_width <- floor(nrow(TF_scores)/5)
#pdf(file=paste0(study_name,"/NIS-plot_top_", target_cell_type, "_", study_name, ".pdf"), height=6, width=pdf_width)
pdf(file=paste0(study_name,"/NIS-plot_", target_cell_type, "_", study_name, ".pdf"), height=6, width=pdf_width)
sample_names <- querySampTab$sra_id
#sample_names <- querySampTab$sample_name
for(sample in sample_names) {
   #descript <- querySampTab$description3[which(querySampTab$sra_id == sample)]
   descript <- querySampTab$description1[which(querySampTab$sra_id == sample)]
   #descript <- querySampTab$description3[which(querySampTab$sample_name == sample)]
   plot_df <- data.frame("TFs" = rownames(TF_scores),
                  "Scores" = as.vector(TF_scores[,sample]))
   sample_TFplot <- ggplot(plot_df, aes(x = TFs , y = Scores)) + geom_bar(stat="identity") + #aes(fill = medVal)) +
      theme_bw() + #scale_fill_gradient2(low = "purple", 
                                                         # mid = "white", 
                                                         # high = "orange") + 
      ggtitle(paste0(sample, ", ", descript, ", ", target_cell_type, " transcription factor scores")) +
      ylab("Network influence score") + xlab("Transcriptional regulator") + 
      theme(legend.position = "none", axis.text = element_text(size = 8)) +
      theme(text = element_text(size=10), 
                   legend.position="none",
                   axis.text.x = element_text(angle = 45, vjust=0.5))
      #coord_flip()
   print(sample_TFplot)
}
dev.off()



