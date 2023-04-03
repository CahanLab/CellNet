# Installation
install.packages("devtools")
library(devtools)
install_github("pcahan1/CellNet", ref="master")
install_github("pcahan1/cancerCellNet@v0.1.1", ref="master")

# Load libraries
library(CellNet)
library(cancerCellNet)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(plotly)
library(igraph)
source("pacnet_utils.R")

# Load training data:
expTrain <- utils_loadObject("Hs_expTrain_Jun-20-2017.rda")
stTrain <- utils_loadObject("Hs_stTrain_Jun-20-2017.rda")

# Load engineered reference data and query data. 
# We need to load these at this point to identify genes found in across all datasets.
liverRefExpDat <- utils_loadObject("liver_engineeredRef_normalized_expDat_all.rda")
liverRefSampTab <- utils_loadObject("liver_engineeredRef_sampTab_all.rda")
queryExpDat <- read.csv("example_data/example_counts_matrix.csv", row.names=1)
querySampTab <- read.csv("example_data/example_sample_metadata_table.csv")
rownames(querySampTab) <- querySampTab$sample_name
study_name <- "liver_example"

# Identify intersecting genes
iGenes <- intersect(rownames(expTrain), rownames(liverRefExpDat))
iGenes <- intersect(iGenes, rownames(queryExpDat))

expTrain <- expTrain[iGenes,]

# Split the training data into a training subset and a validation subset:
set.seed(99) # Setting a seed for the random number generator allows us to reproduce the same split in the future
stList <- splitCommon_proportion(sampTab = stTrain, proportion = 0.66, dLevel = "description1") # Use 2/3 of training data for training and 1/3 for validation
stTrainSubset <- stList$trainingSet
expTrainSubset <- expTrain[,rownames(stTrainSubset)]

# See number of samples of each unique type in description1 in training subset
table(stTrainSubset$description1)

stValSubset <- stList$validationSet
expValSubset <- expTrain[,rownames(stValSubset)]
# See number of samples of each unique type in description1 in validation subset
table(stValSubset$description1)

# Train the random forest classifier, takes 3-10 minutes depending on memory availability:
system.time(my_classifier <- broadClass_train(stTrain = stTrainSubset, 
                                              expTrain = expTrainSubset, 
                                              colName_cat = "description1", 
                                              colName_samp = "sra_id", 
                                              nRand = 70,
                                              nTopGenes = 100, 
                                              nTopGenePairs = 100, 
                                              nTrees = 2000, 
                                              stratify=TRUE, 
                                              sampsize=25, # Must be less than the smallest n in table(stTrainSubset$description1)
                                              quickPairs=TRUE)) # Increasing the number of top genes and top gene pairs increases the resolution of the classifier but increases the computing time
save(my_classifier, file="example_outputs/cellnet_classifier_100topGenes_100genePairs.rda")

# Classifier Validation

stValSubsetOrdered <- stValSubset[order(stValSubset$description1), ] #order samples by classification name
expValSubset <- expValSubset[,rownames(stValSubsetOrdered)]
cnProc <- my_classifier$cnProc #select the cnProc from the earlier class training

classMatrix <- broadClass_predict(cnProc, expValSubset, nrand = 60) 
stValRand <- addRandToSampTab(classMatrix, stValSubsetOrdered, desc="description1", id="sra_id")

grps <- as.vector(stValRand$description1)
names(grps)<-rownames(stValRand)

# Create validation heatmap
png(file="classification_validation_hm.png", height=6, width=10, units="in", res=300)
ccn_hmClass(classMatrix, grps=grps, fontsize_row=10)
dev.off()

# Plot validation precision-recall curves:
assessmentDat <- ccn_classAssess(classMatrix, stValRand, classLevels="description1", dLevelSID="sra_id")
png(file="example_outputs/classifier_assessment_PR.png", height=8, width=10, units="in", res=300)
plot_class_PRs(assessmentDat)
dev.off()

# Gene pair validation
genePairs <- cnProc$xpairs
# Get gene to gene comparison of each gene pair in the expression table
expTransform <- query_transform(expTrainSubset, genePairs)
avgGenePair_train <- avgGeneCat(expDat = expTransform, sampTab = stTrainSubset, 
                                dLevel = "description1", sampID = "sra_id")

genePairs_val <- query_transform(expValSubset, genePairs)
geneCompareMatrix <- makeGeneCompareTab(queryExpTab = genePairs_val,
                                        avgGeneTab = avgGenePair_train, geneSamples = genePairs)
val_grps <- stValSubset[,"description1"]
val_grps <- c(val_grps, colnames(avgGenePair_train))
names(val_grps) <- c(rownames(stValSubset), colnames(avgGenePair_train))

png(file="example_outputs/validation_gene-pair_comparison.png", width=10, height=80, units="in", res=300)
plotGeneComparison(geneCompareMatrix, grps = val_grps, fontsize_row = 6)
dev.off()

# Create and save xpairs_list object
xpairs_list <- vector("list", 14) 
for (pair in rownames(avgGenePair_train)) {
   for (j in 1:ncol(avgGenePair_train)) {
      if (avgGenePair_train[pair,j] >= 0.5) {
         if (is.null(xpairs_list[[j]])) {
            xpairs_list[[j]] <- c(pair)
         } else { 
            xpairs_list[[j]] <- c(xpairs_list[[j]], pair)
         }
      }  
   }
}
xpair_names <- colnames(avgGenePair_train)
xpair_names <- sub(pattern="_Avg", replacement="", x=xpair_names)
names(xpairs_list) <- xpair_names

for (type in names(xpairs_list)) {
   names(xpairs_list[[type]]) <- xpairs_list[[type]]
}
save(xpairs_list, file="example_outputs/Hs_xpairs_list.rda")

################################################
# Querying the classifier

# Classify engineered reference panel samples
classMatrixLiverRef <- broadClass_predict(cnProc = cnProc, expDat = liverRefExpDat, nrand = 10) 
grp_names1 <- c(as.character(liverRefSampTab$description1), rep("random", 10))
names(grp_names1) <- c(as.character(rownames(liverRefSampTab)), paste0("rand_", c(1:10)))
# Re-order classMatrixQuery to match order of rows in querySampTab
classMatrixLiverRef <- classMatrixLiverRef[,names(grp_names1)]

png(file="example_outputs/heatmapLiverRef.png", height=12, width=9, units="in", res=300)
heatmapRef(classMatrixLiverRef, liverRefSampTab) # This function can be found in pacnet_utils.R
dev.off()

# Alternatively, for an interactive plotly version:
heatmapPlotlyRef(classMatrixLiverRef, liverRefSampTab)


# Classify query samples
# Perform log transform:
queryExpDat <- log(1+queryExpDat)

classMatrixQuery <- broadClass_predict(cnProc = cnProc, expDat = queryExpDat, nrand = 3) 
grp_names <- c(as.character(querySampTab$description1), rep("random", 3))
names(grp_names) <- c(as.character(rownames(querySampTab)), paste0("rand_", c(1:3)))
# Re-order classMatrixQuery to match order of rows in querySampTab
classMatrixQuery <- classMatrixQuery[,names(grp_names)]
save(classMatrixQuery, file="example_outputs/example_classificationMatrix.rda")

#Plot classification heatmap:
png(file="example_outputs/query_classification_heatmap.png", height=4, width=8, units="in", res=300)
# This function can be found in pacnet_utils.R
acn_queryClassHm(classMatrixQuery, main = paste0("Classification Heatmap, ", study_name), 
                 grps = grp_names, 
                 fontsize_row=10, fontsize_col = 10, isBig = FALSE)
dev.off()


# Compute GRN Status
#Subset `grnAll` and `trainNormParam` objects based on intersecting genes.
grnAll <- utils_loadObject("liver_grnAll.rda")
trainNormParam <- utils_loadObject("liver_trainNormParam.rda")
# These two functions can be found in pacnet_utils.R
grnAll <- subsetGRNall(grnAll, iGenes)
trainNormParam <- subsetTrainNormParam(trainNormParam, grnAll, iGenes)

#Compute GRN statuses and save:
queryExpDat_ranked <- logRank(queryExpDat, base = 0)
queryExpDat_ranked <- as.data.frame(queryExpDat_ranked)
system.time(GRN_statusQuery <- ccn_queryGRNstatus(expQuery = queryExpDat_ranked, grn_return = grnAll, 
                                                  trainNorm = trainNormParam, classifier_return = my_classifier, prune = TRUE))
save(GRN_statusQuery, file="example_outputs/my_study_GRN_status.rda")

# Plot GRN status bar plots:
cell_types <- rownames(GRN_statusQuery)
# GRN_statusQuery <- GRN_statusQuery[,rownames(querySampTab)]
GRN_statusQuery <- GRN_statusQuery[,querySampTab$sample_name]
pdf_width <- ceiling(ncol(queryExpDat)/3) + 1

pdf(file="example_outputs/my_study_GRN_status_plots.pdf", height=8, width=pdf_width)
plot_list <- list()
i <- 1
for(type in cell_types) {
   plot_df <-  data.frame("SampleNames" = paste(colnames(GRN_statusQuery), querySampTab$description1),
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
}
dev.off()


# Compute Network Influence Score (NIS) for transcriptional regulators
target_cell_type <- "liver" # CHANGE FOR SPECIFIC CONTEXT
system.time(TF_scores <- pacnet_nis(expDat = queryExpDat_ranked, stQuery=querySampTab, iGenes=iGenes,
                                    grnAll = grnAll, trainNorm = trainNormParam,
                                    subnet = target_cell_type, ctt=target_cell_type,
                                    colname_sid="sample_name", relaWeight=0))

save(TF_scores, file="example_outputs/my_study_TF_scores.rda")

# Choose top-scoring 25 TFs for plotting:
TFsums <- rowSums(abs(TF_scores))
ordered_TFsums <- TFsums[order(TFsums, decreasing = TRUE)]
if(length(TFsums) > 25) {
   top_display_TFs <- names(ordered_TFsums)[1:25]    
} else {
   top_display_TFs <- names(ordered_TFsums)
}
TF_scores <- TF_scores[top_display_TFs,]

#Plot TF scores:
sample_names <- rownames(querySampTab)

pdf(file="example_outputs/my_study_TF_scores_my_cell_type.pdf", height=6, width=8)
for(sample in sample_names) {
   descript <- querySampTab$description1[which(rownames(querySampTab) == sample)]
   plot_df <- data.frame("TFs" = rownames(TF_scores),
                         "Scores" = as.vector(TF_scores[,sample]))
   sample_TFplot <- ggplot(plot_df, aes(x = reorder(TFs,Scores,mean) , y = Scores)) + 
     geom_bar(stat="identity") + #aes(fill = medVal)) +
      theme_bw() + 
      ggtitle(paste0(sample, ", ", descript, ", ", target_cell_type, " transcription factor scores")) +
      ylab("Network influence score") + xlab("Transcriptional regulator") + 
      theme(legend.position = "none", axis.text = element_text(size = 8)) +
      theme(text = element_text(size=10), 
            legend.position="none",
            axis.text.x = element_text(angle = 45, vjust=0.5))
   print(sample_TFplot)
}
dev.off()

