# install.packages("devtools")
library(devtools)
# install_github("pcahan1/CellNet", ref="master")
# install_github("pcahan1/cancerCellNet@v0.1.1", ref="master")
library(CellNet)
# install.packages("plyr")
# install.packages("ggplot2")
# install.packages("RColorBrewer")
# install.packages("pheatmap")
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
source("plotting.R")


expTrain <- utils_loadObject("Hs_expTrain_Jun-20-2017.rda")
stTrain <- utils_loadObject("Hs_stTrain_Jun-20-2017.rda")

# Prepare training data
set.seed(99) # Setting a seed for the random number generator allows us to reproduce the same split in the future
stList <- splitCommon_proportion(sampTab = stTrain, proportion = 0.66, dLevel = "description1") # Use 2/3 of training data for training and 1/3 for validation
stTrainSubset <- stList$trainingSet
expTrainSubset <- expTrain[,rownames(stTrainSubset)]
stValSubset <- stList$validationSet
expValSubset <- expTrain[rownames(stValSubset)]

# Train classifier
system.time(my_classifier <- broadClass_train(stTrain = stTrainSubset, 
                                expTrain = expTrainSubset, 
                                colName_cat = "description1", 
                                colName_samp = "sra_id", 
                                nRand = 70,
                                nTopGenes = 100, 
                                nTopGenePairs = 100, 
                                nTrees = 2000, 
                                stratify=TRUE, 
                                sampsize=25, 
                                quickPairs=TRUE)) # Increasing the number of top genes and top gene pairs increases the resolution of the classifier but increases the computing time
save(my_classifier file="cellnet_classifier_100topGenes_100genePairs.rda")

# Classifier Validation
stValSubsetOrdered <- stValSubset[order(stValSubset$description1), ] #order samples by classification name
expValSubset <- expValSubset[rownames(stValSubsetOrdered)]
cnProc <- my_classifier$cnProc #select the cnProc from the earlier class training

classMatrix <- broadClass_predict(cnProc, expValSubset, nrand = 60)
stValRand <- addRandToSampTab(classMatrix, stValSubsetOrdered, desc="description1", id="sra_id")

grps <- as.vector(stValRand$description1)
names(grps)<-rownames(stValRand)

# Plot validation heatmap
ccn_hmClass(classMatrix, grps=grps, fontsize_row=10)
dev.off()

# Plot validation PR curves
assessmentDat <- ccn_classAssess(classMatrix, stValRand, classLevels="description1", dLevelSID="sra_id")
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

pdf(file="validation_gene-pair_comparison.pdf", width=10, height=80)
plotGeneComparison(geneCompareMatrix, grps = val_grps, fontsize_row = 6)
dev.off()

# Create and save xpairs_list object for grn reconstruction and training normalization parameters:
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

save(xpairs_list, file="Hs_xpairs_list.rda")


### Classify query samples
queryExpDat <- utils_loadObject("example_queryExpMat_counts.rda")
querySampTab <- utils_loadObject("example_sampTab.rda")
queryExpDat <- log(1+queryExpDat)

classMatrixQuery <- broadClass_predict(cnProc = cnProc, expDat = queryExpDat, nrand = 3) 
grp_names <- c(as.character(querySampTab$description1), rep("random", 3))
names(grp_names) <- c(as.character(rownames(querySampTab)), "rand_1", "rand_2", "rand_3")

# Re-order classMatrixQuery to match order of rows in querySampTab
classMatrixQuery <- classMatrixQuery[,names(grp_names)]
save(classMatrixQuery, file="example_classificationMatrix.rda"))

# Plot classification heatmap
pdf(file="query_classification_heatmap.pdf"), height=4)
# This function is in plotting.R
acn_queryClassHm(classMatrixQuery, main = paste0("Classification Heatmap, ", study_name), 
            grps = grp_names, 
            fontsize_row=10, fontsize_col = 10, isBig = FALSE)
dev.off()


### GRN Status
grnAll <- utils_loadObject("Hs_grnAll_curatedTFs_Apr-22-2020.rda")
trainNormParam <- utils_loadObject("Hs_trainingNormalization_Apr-22-2020.rda")


# Compute GRN statuses
queryExpDat_ranked <- logRank(queryExpDat, base = 0)

system.time(GRN_statusQuery <- ccn_queryGRNstatus(expQuery = queryExpDat_ranked, grn_return = grnAll, 
                                      trainNorm = trainNormParam, classifier_return = my_classifier, prune = TRUE))
save(GRN_statusQuery, file="my_study_GRN_status.rda")

# Plot GRN status barplots
cell_types <- rownames(GRN_statusQuery)
GRN_statusQuery <- GRN_statusQuery[,rownames(querySampTab)]
#GRN_statusQuery <- GRN_statusQuery[,querySampTab$sample_name]

pdf_width <- ceiling(ncol(queryExpDat)/3) + 1
pdf(file="my_study_GRN_status_plots.pdf", height=8, width=pdf_width)

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


### Network Influence Score (Transcriptional regulator scoring)

target_cell_type <- "my_cell_type" # CHANGE FOR SPECIFIC CONTEXT
system.time(TF_scores <- ccn_tfScores(expQuery = queryExpDat_ranked, grnAll = grnAll, trainNorm = trainNormParam,
                          classifier_return = broad_return, subnetName = target_cell_type,
                          exprWeight = FALSE, normTFscore = TRUE))
save(TF_scores, file="my_study_TF_scores.rda")

# Choose top scoring 25 TFs for plotting
TFsums <- abs(rowSums(TF_scores))
ordered_TFsums <- TFsums[order(TFsums, decreasing = TRUE)]
if(length(TFsums) > 25) {
    top_display_TFs <- names(ordered_TFsums)[1:25]    
} else {
    top_display_TFs <- names(ordered_TFsums)
}
TF_scores <- TF_scores[top_display_TFs,]


# Plot TF scores
sample_names <- rownames(querySampTab)

pdf(file="my_study_TF_scores_my_cell_type.pdf", height=6, width=8)
for(sample in sample_names) {
   descript <- querySampTab$description1[which(rownames(querySampTab) == sample)]
   plot_df <- data.frame("TFs" = rownames(TF_scores),
                         "Scores" = as.vector(TF_scores[,sample]))
   sample_TFplot <- ggplot(plot_df, aes(x = TFs , y = Scores)) + geom_bar(stat="identity") + #aes(fill = medVal)) +
      theme_bw() + 
      #scale_fill_gradient2(low = "purple", 
      # mid = "white", 
      # high = "orange") + 
      ggtitle(paste0(sample, ", ", descript, ", ", target_cell_type, " transcription factor scores")) +
      ylab("Network influence score") + xlab("Transcriptional regulator") + 
      theme(legend.position = "none", axis.text = element_text(size = 8)) +
      theme(text = element_text(size=10), 
            legend.position="none",
            axis.text.x = element_text(angle = 45, vjust=0.5))
   print(sample_TFplot)
}
dev.off()


