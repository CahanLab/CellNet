# Training and Running Platform-Agnostic CellNet

### Table of contents

[Introduction](#introduction)

[Data](#data)

[Installation](#installation)

[Training](#training)

[Query](#query)

---

### Introduction <a name="introduction"></a>
This is a walk-through tutorial on 1) how to train a PACNet classifier using a preprocessed training expression matrix and 2) how to apply the classifier to a preprocessed query study. All the code below is also available in a single R file, `agnosticCellNet_example.R`.

---

### Data <a name="data"></a>

Training Data

| Species | Metadata table | Expression matrix | Grn Object | Training Parameters |
|---------|----------------|-------------------|------------|---------------------|
| Human   | [Hs_stTrain_Jun-20-2017.rda](https://s3.console.aws.amazon.com/s3/object/cellnet-rnaseq?region=us-east-1&prefix=ref/cnproc/HS/Hs_stTrain_Jun-20-2017.rda) | [Hs_expTrain_Jun-20-2017.rda](https://s3.console.aws.amazon.com/s3/object/cellnet-rnaseq?region=us-east-1&prefix=ref/cnproc/HS/Hs_expTrain_Jun-20-2017.rda) | [Hs_grnAll_curatedTFs_Apr-22-2020.rda](https://s3.console.aws.amazon.com/s3/object/cellnet-rnaseq?region=us-east-1&prefix=ref/cnproc/HS/Hs_grnAll_curatedTFs_Apr-22-2020.rda) | [Hs_trainingNormalization_Apr-22-2020.rda](https://s3.console.aws.amazon.com/s3/object/cellnet-rnaseq?region=us-east-1&prefix=ref/cnproc/HS/Hs_trainingNormalization_Apr-22-2020.rda) |
| Mouse   | [Mm_stTrain_Oct-24-2016.rda](https://s3.console.aws.amazon.com/s3/object/cellnet-rnaseq?region=us-east-1&prefix=ref/cnproc/MM/Mm_stTrain_Oct-24-2016.rda) | [Mm_expTrain_Oct-24-2016.rda](https://s3.console.aws.amazon.com/s3/object/cellnet-rnaseq?region=us-east-1&prefix=ref/cnproc/MM/Mm_expTrain_Oct-24-2016.rda) | [Mm_grnAll_curatedTFs_Apr-22-2020.rda](https://s3.console.aws.amazon.com/s3/object/cellnet-rnaseq?region=us-east-1&prefix=ref/cnproc/MM/Mm_grnAll_curatedTFs_Apr-22-2020.rda) | [Mm_trainingNormalization_Apr-22-2020.rda](https://s3.console.aws.amazon.com/s3/object/cellnet-rnaseq?region=us-east-1&prefix=ref/cnproc/MM/Mm_trainingNormalization_Apr-22-2020.rda) |


Human Engineered Reference Panels

| Cell type | Metadata table | Expression matrix | Trained Classifier | Grn Object Subset | Training Parameter Subset |
|-----------|----------------|-------------------|--------------------|-------------------|---------------------------|
| Heart | [heart_engineeredRef_sampTab_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/heart_engineeredRef_sampTab_all.rda) | [heart_engineeredRef_normalized_expDat_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/heart_engineeredRef_normalized_expDat_all.rda) | [heart_broadClassifier100.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/heart_broadClassifier100.rda) | [heart_grnAll.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/heart_grnAll.rda) | [heart_trainNormParam.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/heart_trainNormParam.rda) |
| HSPC | [hspc_engineeredRef_sampTab_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/hspc_engineeredRef_sampTab_all.rda) | [hspc_engineeredRef_normalized_expDat_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/hspc_engineeredRef_normalized_expDat_all.rda) | [hspc_broadClassifier100.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/hspc_broadClassifier100.rda) | [hspc_grnAll.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/hspc_grnAll.rda) | [hspc_trainNormParam.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/hspc_trainNormParam.rda) |
| Intestine/colon | [intestine_colon_engineeredRef_sampTab_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/intestine_colon_engineeredRef_sampTab_all.rda) | [intestine_colon_engineeredRef_normalized_expDat_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/intestine_colon_engineeredRef_normalized_expDat_all.rda) | [intestine_colon_broadClassifier100.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/intestine_colon_broadClassifier100.rda) | [intestine_colon_grnAll.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/intestine_colon_grnAll.rda) | [intestine_colon_trainNormParam.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/intestine_colon_trainNormParam.rda) |
| Liver | [liver_engineeredRef_sampTab_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/liver_engineeredRef_sampTab_all.rda) | [liver_engineeredRef_normalized_expDat_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/liver_engineeredRef_normalized_expDat_all.rda) | [broadClassifier100.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/liver_broadClassifier100.rda) | [liver_grnAll.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/liver_grnAll.rda) | [liver_trainNormParam.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/liver_trainNormParam.rda) |
| Lung | [lung_engineeredRef_sampTab_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/lung_engineeredRef_sampTab_all.rda) | [lung_engineeredRef_normalized_expDat_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/lung_engineeredRef_normalized_expDat_all.rda) | [lung_broadClassifier100.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/lung_broadClassifier100.rda) | [lung_grnAll.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/lung_grnAll.rda) | [lung_trainNormParam.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/lung_trainNormParam.rda) |
| Neuron | [neuron_engineeredRef_sampTab_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/neuron_engineeredRef_sampTab_all.rda) | [neuron_engineeredRef_normalized_expDat_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/neuron_engineeredRef_normalized_expDat_all.rda) | [neuron_broadClassifier100.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/neuron_broadClassifier100.rda) | [neuron_grnAll.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/neuron_grnAll.rda) | [neuron_trainNormParam.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/neuron_trainNormParam.rda) |
| Skeletal muscle | [skeletal_muscle_engineeredRef_sampTab_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/skeletal_muscle_engineeredRef_sampTab_all.rda) | [skeletal_muscle_engineeredRef_normalized_expDat_all.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/skeletal_muscle_engineeredRef_normalized_expDat_all.rda) | [skeletal_muscle_broadClassifier100.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/skeletal_muscle_broadClassifier100.rda) | [skeletal_muscle_grnAll.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/skeletal_muscle_grnAll.rda) | [skeletal_muscle_trainNormParam.rda](https://s3.console.aws.amazon.com/s3/object/cnobjects?region=us-east-1&prefix=webApps/agnosticCellNet_web/skeletal_muscle_trainNormParam.rda) |

---

### Installation <a name="installation"></a>
```R
install.packages("devtools")
library(devtools)
install_github("pcahan1/CellNet", ref="master")
install_github("pcahan1/cancerCellNet@v0.1.1", ref="master")
```
Other required packages: plyr, ggplot2, RColorBrewer, pheatmap, plotly

#### Prerequisites 
Load required R packages.
```R
library(CellNet)
library(cancerCellNet)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(plotly)
source("pacnet_utils.R")
```

### Training <a name="training"></a>

#### Confirm correct format of training data expression matrix and training sample metadata table
Expression matrix should have gene symbols as row names and sample names as column names. Sample metadata table should have sample names as row names and sample features as column names. Column names of expression matrix must match row names of metadata table. (See the [example_data folder](example_data/) for a small example of an expression matrix and metadata table.

For classifier training to be robust, there should be at least 60 independent replicates per training type.

#### Begin Training

Load training data:
```R
expTrain <- utils_loadObject("Hs_expTrain_Jun-20-2017.rda") 
stTrain <- utils_loadObject("Hs_stTrain_Jun-20-2017.rda")
```

Load engineered reference data and query data. We need to load these at this point to identify genes found in across all datasets.
```R
liverRefExpDat <- utils_loadObject("liver_engineeredRef_normalized_expDat_all.rda")
liverRefSampTab <- utils_loadObject("liver_engineeredRef_sampTab_all.rda")
queryExpDat <- read.csv("example_counts_matrix.csv", row.names=1)
querySampTab <- read.csv("example_sample_metadata_table.csv", row.names=1)
```

Identify intersecting genes:
```R
iGenes <- intersect(rownames(expTrain), rownames(liverRefExpDat))
iGenes <- intersect(iGenes, rownames(queryExpDat))

expTrain <- expTrain[iGenes,]
```

Split the training data into a training subset and a validation subset:
```R
set.seed(99) # Setting a seed for the random number generator allows us to reproduce the same split in the future
stList <- splitCommon_proportion(sampTab = stTrain, proportion = 0.66, dLevel = "description1") # Use 2/3 of training data for training and 1/3 for validation
stTrainSubset <- stList$trainingSet
expTrainSubset <- expTrain[,rownames(stTrainSubset)]

#See number of samples of each unique type in description1 in training subset
table(stTrainSubset$description1)

stValSubset <- stList$validationSet
expValSubset <- expTrain[rownames(stValSubset)]
#See number of samples of each unique type in description1 in validation subset
table(stValSubset$description1)

```

Train the random forest classifier, takes 3-10 minutes depending on memory availability:
```R
system.time(my_classifier <- broadClass_train(stTrain = stTrainSubset, 
                                expTrain = expTrainSubset, 
                                colName_cat = "description1", 
                                colName_samp = "sra_id", 
                                nRand = 70, # Must be less than the smallest number in table(stTrainSubset$description1)
                                nTopGenes = 100, 
                                nTopGenePairs = 100, 
                                nTrees = 2000, 
                                stratify=TRUE, 
                                sampsize=25, # Must be less than the smallest number in table(stTrainSubset$description1)
                                quickPairs=TRUE)) # Increasing the number of top genes and top gene pairs increases the resolution of the classifier but increases the computing time
save(my_classifier, file="cellnet_classifier_100topGenes_100genePairs.rda")
```

#### Classifier Validation
 
Plot validation heatmap:
```R
stValSubsetOrdered <- stValSubset[order(stValSubset$description1), ] #order samples by classification name
expValSubset <- expValSubset[rownames(stValSubsetOrdered)]
cnProc <- my_classifier$cnProc #select the cnProc from the earlier class training

classMatrix <- broadClass_predict(cnProc, expValSubset, nrand = 60) #nrand must be less than the smallest number in table(stValSubset$description1)
stValRand <- addRandToSampTab(classMatrix, stValSubsetOrdered, desc="description1", id="sra_id")

grps <- as.vector(stValRand$description1)
names(grps)<-rownames(stValRand)

# Create pdf of validation heatmap
pdf(file="classification_validation_hm.pdf", height=6, width=10)
ccn_hmClass(classMatrix, grps=grps, fontsize_row=10)
dev.off()
```

![Example validation heatmap](example_plots/classifier_validation_heatmap.pdf)


Plot validation precision-recall curves:
```R
assessmentDat <- ccn_classAssess(classMatrix, stValRand, classLevels="description1", dLevelSID="sra_id")

pdf(file="classifier_assessment_PR.pdf", height=8, width=10)
plot_class_PRs(assessmentDat)
dev.off()
```

![Example PR plots](example_plots/classifier_precision_recall.pdf)


#### Gene pair validation
 
```R
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
```


Create and save xpairs_list object for grn reconstruction and training normalization parameters:
```R
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
```

---

### Querying the classifier <a name="query"></a>

#### Classify engineered reference panel samples

```R
classMatrixLiverRef <- broadClass_predict(cnProc = cnProc, expDat = liverRefExpDat, nrand = 10) 
grp_names1 <- c(as.character(liverRefSampTab$description1), rep("random", 10))
names(grp_names1) <- c(as.character(rownames(liverRefSampTab)), paste0("rand_", c(1:10)))

# Re-order classMatrixQuery to match order of rows in querySampTab
classMatrixLiverRef <- classMatrixLiverRef[,names(grp_names1)]
```

Plot classification heatmap:
```R
pdf("heatmapLiverRef.pdf", height=12, width=9)
# This function can be found in pacnet_utils.R
heatmapRef(classMatrixLiverRef, liverRefSampTab)
dev.off()

# Alternatively, for an interactive plotly version:
heatmapPlotlyRef(classMatrixLiverRef, liverRefSampTab)

```

![Example engineered reference panel plotly heatmap](example_plots/heatmapLiverRef.pdf)

#### Classify query samples

Perform log transform:
```R
queryExpDat <- log(1+queryExpDat)
```

Classify query samples:
```R
classMatrixQuery <- broadClass_predict(cnProc = cnProc, expDat = queryExpDat, nrand = 3) 
grp_names <- c(as.character(querySampTab$description1), rep("random", 3))
names(grp_names) <- c(as.character(rownames(querySampTab)), paste0("rand_", c(1:10)))

# Re-order classMatrixQuery to match order of rows in querySampTab
classMatrixQuery <- classMatrixQuery[,names(grp_names)]

save(classMatrixQuery, file="example_classificationMatrix.rda"))
```

Plot classification heatmap:
```R
pdf(file="query_classification_heatmap.pdf"), height=4)

# This function can be found in pacnet_utils.R
acn_queryClassHm(classMatrixQuery, main = paste0("Classification Heatmap, ", study_name), 
            grps = grp_names, 
            fontsize_row=10, fontsize_col = 10, isBig = FALSE)
dev.off()
```

#### Compute GRN Status

Subset `grnAll` and `trainNormParam` objects based on intersecting genes.
```R
grnAll <- utils_loadObject("liver_grnAll.rda")
trainNormParam <- utils_loadObject("liver_trainNormParam.rda")


# These two functions can be found in pacnet_utils.R
grnAll <- subsetGRNall(grnAll, iGenes)
trainNormParam <- subsetTrainNormParam(trainNormParam, grnAll, iGenes)
```

Compute GRN statuses and save:
```R
queryExpDat_ranked <- logRank(queryExpDat, base = 0)

system.time(GRN_statusQuery <- ccn_queryGRNstatus(expQuery = queryExpDat_ranked, grn_return = grnAll, 
                                      trainNorm = trainNormParam, classifier_return = my_classifier, prune = TRUE))
save(GRN_statusQuery, file="my_study_GRN_status.rda")
```

Plot GRN status bar plots:
```R
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
```

[Example GRN status plots](example_plots/example_GRN_status.pdf)


#### Compute Network Influence Score (NIS) for transcriptional regulators

Compute and save TF scores:
```R
target_cell_type <- "my_cell_type" # CHANGE FOR SPECIFIC CONTEXT
system.time(TF_scores <- ccn_tfScores(expQuery = queryExpDat_ranked, grnAll = grnAll, trainNorm = trainNormParam,
                          classifier_return = broad_return, subnetName = target_cell_type,
                          exprWeight = FALSE, normTFscore = TRUE))
save(TF_scores, file="my_study_TF_scores.rda")
```

Choose top scoring 25 TFs for plotting:
```R
TFsums <- rowSums(abs(TF_scores))
ordered_TFsums <- TFsums[order(TFsums, decreasing = TRUE)]
if(length(TFsums) > 25) {
    top_display_TFs <- names(ordered_TFsums)[1:25]    
} else {
    top_display_TFs <- names(ordered_TFsums)
}
TF_scores <- TF_scores[top_display_TFs,]
```

Plot TF scores:
```R
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
```

[Example NIS plots](example_plots/example_NIS.pdf)

Fin.

---
