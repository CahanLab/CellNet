# Training and Running Platform-Agnostic CellNet

---

## Introduction
This is a walk-through tutorial on 1) how to train a platform-agnostic CellNet classifier using an already preprocessed (i.e. converted from fastqs to an expression matrix) bulk training dataset and 2) how to apply the classifier to an already preprocessed bulk query study. All the code below is also available in a single R file, `agnosticCellNet_example.R`.

---

#### Prerequisite installations
First, install and/or load required R packages. (Uncomment lines for packages that are not yet installed.)
```R
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
```

### Training

#### Confirm correct format of training data expression matrix and training sample metadata table
Expression matrix should have gene symbols as row names and sample names as column names. Sample metadata table should have sample names as row names and sample features as column names. Column names of expression matrix must match row names of metadata table. (See the [example_data folder](example_data/) for a small example of an expression matrix and metadata table. Load .rda files in R with `load("my_example_file.rda")` ).

For classifier training to be robust, there should be at least 60 independent replicates per training type.

#### Begin Training

Load training data:
```R
expTrain <- utils_loadObject("Hs_expTrain_Jun-20-2017.rda")
stTrain <- utils_loadObject("Hs_stTrain_Jun-20-2017.rda")
```

Split the training data into a training subset and a validation subset:
```R
set.seed(99) # Setting a seed for the random number generator allows us to reproduce the same split in the future
stList <- splitCommon_proportion(sampTab = stTrain, proportion = 0.66, dLevel = "description1") # Use 2/3 of training data for training and 1/3 for validation
stTrainSubset <- stList$trainingSet
expTrainSubset <- expTrain[,rownames(stTrainSubset)]
stValSubset <- stList$validationSet
expValSubset <- expTrain[rownames(stValSubset)]
```

Train the random forest classifier, takes 3-10 minutes depending on memory availability:
```R
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
```

#### Classifier Validation
 
Plot validation heatmap:
```R
stValSubsetOrdered <- stValSubset[order(stValSubset$description1), ] #order samples by classification name
expValSubset <- expValSubset[rownames(stValSubsetOrdered)]
cnProc <- my_classifier$cnProc #select the cnProc from the earlier class training

classMatrix <- broadClass_predict(cnProc, expValSubset, nrand = 60)
stValRand <- addRandToSampTab(classMatrix, stValSubsetOrdered, desc="description1", id="sra_id")

grps <- as.vector(stValRand$description1)
names(grps)<-rownames(stValRand)

# Create pdf of validation heatmap
pdf(file="classification_validation_hm.pdf", height=6, width=10)
ccn_hmClass(classMatrix, grps=grps, fontsize_row=10)
dev.off()
```

![Validation heatmap](example_plots/classifier_validation_heatmap.pdf)


Plot validation precision-recall curves:
```R
assessmentDat <- ccn_classAssess(classMatrix, stValRand, classLevels="description1", dLevelSID="sra_id")

pdf(file="classifier_assessment_PR.pdf", height=8, width=10)
plot_class_PRs(assessmentDat)
dev.off()
```

![PR plots](example_plots/classifier_precision_recall.pdf)


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

### Querying the classifier

#### Check format of query expression matrix and query sample table
Similarly to training data, expression matrix should have gene symbols as row names and sample names as column names. Sample metadata table should have sample names as row names and sample features as column names. Column names of expression matrix must match row names of metadata table.

#### Classify query samples

Load query data:
```R
queryExpDat <- utils_loadObject("example_queryExpMat_counts.rda")
querySampTab <- utils_loadObject("example_sampTab.rda")
```

If query expression matrix is in the form of counts, perform log transform:
```R
queryExpDat <- log(1+queryExpDat)
```

Classify query samples:
```R
classMatrixQuery <- broadClass_predict(cnProc = cnProc, expDat = queryExpDat, nrand = 3) 
grp_names <- c(as.character(querySampTab$description1), rep("random", 3))
names(grp_names) <- c(as.character(rownames(querySampTab)), "rand_1", "rand_2", "rand_3")

# Re-order classMatrixQuery to match order of rows in querySampTab
classMatrixQuery <- classMatrixQuery[,names(grp_names)]

save(classMatrixQuery, file="example_classificationMatrix.rda"))
```

Plot classification heatmap:
```R
pdf(file="query_classification_heatmap.pdf"), height=4)

# This function can be found in plotting.R
acn_queryClassHm(classMatrixQuery, main = paste0("Classification Heatmap, ", study_name), 
            grps = grp_names, 
            fontsize_row=10, fontsize_col = 10, isBig = FALSE)
dev.off()
```

#### Compute GRN Status

Load `grnAll` and `trainNormParam` objects. If there is no relevant grnAll object available for your purposes, you will need to train a new one on AWS. See the [Appendix](#Appendix) for more. 
```R
grnAll <- utils_loadObject("Hs_grnAll_curatedTFs_Apr-22-2020.rda")
trainNormParam <- utils_loadObject("Hs_trainingNormalization_Apr-22-2020.rda")
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
---

## Appendix

### Reconstructing GRNs on AWS

GRN reconstruction requires a lot of RAM and thus makes more sense on an AWS instance. Once set up in the instance, GRN construction will take 15-60 min depending on size of expression matrix.

#### Starting an AWS instance

Start an EC2 c5.18xlarge instance with the following AMI:
Name: bigmomma_R3.5.2
AMI ID: ami-01473625b196bb951

Make sure you have 2 instance stores in addition to the root volume in the "Add Storage" tab. These allow you to mount the ephemeral0 and ephemeral1 drives. 

Select the "launch-wizard-1" security group.

Wait for 2/2 checks to be finished on the AWS GUI.

#### Upload necessary files to your instance

```
scp -i replace_with_key_name.pem Hs_stTrain_Jun-20-2017.rda ec2-user@<REPLACE-WITH-PUBLIC-DNS>.amazonaws.com:~
scp -i replace_with_key_name.pem Hs_expTrain_Jun-20-2017.rda ec2-user@<REPLACE-WITH-PUBLIC-DNS>.amazonaws.com:~
scp -i replace_with_key_name.pem cellnet_classifier_100topGenes_100genePairs.rda ec2-user@<REPLACE-WITH-PUBLIC-DNS>.amazonaws.com:~
scp -i replace_with_key_name.pem Hs_xpairs_list.rda ec2-user@<REPLACE-WITH-PUBLIC-DNS>.amazonaws.com:~
```

#### Begin 

ssh into your instance:
`ssh -i key_name.pem ec2-user@<REPLACE-WITH-PUBLIC-DNS>.amazonaws.com`

`screen` to preserve your session if disconnected.

Mount the storage drives and move .rda files:
```
$ sudo mkfs /dev/xvdb
$ sudo mount /dev/xvdb /media/ephemeral1
$ cd /media/ephemeral1
$ sudo mkdir analysis
$ sudo chown ec2-user analysis
$ cd analysis
$ mv ~/*.rda .
```

Start an R session:
`R`

Load dependencies:
```R
# library(devtools)
# install_github("pcahan1/cancerCellNet@v0.1.1", ref="master")
library(CellNet)
library(cancerCellNet)
library(igraph)
```

Load training data:
```R
expTrain <- utils_loadObject("Hs_expTrain_Jun-20-2017.rda")
expTrain_transformed <- trans_prop(weighted_down(expTrain, 5e5, dThresh=0.25), 1e5)

stTrain <- utils_loadObject("Hs_stTrain_Jun-20-2017.rda")

my_classifier <- utils_loadObject("cellnet_classifier_100topGenes_100genePairs.rda")
cnProc <- my_classifier$cnProc

xpairs_list <- utils_loadObject("Hs_xpairs_list.rda")
```

Construct GRNs (will take 15-60 min. Since you are in `screen`, you can leave your instance and come back anytime):
```R
system.time(grnAll <- ccn_makeGRN(expTrain_transformed, stTrain, "description1",
                                  zThresh = 4, dLevelGK = NULL, prune = TRUE,
                                  holm = 1e-4, cval=0.3)
            )
```

Save the grnAll object:
```R
save(grnAll, file="Hs_grnAll_todays_date.rda")
```

Train normalization parameters:
```R
expTrain_ranked <- logRank(expTrain2, base = 0)
# Extract the importance of genes based on the classifier
geneImportance <- processImportance(classifier = cnProc$classifier,
                                      xpairs = xpairs_list, prune = TRUE)
subNets <- grnAll$ctGRNs$geneLists
   

system.time(trainNormParam <- ccn_trainNorm(expTrain_ranked, stTrain,
                                            subNets=subNets,
                                            classList = geneImportance,
                                            dLevel = "description1", sidCol = "sra_id",
                                            classWeight = TRUE, exprWeight = FALSE,
                                            meanNorm = TRUE)
```

Save and exit:
```R
save(trainNormParam, file="Hs_trainingNormalization_todays_date.rda")
q(save="no")
```

To download the object, there are 2 options:
A. In a separate terminal tab on your local computer, scp the file from the instance to your local directory:
```
$ scp -i replace_with_key_name.pem ec2-user@<REPLACE-WITH-PUBLIC-DNS>.amazonaws.com:/media/ephemeral1/analysis/Hs_grnAll_todays_date.rda .
$ scp -i replace_with_key_name.pem ec2-user@<REPLACE-WITH-PUBLIC-DNS>.amazonaws.com:/media/ephemeral1/analysis/Hs_trainingNormalization_todays_date.rda .
```

OR

B. Use the AWS CLI to copy the object to S3, then download using the S3 GUI.
On your instance:
```
$ aws configure
AWS Access Key ID [None]: REPLACE_WITH_MY_ACCESS_KEY_ID
AWS Secret Access Key [None]: replace_with_my_secret_access_key
Default region name [None]: us-east-1e
Default output format [None]: 

$ aws s3 cp Hs_grnAll_todays_date.rda s3://cahanlab/my.folder/my-subfolder/
$ aws s3 cp Hs_trainingNormalization_todays_date.rda s3://cahanlab/my.folder/my-subfolder/
```

`exit` screen, then `exit` your instance.

TERMINATE YOUR INSTANCE ON THE AWS EC2 GUI. HOURLY RATES FOR LARGE INSTANCES ARE EXPENSIVE.




