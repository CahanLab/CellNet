
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


After you are finished downloading the final objects, `exit` screen, then `exit` your instance.

TERMINATE YOUR INSTANCE ON THE AWS EC2 GUI. HOURLY RATES FOR LARGE INSTANCES ARE EXPENSIVE.




