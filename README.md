# CellNet

[Shortcut to Platform-Agnostic CellNet (PACNet) Web Application](http://cahanlab.org/resources/agnosticCellNet_web/)

[Shortcut to tutorial for running PACNet locally](https://github.com/pcahan1/CellNet/tree/master/agnostic_cellnet)

[Shortcut to bulk rna-seq protocol](#bulk_protocol)

[Cloud-based RNA-Seq web application](https://github.com/pcahan1/CellNet_Cloud)

[Microarray CellNet web application](http://cellnet.hms.harvard.edu/)

[Microarray CellNet code](https://pcahan1.github.io/cellnetr/)



### Introduction
CellNet is a network-biology-based, computational platform that assesses the fidelity of cellular engineering and generates hypotheses for improving cell derivations. CellNet is based on the reconstruction of cell type-specific gene regulatory networks (GRNs), which we performed using publicly available **RNA-Seq** data of 16 mouse and 16 human cell and tissue types. For now, there are two ways to run CellNet for RNA-Seq data. The easiest way to perform CellNet analysis is to use our web app. You can also run it as a command line tool on the cloud through Amazon Web Services, or you can run it locally. Below, we describe how to apply CellNet to your RNA-Seq data. 

For more, see our relevant publications:
[RNASeq CellNet Protocol](https://www.nature.com/articles/nprot.2017.022)
[Original CellNet Publication](https://www.sciencedirect.com/science/article/pii/S0092867414009349)

## Ways to Run CellNet

#### Web application
The web application takes as input an expression matrix (counts, TPM, or FPKM), and sample meta-data. The application performs CellNet analysis. Additionally, this tool includes analysis of many state-of-the-art differentiation protocols, so that you can benchmark your results against those commonly used methods:

[CellNet web app](https://cahanlab.org/resources/agnosticCellNet_web/)

#### Running CellNet in the Cloud 
The public CellNet Amazon Machine Image (AMI), available on Amazon Web Services (AWS), has all of the prerequisite software and libraries pre-installed. Because of this and the scalable computing capacity of AWS, **we highly recommend that you use AWS to run CellNet for RNA-Seq data** instead of running it locally. If you are unfamiliar with AWS or cloud computing in general, we recommend the following links for further information:
* [Cloud Computing concepts](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/concepts.html)
* [Amazon's Elastic Compute Cloud (EC2)](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/get-set-up-for-amazon-ec2.html)
* [Amazon Machine Images](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AMIs.html)
* [Amazon's Simple Storage Service (S3)](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AmazonS3.html)

 The current CellNet AMI (CellNet_v_0.1.1 ami-2ab59855, as of July 2018) is available in the AWS US East 1 region (N. Virginia region). Running CellNet on AWS requires uploading your raw data (in the form of .fastq files) either directly to your running instance on AWS EC2, or to S3 and then to your instance. To learn about transferring your fastq files directly to your instance, see [Transferring files to Linux machines using SCP](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AccessingInstancesLinux.html). **Note that Amazon charges by the hour for compute resources ($1.68/hour for a c3.8xlarge EC2 instance type)**. On average, it takes up to 2 hours to run a complete CellNet analysis for 144GB of raw data (9 samples of 16GB each).

#### Running CellNet Locally
Alternatively, you can run CellNet locally. The steps to do this are covered in our [Nature Protocol](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5765439/).

You will need to install the following command line software:
* [Cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html) 
* [Salmon](https://combine-lab.github.io/salmon/)
* [GNU Parallel](https://www.gnu.org/software/parallel/)
          
If you are using Mac OS, this can be done easily with PIP and Homebrew.


## Background information

#### Trained CellNet Objects (*cnProc*)
>At the heart of CellNet is the [Random Forest Classifier](https://en.wikipedia.org/wiki/Random_forest). This is the algorithm that will classify the results of a cell fate experiment. To analyze your own expression data with CellNet, you need a trained CellNet classifier object, which we refer to as a **cnProc** (CellNet Processor). You can select and use the appropriate cnProc that we have generated from the list below. You can also make your own using the code we provide [here](http://rdcu.be/rEmP). This is useful if you want to add more cell types, or if you want to train up a cnProc for a different species. **Note: generating a human cnProc requires a lot of computing power.** In general, it should be generated using an EC2 instance - it is probably not a good idea to try performing this locally.

The main ingredients of a cnProc are:
* An R matrix giving the expression levels of a number of genes across all the samples used to train CellNet
* An R dataframe providing metadata on the samples in the expression data matrix (things like cell-type, alignment metrics, SRA accession numbers...)

| SPECIES | DATE | CELL & TISSUE TYPES(# of profiles) | cnProc | raw training data |
|---------|------|------------------------------------|--------|-------------------|
| HS | Oct_25_2016 | b_cell (83), dendritic_cell (75), endothelial_cell (53), esc (52), fibroblast (79), heart (30), hspc (27), intestine_colon (64), kidney (29), liver (33), lung (95), macrophage (254), monocyte (207), neuron (109), skeletal_muscle (189), t_cell (53) | [Download](https://s3.amazonaws.com/cellnet-rnaseq/ref/cnproc/HS/cnProc_RS_hs_Oct_25_2016.rda) | |
| Mouse | Oct_24_2016 | b_cell (193), dendritic_cell (134), esc (134), fibroblast (182), heart (189), hspc (75), intestine_colon (149), kidney (109), liver (265), lung (116), macrophage (176), neuron (188), nk_cell (53), skeletal_muscle (130), t_cell (87), wat (64) | [Download](https://s3.amazonaws.com/cellnet-rnaseq/ref/cnproc/MM/cnProc_MM_RS_Oct_24_2016.rda) | [Download](https://s3.amazonaws.com/cellnet-rnaseq/ref/cnproc/MM/expTrain_MM_rawcounts_Dec_29_16.rda) | 
| Human | Apr_05_2017 | b_cell (83), dendritic_cell (55), endothelial_cell (51), esc (52), fibroblast (46), heart (60), hspc (192), intestine_colon (85), kidney (62), liver (107), lung (94), monocyte_macrophage (206), neuron (90), skeletal_muscle (187), t_cell (43) | [Download](https://s3.amazonaws.com/cellnet-rnaseq/ref/cnproc/HS/cnProc_HS_RS_Apr_05_2017.rda) | |
| Human | Jun_20_2017 | | [Download](https://s3.amazonaws.com/cellnet-rnaseq/ref/cnproc/HS/cnProc_HS_RS_Jun_20_2017.rda) |  [Download](https://s3.amazonaws.com/cellnet-rnaseq/ref/cnproc/HS/expTrain_HS_rawcounts_8_31_2017.rda) |


#### Example Data

These are some datasets you can use to test-drive applying CellNet to RNA-Seq data:

| SPECIES | DATE | SRA ID | DESCRIPTION | METADATA | EXPRESSION |
|---------|------|--------|-------------|----------|------------|
| Human   | Oct 30, 2015 | SRP043684 | Engineered Neurons | [metadata](https://s3.amazonaws.com/cellnet-rnaseq/ref/examples/st_SRP043684_example.rda) | [expression data](https://s3.amazonaws.com/cellnet-rnaseq/ref/examples/expList_SRP043684_example.rda) |
| Mouse | Mar 15, 2016 | SRP059670 | Reprogramming to Pluripotency | [metadata](https://s3.amazonaws.com/cellnet-rnaseq/ref/examples/st_SRP059670_example.rda) | [expression data](https://s3.amazonaws.com/cellnet-rnaseq/ref/examples/expList_SRP059670_example.rda) |

#### Salmon Index Table

If you are running CellNet locally, you will need to have salmon installed on your machine. Below are a few indexes that we have created from our transcriptome and know to work.

| SPECIES | SALMON | INDEX DOWNLOAD | NOTE/USAGE |
|---------|----------------|----------------|------------|
| Human | 0.6.0 | [salmon.index.human.050316.tgz](https://s3.amazonaws.com/cellnet-rnaseq/ref/salmon.index.human.050316.tgz) | Default for AWS workflow |
| Mouse | 0.6.0 | [salmon.index.mouse.050316.tgz](https://s3.amazonaws.com/cellnet-rnaseq/ref/salmon.index.mouse.050316.tgz) | Default for AWS workflow |
| Human | 0.7.3 | [salmon.index.human.122116.tgz](https://s3.amazonaws.com/cellnet-rnaseq/ref/salmon.index.human.122116.tgz) | Protocol for local** |
| Mouse | 0.7.3 | [salmon.index.mouse.122116.tgz](https://s3.amazonaws.com/cellnet-rnaseq/ref/salmon.index.mouse.122116.tgz) | Protocol for local** |
| Human | 0.8.2 | [salmon.index.human.052617.tgz](https://s3.amazonaws.com/cellnet-rnaseq/ref/salmon.index.human.052617.tgz) | Uses latest version of Salmon to date |
| Mouse | 0.8.2 | [salmon.index.mouse.052617.tgz](https://s3.amazonaws.com/cellnet-rnaseq/ref/salmon.index.mouse.052617.tgz) | Uses latest version of Salmon to date |

** Here's the [binary Salmon-0.7.3 Mac OSX link](https://github.com/COMBINE-lab/salmon/files/581546/Salmon-0.7.3-pre_OSX_10.11.tar.gz). Salmon-0.8.2 is a stable update and will work for either MacOSX or Linux.

#### <a name="bulk_protocol">Running CellNet on AWS</a>

The steps below demonstrate how to run RNA-Seq CellNet on AWS. You need to log in to the [AWS console](https://console.aws.amazon.com), then click on EC2, and launch the CellNet_v_0.1.1 image (ami-2ab59855) on a c3.4×large or c3.8×large instance type. 

To log in to the running instance, type the following command in the shell/terminal, but _aws_private_key_ with the full path of the AWS key that you used to launch the instance. And replace _instance_public_dns_ with the public DNS of your instance that can be found in the AWS console

```Bash
ssh -i aws_private_key ec2-user@instance_public_dns
```

Once, you have logged in to the instance, you should launch screen. 


 Then, you need to load the latest version of CellNet (v0.1.1), which is pre-installed on this image. (However, you can follow [these steps to install CellNet](#install_cellnet), if you ever need to do so.) You also need to configure the instance for pre-processing RNA-Seq data, including fetching the mouse transcriptome index. It is important to work in the same R session as where you call _cn_setup()_.
```R
screen
R
library(CellNet)
cn_setup() 
```

Now fetch the demonstration data, and the mouse cnProc
```R
fetch_salmon_indices(species="mouse")
download.file("https://s3.amazonaws.com/CellNet/rna_seq/mouse/examples/SRP059670/st_SRP059670_example.rda", "st_SRP059670_example.rda")
stQuery <- utils_loadObject("st_SRP059670_example.rda")
stQuery <- cn_s3_fetchFastq("CellNet","rna_seq/mouse/examples/SRP059670",stQuery,fname="fname", compressed="gz")

download.file("https://s3.amazonaws.com/CellNet/rna_seq/mouse/cnProc_MM_RS_Oct_24_2016.rda", dest="./cnProc_MM_RS_Oct_24_2016.rda")
cnProc<-utils_loadObject("cnProc_MM_RS_Oct_24_2016.rda")
```

Pre-process the fastq files. This runs Salmon to estimate transcript abundances:
```R
expList <- cn_salmon(stQuery) ## Assumes your fastq files are in the working directory. This takes ~15 minutes on the demo data
```
    
Applying CellNet:
```R
cnRes <- cn_apply(expList[['normalized']], stQuery, cnProc)
```

#### Interpreting Output
CellNet produces a number of outputs, the most commonly used is the cnRes Object (CellNet Result). There are three figures that can be created from this:

Classification Heat Map: Displays classification score of each sample (column) to each of the cell and tissue types in the training data (rows):
```R
pdf(file='hmclass_example.pdf', width=7, height=5)
cn_HmClass(cnRes)
dev.off()
```

![](md_img/hm.png)

To fetch this and other files that you save on AWS/EC2, you can use _scp_ as shown below, replacing _aws_private_key_ and _instance_public_dns_ with your values:
```Bash
scp -i aws_private_key ec2-user@instance_public_dns:/media/ephemeral0/analysis/*.pdf ./

```

**G**ene **R**egulatory **N**etwork Status Bar Plot: A more sensitive measure of the degree to which a particular cell type's GRN has been established in your experimental data.
```R
fname<-'grnstats_esc_example.pdf'
bOrder<-c("fibroblast_train", unique(as.vector(stQuery$description1)), "esc_train")
cn_barplot_grnSing(cnRes,cnProc,"esc", c("fibroblast","esc"), bOrder, sidCol="sra_id")
ggplot2::ggsave(fname, width=5.5, height=5)
dev.off()
```
<img src="md_img/grnStat.png" style="height: 400px;"/>

**N**etwork **I**nfluence **S**core Box and Whisker Plot: A suggestion of transcription factors that could be better regulated, ranked by their potential impact
```R
rownames(stQuery)<-as.vector(stQuery$sra_id)
tfScores<-cn_nis_all(cnRes, cnProc, "esc") 

fname<-'nis_esc_example_Day0.pdf'
plot_nis(tfScores, "esc", stQuery, "Day0", dLevel="description1", limitTo=0) 
ggplot2::ggsave(fname, width=4, height=12)
dev.off()
```
<img src="md_img/nis.png" style="height: 400px;"/>



#### <a name="install_cellnet">Installing CellNet</a>

If, for some reason, you need to install CellNet anew, you can do so by using devtools:
```R
sudo R
library(devtools)
install_github("pcahan1/CellNet", ref="master")
q(save='no')
```


