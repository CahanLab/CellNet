# CellNet WebApp for top-pairs, rank-based agnostic analysis of expression profiles
suppressPackageStartupMessages({
   library(shiny)
   library(shinythemes)
   #library(shinycssloaders)
   library(shinyWidgets)
   library(shinyjs)
   library(CellNet)
   library(ggplot2)
   library(scater)
   library(cancerCellNet)
   library(RColorBrewer)
   library(pheatmap)
   library(igraph)
   source("plotting.R")
})

options(shiny.maxRequestSize = 50*1024^2)

ui <- fluidPage(
   theme = shinytheme("cosmo"),
   #theme = shinytheme("flatly"),
   
   useShinyjs(),
   
   titlePanel(windowTitle="CellNet Web",
              title = div("CellNet Web App", img(src="logo.png", height=60))),
   hr(),
   h2("Introduction"),
   p("CellNet is a network-biology-based, computational platform that assesses the fidelity of 
     cellular engineering and generates hypotheses for improving cell derivations. 
     CellNet is based on the reconstruction of cell type-specific gene regulatory networks (GRNs), 
     which we performed using publicly available RNA-Seq data of 16 mouse and 16 human cell and tissue types.
     Below, we describe how to apply CellNet to your RNA-Seq data."),
   p("See more at ", tags$a("https://github.com/pcahan1/CellNet")),
   hr(),
   
   h2("How to use this application"),
   h3("Inputs and outputs"),
   h4("Inputs: "),
   p("1) ", strong("Study name"), ": whatever you'd like to name your study"),
   p("2) ", strong("Sample metadata table"), "with the following columns at minimum: sample_name, description1"),
   p("3) ", strong("non-normalized"), "counts, TPM, or FPKM expression matrix with gene", strong("symbols"), "as row names and sample names as column names. 
     Sample names in counts matrix must match sample names in sample metadata table (but column order does not matter)."),

   h4("Outputs: "),
   p("1) ", strong("Classification heatmap"),": Columns represent query samples, and rows represent tissue types of the training data. 
     Each square is colored by the query sample's classification score for the given tissue type. 
     Scores range from 0 (distinct from the tissue type of the training data) to 1 (indistinguishable from the tissue type of the training data)."),
   p("2) ", strong("GRN Status"),": GRN status indicates the extent to which a tissue GRN is established in the training and query samples. 
     The raw GRN status is computed as the mean z-score of all genes in a tissue GRN, weighted by their importance to the associated tissue classifier. 
     The raw GRN status is then normalized to the mean raw GRN status of the training data samples of the given tissue. 
     Error bars represent mean ± 1 s.d."),
   p("3) ", strong("Transcription Factor Scores"),": The transcriptional regulators of the tissue GRN are shown on the y axis, 
     with the Network Influence Score (NIS) of each regulator on the x axis. 
     The NIS prioritizes transcription factors (TFs) such that their experimental perturbation is predicted to improve the target tissue classification. 
     The NIS of a TF is computed based on three components. 
     The first component is the extent to which the TF is dysregulated as compared with its expected value in the target tissue type. 
     The second component is the number of predicted targets of the TF. The third component is the extent to which the target genes are dysregulated."),
   
   
   #h3("Protocol"),
   #p("Step 1: "),
   hr(),
   
   
   sidebarLayout(
      sidebarPanel(
                     h4("Uploading your data"),
                     textInput(inputId="studyName", 
                               label="Study name"),
                     fileInput(inputId = "sampTabUploadCSV", label= "Upload .csv sample metadata table", 
                               accept = ".csv", multiple=FALSE),
                     fileInput(inputId = "expMatUpload", label = "Upload .csv counts matrix", 
                               accept = ".csv", multiple=FALSE),
                     selectInput(inputId="species", label="Species", choices = c("Human", "Mouse")),
                     uiOutput("tissueType"),
                     actionButton(inputId = "submit", label= "Submit")
                  ),
      
      mainPanel(
                  h3("CellNet Results"),
                  hr(),
                  #shinyjs for hidden, shinyWidgets for progressBar
                  hidden(div(id="progressDiv", progressBar(id = "progress", value = 0, total = 100, status = "info", 
                                                           display_pct = TRUE, striped = TRUE, title = "Progress..."))),
               
                  br(),
                  plotOutput("classHm", height="300px"),
                  br(),
                  textOutput("engineeredDescrip"),
                  br(),
                  plotOutput("engineeredRef", height="300px"),
                  br(),
                  plotOutput("GRNstatus", height="4000px"),
                  br(),
                  plotOutput("NIS", height="4000px"),
                  br(),
                  hidden(div(id="validat_plots", h4("Newly trained classifier performance"))),
                  textOutput("iGenes"),
                  plotOutput("newTrainHm"),
                  br(),
                  plotOutput("newTrainPR")
         )
   )
)

server <- function(input, output, session) {

   ########################################################################################
   ########################################################################################
   ##### FUNCTION DEFINITIONS #####
   
   trainClassifier <- function(expTrain, stTrain, querySampTab, queryExpDat, iGenes) {
      # Find intersecting genes between query sample and training samples
      iGenes <- Reduce(intersect, list(rownames(queryExpDat), iGenes))
      
      # Split expTrain into training and validation
      set.seed(39) # Set same seed every time for random number generator
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
      
      updateProgressBar(session = session, id = "progress", 
                        title = "Creating validation plots...",
                        value = 35, total = 100)
      
      # Call plotting function to make assessment plots
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
      
      output$iGenes <- renderText({
         paste("Number of intersecting genes: ", length(iGenes))
      })
      
      # Plot
      Sys.setlocale("LC_COLLATE", "C") # So that sort() is case sensitive
      output$newTrainHm <- renderPlot({
         ccn_hmClass(classMatrix_broad, grps=grps, fontsize_row=10)
      })
      
      # Classifier Assessment
      assessmentDat <- ccn_classAssess(classMatrix_broad, stValRand_broad, "description1","sra_id")
      
      output$newTrainPR <- renderPlot({
         plot_class_PRs(assessmentDat)
      }, height=720)
      
      return(broad_return)
   }
   
   queryClassifier <- function(broadReturn, studyName="Test", querySampTab, queryExpDat, groupBy="description3") {
      cnProc_broad <- broadReturn$cnProc
      
      # Query the classifier:
      classMatrixQuery <- broadClass_predict(cnProc = cnProc_broad, expDat = queryExpDat, nrand = 3)
      
      # Adjust labels
      grp_names <- c(as.character(querySampTab[,groupBy]), rep("random", 3))
      names(grp_names) <- c(as.character(querySampTab$sra_id), "rand_1", "rand_2", "rand_3")
      #names(grp_names) <- c(as.character(querySampTab$sample_name), "rand_1", "rand_2", "rand_3")
      # Re-order classMatrixQuery to match order of rows in querySampTab
      classMatrixQuery <- classMatrixQuery[,names(grp_names)]
      
      
      acn_queryClassHm(classMatrixQuery, main = paste0("Classification Heatmap, ", studyName),
                       grps = grp_names, is.Big = TRUE,
                       fontsize_row=9, fontsize_col = 10)
   }
   
   queryGRN <- function(grnAll, trainNormParam, broadReturn, queryExpDat_ranked, querySampTab) {
      updateProgressBar(session = session, id = "progress", title="Now analyzing GRN status...",
                        value = 70, total = 100)
      
      system.time(GRN_statusQuery <- ccn_queryGRNstatus(expQuery = queryExpDat_ranked, grn_return = grnAll, 
                                                        trainNorm = trainNormParam, classifier_return = broadReturn, prune = TRUE))
      
      cell_types <- rownames(GRN_statusQuery)
      GRN_statusQuery <- GRN_statusQuery[,querySampTab$sra_id]
      
      plot_list <- list()
      i <- 1
      for(type in cell_types) {
         plot_df <-  data.frame("SampleNames" = paste(colnames(GRN_statusQuery), querySampTab$description3),
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
         
         plot_list[[i]] <- type_plot
         i <- i+1
      }
      
      output$GRNstatus <- renderPlot({
         multiplot(plotlist=plot_list, cols=2)
      }, height=4000)
      

   }
   
   queryNIS <- function(grnAll, trainNormParam, broadReturn, queryExpDat_ranked, querySampTab, tissueType) {
      updateProgressBar(session = session, id = "progress", title="Now analyzing transcriptional regulators...",
                        value = 85, total = 100)
      system.time(TF_scores <- ccn_tfScores(expQuery = queryExpDat_ranked, grnAll = grnAll, trainNorm = trainNormParam,
                                            classifier_return = broadReturn, subnetName = tissueType,
                                            exprWeight = FALSE, normTFscore = TRUE))
      updateProgressBar(session = session, id = "progress", title="Finishing up...",
                        value = 97, total = 100)
      
      sample_names <- querySampTab$sra_id
      #sample_names <- querySampTab$sample_name
      
      plot_list <- list()
      i <- 1
      for(sample in sample_names) {
         descript <- querySampTab$description3[which(querySampTab$sra_id == sample)]
         #descript <- querySampTab$description3[which(querySampTab$sample_name == sample)]
         plot_df <- data.frame("TFs" = rownames(TF_scores),
                               "Scores" = as.vector(TF_scores[,sample]))
         sample_TFplot <- ggplot(plot_df, aes(x = TFs , y = Scores)) + geom_bar(stat="identity") + #aes(fill = medVal)) +
            theme_bw() + #scale_fill_gradient2(low = "purple", 
            # mid = "white", 
            # high = "orange") + 
            ggtitle(paste0(sample, ", ", descript, ", ", tissueType, " transcription factor scores")) +
            ylab("Network influence score") + xlab("Transcriptional regulator") + 
            theme(legend.position = "none", axis.text = element_text(size = 8)) +
            theme(text = element_text(size=10), 
                  legend.position="none",
                  axis.text.x = element_text(angle = 45, vjust=0.5))
         #coord_flip()
         plot_list[[i]] <- sample_TFplot
         i <- i+1
      }
      
      output$NIS <- renderPlot({
         multiplot(plotlist=plot_list, cols=1)
      }, height=4000)
   }
   
   
   
   ##### END OF FUNCTION DEFINITIONS #####
   ########################################################################################
   ########################################################################################
   
   
   
   ########################################################################################
   ########################################################################################
   ##### PROCESS INPUTS #####
   
   queryStudyName <- eventReactive(input$studyName, { 
      return(input$studyName)
   })
   
   # Read in file input for sample metadata table
   querySampTab <- eventReactive(input$sampTabUploadCSV, { 
      st <- read.csv(input$sampTabUploadCSV$datapath, stringsAsFactors = FALSE, header=TRUE)
      rownames(st) <- st$sra_id
      #rownames(st) <- st$sample_name
      return(st)
   })
   
   
   # Read in file input for expression matrix
   queryExpDat <- eventReactive(input$expMatUpload, { 
      tempDat <- read.csv(input$expMatUpload$datapath, stringsAsFactors = FALSE, 
                          header=TRUE, row.names=1) 
      return(tempDat)
   })
   
   species <- eventReactive(input$species, {
      input$species
   })
   
   
   output$tissueType <- renderUI({
      if (species() == "Human") {
         selectInput(inputId = "tissueType", label = "Target Cell/Tissue Type",
                  choices = c("lung","neuron","intestine_colon","kidney","heart","liver",
                              "skeletal_muscle","esc","endothelial_cell","hspc","b_cell",
                              "monocyte_macrophage","t_cell","fibroblast"))
         
      } else if (species() == "Mouse") {
         selectInput(inputId = "tissueType", label = "Target Cell/Tissue Type",
                     choices = c("t_cell","lung","nk_cell","hspc","heart",
                                 "macrophage","skeletal_muscle","neuron","wat",
                                 "kidney","intestine_colon","dendritic_cell",
                                 "b_cell","fibroblast","esc","liver"))
      }   
   })
   
   tissueType <- eventReactive(input$tissueType, { 
      return(input$tissueType)
   })
   
   ########################################################################################
   ########################################################################################
   
   
   
   
   ########################################################################################
   ########################################################################################
   ####### BEGIN ANALYSIS #######
   
   # Upon clicking submit:
   observeEvent(input$submit, {
      if (is.null(querySampTab()) || is.null(queryExpDat())) {
         shinyalert("Error!", "Please submit query sample table and expression matrix.", type = "error")
         return(NULL)
      }
      
      shinyjs::show("progressDiv") #shinyjs
      
      queryExpDat <- apply(queryExpDat(), 2, downSampleW, 1e5)
      queryExpDat <- log(1 + queryExpDat())
      
      updateProgressBar(session = session, id = "progress", title = "Loaded in query data...",
                        value = 10, total = 100)
      
      # Load in TISSUE-SPECIFIC broadReturn and iGenes
      ##### Fix naming scheme for classifiers and iGenes
      broadReturn <- utils_loadObject(paste0(tissueType(), "_broadClassifier100.rda"))
      iGenes <- utils_loadObject(paste0(tissueType(), "_studies_iGenes.rda"))
      ######
      
      # Check if retraining is necessary
      if (all(broadReturn$cnProc$cgenes %in% rownames(queryExpDat()))) {
         ## If no need to retrain:
         updateProgressBar(session = session, id = "progress", 
                           title = "No need to retrain! Starting classification...",
                           value = 40, total = 100)
      } else {
         ## If need to retrain:
         updateProgressBar(session = session, id = "progress", 
                           title = "Retraining classifier based on intersecting genes...This may take several minutes",
                           value = 20, total = 100)
         expTrain <- utils_loadObject("Hs_expTrain_Jun-20-2017.rda")
         stTrain <- utils_loadObject("Hs_stTrain_Jun-20-2017.rda")
         broadReturn <- trainClassifier(expTrain, stTrain, querySampTab(), queryExpDat(), iGenes)
         show("validat_plots")
         updateProgressBar(session = session, id = "progress", 
                           title = "Done retraining! Starting classification...",
                           value = 40, total = 100)
      }
         
      
      #Load tissueType-specific engineered reference data
      
      tissueTypeRefDat <- utils_loadObject(paste0(tissueType(), "_engineeredRef_expDat_all.rda"))
      tissueTypeRefSt <- utils_loadObject(paste0(tissueType(), "_engineeredRef_sampTab_all.rda"))
      
      
      ### Plot classification heatmap for user-inputted study
      output$classHm <- renderPlot({
         queryClassifier(broadReturn, queryStudyName(), querySampTab(), queryExpDat(), groupBy="description3")
      }, height=300)
      
      
      ### Plot classification heatmap for engineered reference studies
      output$engineeredDescrip <- renderText({
         "How does your protocol compare to publicly available data from related studies? \n \n"   
      })
      
      output$engineeredRef <- renderPlot({
         queryClassifier(broadReturn, studyName=paste0(input$tissueType, " engineered reference studies"), 
                         tissueTypeRefSt, tissueTypeRefDat, groupBy="study_id")
      }, height=300)
      
      updateProgressBar(session = session, id = "progress", title="Now starting GRN analysis...",
                        value = 55, total = 100)
      
      
      ###### GRN Analysis ######
      
      # Load in TISSUE-SPECIFIC grnAll and trainNormParam
      grnAll <- utils_loadObject(paste0(tissueType(), "_grnAll.rda"))
      trainNormParam <- utils_loadObject(paste0(tissueType(), "_trainNormParam.rda"))
      
      ### Check if subsetting grnAll and trainNormParam is necessary:
      updateProgressBar(session = session, id = "progress", title="Updating GRNs based on available genes...",
                        value = 60, total = 100)
      
      
      vertex_names <- V(grnAll$overallGRN$graph)$name
      
      if (!all(vertex_names %in% iGenes)) {
         allTargets <- grnAll$overallGRN$grnTable$TG
         newGRNTable <- grnAll$overallGRN$grnTable[which(allTargets %in% iGenes),]
         newTFsAll <- newGRNTable$TF
         newGRNTable <- newGRNTable[which(newTFsAll %in% iGenes),]
         grnAll$overallGRN$grnTable <- newGRNTable

         # Subset overallGRN graph based on iGenes
         vertex_names <- V(grnAll$overallGRN$graph)$name
         graph_iGenes <- which(vertex_names %in% iGenes)
         newGraph <- induced_subgraph(graph=grnAll$overallGRN$graph, vids=graph_iGenes, impl="copy_and_delete")
         grnAll$overallGRN$graph <- newGraph

         # Subset specGenes based on iGenes and tissue type
         tissueTypes <- names(grnAll$specGenes$context$general)
         newGeneral <- grnAll$specGenes$context$general
         for (tissue in tissueTypes) {
            tissueSpecGenes <- newGeneral[[tissue]]
            tissueSpecGenes <- tissueSpecGenes[which(names(tissueSpecGenes) %in% iGenes)]
            newGeneral[[tissue]] <- tissueSpecGenes
         }
         grnAll$specGenes$context$general <- newGeneral

         # Subset ctGRN geneLists, graphLists, and tfTargets  based on iGenes and tissue type
         grnAll$ctGRNs$geneLists <- newGeneral

         newGraphLists <- grnAll$ctGRNs$graphLists
         for (tissue in tissueTypes) {
            tissueGRN <- newGraphLists[[tissue]]
            iVertices <- vertex_attr(tissueGRN, name="name")
            iVertices <- iVertices[which(iVertices %in% iGenes)]
            tissueGRN <- induced_subgraph(graph=tissueGRN, vids=iVertices, impl="copy_and_delete")
            newGraphLists[[tissue]] <- tissueGRN
         }
         grnAll$ctGRNs$graphLists <- newGraphLists

         newTFTargets <- grnAll$ctGRNs$tfTargets
         for (tissue in tissueTypes) {
            tissueTFTargets <- newTFTargets[[tissue]]
            tissueTFTargets <- tissueTFTargets[which(names(tissueTFTargets) %in% iGenes)]
            for (TF in names(tissueTFTargets)) {
               newTargets <- tissueTFTargets[[TF]]
               newTargets <- newTargets[which(newTargets %in% iGenes)]
               tissueTFTargets[[TF]] <- newTargets
            }
            newTFTargets[[tissue]] <- tissueTFTargets
         }
         grnAll$ctGRNs$tfTargets <- newTFTargets

         # Subset trainNormParam
         newTVals <- trainNormParam$tVals
         for (tissue in tissueTypes) {
            newIndices <- which(names(newTVals[[tissue]][["mean"]]) %in% iGenes)
            newTVals[[tissue]][["mean"]] <- newTVals[[tissue]][["mean"]][newIndices]
            newTVals[[tissue]][["sd"]] <- newTVals[[tissue]][["sd"]][newIndices]
         }
         trainNormParam$tVals <- newTVals
      }
      
      updateProgressBar(session = session, id = "progress", title="Done updating GRNs!",
                        value = 65, total = 100)
      
      #####
      # GRN Status
      queryExpDat_ranked <- logRank(queryExpDat(), base = 0)
      
      queryGRN(grnAll, trainNormParam, broadReturn, queryExpDat_ranked, querySampTab())
      
      
      #####
      # NIS
      queryNIS(grnAll, trainNormParam, broadReturn, queryExpDat_ranked, querySampTab(), tissueType())
      

      updateProgressBar(session = session, id = "progress", title="Analysis Complete!",
                        value = 100, total = 100)

   })
   
}
shinyApp(ui = ui, server = server)
