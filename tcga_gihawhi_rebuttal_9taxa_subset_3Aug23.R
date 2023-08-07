# tcga_gihawi_rebuttal_RT_subset_31July23.R
# Author: Greg Poore
# Date: Aug 3, 2023
# Purpose: To explore cancer type differences in data released by 
# Gihawi et al. 2023 bioRxiv: https://www.biorxiv.org/content/10.1101/2023.07.28.550993v1
# only using bacterial and fungal genera found in tumors by the Weizmann Institute of Science (RT)

# Load dependencies
require(ggplot2) # 3.4.2
require(ggpubr) # 0.4.0
require(ggsci) # 2.9
require(dplyr) # 1.1.1
require(tibble) # 3.2.1
require(doMC) # 1.3.7
require(gtools) # 3.9.2
require(biomformat) # 1.22.0
require(Rhdf5lib) # 1.16.0
require(EnvStats) # 2.4.0
require(ggbeeswarm) # 0.7.1
require(phyloseq) # 1.38.0
require(reshape2) # 1.4.4
require(forcats) # 0.5.1
require(rstatix) # 0.7.0
require(pROC) # 1.18.0

numCores <- detectCores()
registerDoMC(cores=numCores)

#------------------------------------------------------#
# Input data
#------------------------------------------------------#

metadataSamplesAll <- read.csv("tcga_metadata_poore_et_al_2020_1Aug23.csv", row.names = 1)

# Load tables from Gihawi et al. listed here:
# https://github.com/yge15/Cancer_Microbiome_Reanalyzed
blcaKraken <- read.csv("TableS8_BLCA.all.csv",
                       row.names = 1, stringsAsFactors = FALSE)
hnscKraken <- read.csv("TableS9_HNSC_all.csv",
                       row.names = 1, stringsAsFactors = FALSE)
colnames(hnscKraken) <- gsub("^g_","",colnames(hnscKraken)) # Uniformize genera names
brcaKraken <- read.csv("TableS10_BRCA_WGS.csv",
                       row.names = 1, stringsAsFactors = FALSE)
colnames(brcaKraken) <- gsub("^g_","",colnames(brcaKraken))  # Uniformize genera names

#------------------------------------------------------#
# Select 9 "well-known" common genera that were also in the RT data
#------------------------------------------------------#
# 9 taxa
taxa2Keep <- c("Prevotella","Staphylococcus","Fusobacterium",
               "Streptococcus","Veillonella","Pseudomonas",
               "Klebsiella","Acinetobacter","Clostridium")

## Subset tables --> all 9 taxa are found in each cancer type
sum(colnames(blcaKraken) %in% taxa2Keep) # 9
sum(colnames(hnscKraken) %in% taxa2Keep) # 9
sum(colnames(brcaKraken) %in% taxa2Keep) # 9

# "RT" = "restricted taxa"
blcaKrakenRT <- blcaKraken[,colnames(blcaKraken) %in% taxa2Keep]
hnscKrakenRT <- hnscKraken[,colnames(hnscKraken) %in% taxa2Keep]
brcaKrakenRT <- brcaKraken[,colnames(brcaKraken) %in% taxa2Keep]

#------------------------------------------------------#
# Merge tables
# - Selected conservative method to retain only overlapping features
# (although that reduce some effect size since the files are
# separated by cancer type)
# - Note that RT did not contain "Homo" genus calls, so they are inherently excluded
#------------------------------------------------------#
sharedFeat <- Reduce(intersect, list(colnames(blcaKrakenRT),
                                     colnames(hnscKrakenRT),
                                     colnames(brcaKrakenRT)))

countMergedRT <- smartbind(blcaKrakenRT[,sharedFeat],
                            hnscKrakenRT[,sharedFeat],
                            brcaKrakenRT[,sharedFeat])
# Missing values after the merge should be converted to 0s
countMergedRT[is.na(countMergedRT)] <- 0 
rownames(countMergedRT) <- c(rownames(blcaKrakenRT),
                           rownames(hnscKrakenRT),
                           rownames(brcaKrakenRT))
dim(countMergedRT) # 728 9

# Subset metadata to samples reflected by Gihawi et al. 2023
metadataSamplesMergedRT <- metadataSamplesAll[rownames(countMergedRT),]

# Subset to primary tumor ("PT") and blod derived normal ("BDN") samples
metadataSamplesMergedRT_PT <- metadataSamplesMergedRT %>% 
  filter(sample_type == "Primary Tumor") %>% droplevels()
metadataSamplesMergedRT_BDN <- metadataSamplesMergedRT %>% 
  filter(sample_type == "Blood Derived Normal") %>% droplevels()

#------------------------------------------------------#
# Examine PT sequencing center subsets
#
# NOTE: The goal is to have individual seq-center subsets
# to avoid needing to batch correction. This allows us to use
# the raw data for comparisons.
#------------------------------------------------------#

# Count PT subsets
metadataSamplesMergedRT_PT %>% count(data_submitting_center_label, investigation)
# data_submitting_center_label investigation  n
# 1                          Baylor College of Medicine     TCGA-HNSC 32
# 2                  Broad Institute of MIT and Harvard     TCGA-BLCA 16
# 3                  Broad Institute of MIT and Harvard     TCGA-HNSC 31
# 4                              Harvard Medical School     TCGA-BLCA 51
# 5                              Harvard Medical School     TCGA-BRCA 19
# 6                              Harvard Medical School     TCGA-HNSC 85
# 7  MD Anderson - Institute for Applied Cancer Science     TCGA-BLCA 62
# 8  MD Anderson - Institute for Applied Cancer Science     TCGA-BRCA  2
# 9  MD Anderson - Institute for Applied Cancer Science     TCGA-HNSC 21
# 10      MD Anderson - RPPA Core Facility (Proteomics)     TCGA-HNSC  1
# 11           Washington University School of Medicine     TCGA-BRCA 93

##--> Harvard Med School (HMS) has 3 cancer types and the most samples by far
##--> Broad has 2 cancer types
##--> MD Anderson (MDA) effectively has 2 cancer types (after removing the 2 breast cancer samples)
##--> Everywhere else cannot be used without batch correction, so they will be excluded

metadataSamplesMergedRT_PT_HMS <- metadataSamplesMergedRT_PT %>%
  filter(data_submitting_center_label == "Harvard Medical School") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  droplevels()
countMergedRT_PT_HMS <- countMergedRT[metadataSamplesMergedRT_PT_HMS$sampleid,]

metadataSamplesMergedRT_PT_Broad <- metadataSamplesMergedRT_PT %>%
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  droplevels()
countMergedRT_PT_Broad <- countMergedRT[metadataSamplesMergedRT_PT_Broad$sampleid,]

metadataSamplesMergedRT_PT_MDA <- metadataSamplesMergedRT_PT %>%
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  filter(investigation != "TCGA-BRCA") %>% # Drop the 2 samples
  droplevels()
countMergedRT_PT_MDA <- countMergedRT[metadataSamplesMergedRT_PT_MDA$sampleid,]

#------------------------------------------------------#
# Examine BDN sequencing center subsets
#
# NOTE: The goal is to have individual seq-center subsets
# to avoid needing batch correction. This allows us to use
# the raw data for comparisons.
#------------------------------------------------------#
metadataSamplesMergedRT_BDN %>% count(data_submitting_center_label, investigation)
# data_submitting_center_label investigation  n
# 1                         Baylor College of Medicine     TCGA-HNSC 28
# 2                 Broad Institute of MIT and Harvard     TCGA-HNSC 14
# 3                             Harvard Medical School     TCGA-BRCA 19
# 4                             Harvard Medical School     TCGA-HNSC 76
# 5 MD Anderson - Institute for Applied Cancer Science     TCGA-HNSC 22
# 6           Washington University School of Medicine     TCGA-BRCA 87

##--> Harvard Med School (HMS) has 2 cancer types
##--> Everywhere else cannot be used without batch correction, so they will be excluded

metadataSamplesMergedRT_BDN_HMS <- metadataSamplesMergedRT_BDN %>%
  filter(data_submitting_center_label == "Harvard Medical School") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  droplevels()
countMergedRT_BDN_HMS <- countMergedRT[metadataSamplesMergedRT_BDN_HMS$sampleid,]

#------------------------------------------------------#
# For completeness - sequencing platform
#
# - NOTE: 1 sample was processed on a MiSeq, but it was at WashU,
# which is not included in the PT or BDN subsets (so ignore)
#------------------------------------------------------#
metadataSamplesMergedRT %>% count(platform, data_submitting_center_label)
# platform                       data_submitting_center_label   n
# 1 Illumina HiSeq                         Baylor College of Medicine  65
# 2 Illumina HiSeq                 Broad Institute of MIT and Harvard  67
# 3 Illumina HiSeq                             Harvard Medical School 284
# 4 Illumina HiSeq MD Anderson - Institute for Applied Cancer Science 113
# 5 Illumina HiSeq      MD Anderson - RPPA Core Facility (Proteomics)   1
# 6 Illumina HiSeq           Washington University School of Medicine 197
# 7 Illumina MiSeq           Washington University School of Medicine   1


#------------------------------------------------------#
# Examine prevalence of these 9 taxa
#------------------------------------------------------#

# Among HMS PTs
prevRT_PT_HMS <- round(colSums(countMergedRT_PT_HMS>0) / 
                                     nrow(countMergedRT_PT_HMS)*100,2)

# Among HMS BDNs
prevRT_BDN_HMS <- round(colSums(countMergedRT_BDN_HMS>0) / 
                              nrow(countMergedRT_BDN_HMS)*100,2)

# Combine the data
prevRT_PT_BDN_HMS <- cbind(as.data.frame(prevRT_PT_HMS),
                           as.data.frame(prevRT_BDN_HMS))

# Calculate average prevalence for these restricted taxa
colMeans(prevRT_PT_BDN_HMS)
# prevRT_PT_HMS prevRT_BDN_HMS 
# 81.28778       72.63111

## Plot correlation with labels
prevRT_PT_BDN_HMS %>%
  rownames_to_column("Microbe") %>%
ggscatter(x = "prevRT_BDN_HMS",
          y = "prevRT_PT_HMS",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Prevalence among HMS blood (%)", 
          ylab = "Prevalence among HMS tumors (%)") +
  theme(aspect.ratio=1) +
  ggrepel::geom_label_repel(aes(label = Microbe), size = 2.5, force = 5)
ggsave(filename = "Figures/prevRT_HMS_PTvsBDN_corr_plot_5Aug23.jpeg",
       dpi = "retina", units = "in", width = 4, height = 4)

## Barplot with labels
prevRT_PT_BDN_HMS %>%
  arrange(desc(prevRT_PT_HMS)) %>%
  rownames_to_column("Microbe") %>%
  melt() %>%
  mutate(variable = gsub("prevRT_PT_HMS","HMS PT",variable),
         variable = gsub("prevRT_BDN_HMS","HMS BDN",variable)) %>%
  rename(`Sample type`=variable) %>%
  ggbarplot(x = "Microbe",
            y = "value",
            fill = "Sample type",
            palette = "nejm",
            label = TRUE,
            lab.size = 2.7,
            xlab = "",
            ylab = "Prevalence",
            position = position_dodge(0.9)) +
  rotate_x_text(30)
ggsave(filename = "Figures/prevRT_HMS_PTvsBDN_barplot_5Aug23.jpeg",
       dpi = "retina", units = "in", width = 6, height = 5)

#------------------------------------------------------#
# Calculate log-ratios of cancer types using these restricted taxa
#------------------------------------------------------#

lrFUN <- function(countTable,metaData,
                  taxaNumerator = c("Prevotella","Fusobacterium",
                                    "Streptococcus","Veillonella",
                                    "Klebsiella","Clostridium","Pseudomonas"),
                  taxaDenominator = c("Staphylococcus","Acinetobacter"),
                  yStatPos = c(10, 13, 16),
                  plotPalette = c("#0072B5FF","#BC3C29FF","#E18727FF")){
  
  print(sprintf("Matched sample order: %s",all(metaData$sampleid == rownames(countTable))))
  metaDataLRs <- metaData %>%
    mutate(LR = log2(rowSums(countTable[,taxaNumerator]) / 
                       rowSums(countTable[,taxaDenominator]))) %>%
    filter(is.finite(LR))
  
  print(sprintf("Number samples with Inf LR and dropped: %d",
                nrow(metaData)-nrow(metaDataLRs)))
  
  metaDataLRs_Formatted <- metaDataLRs %>%
    mutate(investigation = gsub("^TCGA-","",investigation)) %>%
    select(investigation, LR) %>%
    mutate(investigation = fct_reorder(investigation, LR, median))
    
  pairwise.test <- metaDataLRs_Formatted %>% 
    wilcox_test(LR ~ investigation)
  print(pairwise.test)
  
  metaDataLRs_Formatted %>%
    ggboxplot(x = "investigation",
              y = "LR",
              xlab = "",
              ylab = "Log-ratio",
              fill = "investigation",
              palette = plotPalette,
              legend = "none") +
    stat_pvalue_manual(
      pairwise.test, label = ifelse(length(levels(metaDataLRs_Formatted$investigation))>2,
                                    yes = "q={p.adj}", no = "p={p}"),
      y.position = yStatPos
    ) +
    stat_n_text() -> lrPlot
  print(lrPlot)
  
  res <- list(lrPlot=lrPlot,
              metaDataLRs_Formatted=metaDataLRs_Formatted)
}

lrRT_PT_HMS <- lrFUN(countTable = countMergedRT_PT_HMS,
                     metaData = metadataSamplesMergedRT_PT_HMS,
                     taxaNumerator = c("Prevotella","Fusobacterium",
                                       "Streptococcus","Veillonella",
                                       "Klebsiella","Clostridium","Pseudomonas"),
                     taxaDenominator = c("Staphylococcus","Acinetobacter"),
                     plotPalette = c("#BC3C29FF","#0072B5FF","#E18727FF"))
ggsave(filename = "Figures/lrRT_HMS_PT_boxplot_5Aug23.jpeg",
       dpi = "retina", units = "in", width = 3, height = 4)

lrRT_BDN_HMS <- lrFUN(countTable = countMergedRT_BDN_HMS,
                      metaData = metadataSamplesMergedRT_BDN_HMS,
                      yStatPos = c(5),
                      taxaNumerator = c("Prevotella","Fusobacterium",
                                        "Streptococcus","Veillonella",
                                        "Klebsiella","Clostridium"),
                      taxaDenominator = c("Staphylococcus","Acinetobacter"),
                      plotPalette = c("#0072B5FF","#E18727FF"))
ggsave(filename = "Figures/lrRT_HMS_BDN_boxplot_5Aug23.jpeg",
       dpi = "retina", units = "in", width = 2.5, height = 4)

lrRT_PT_MDA <- lrFUN(countTable = countMergedRT_PT_MDA,
                     metaData = metadataSamplesMergedRT_PT_MDA,
                     taxaNumerator = c("Prevotella","Fusobacterium",
                                       "Streptococcus","Veillonella",
                                       "Klebsiella","Clostridium","Pseudomonas"),
                     taxaDenominator = c("Staphylococcus","Acinetobacter"),
                     yStatPos = c(12),
                     plotPalette = c("#BC3C29FF","#E18727FF"))
ggsave(filename = "Figures/lrRT_MDA_PT_boxplot_5Aug23.jpeg",
       dpi = "retina", units = "in", width = 2.5, height = 4)

lrRT_PT_Broad <- lrFUN(countTable = countMergedRT_PT_Broad,
                     metaData = metadataSamplesMergedRT_PT_Broad,
                     taxaNumerator = c("Prevotella","Fusobacterium",
                                       "Streptococcus","Veillonella",
                                       "Klebsiella","Clostridium","Pseudomonas"),
                     taxaDenominator = c("Staphylococcus","Acinetobacter"),
                     yStatPos = c(12),
                     plotPalette = c("#BC3C29FF","#E18727FF"))
ggsave(filename = "Figures/lrRT_Broad_PT_boxplot_5Aug23.jpeg",
       dpi = "retina", units = "in", width = 2.5, height = 4)
  
#----------------------------------------------------------#
# HMS PT multi-class machine learning
#
# NOTE: Since HMS is the only seq-center subset that had 3 cancer types,
# it will be the only multi-class ML instance
#----------------------------------------------------------#
require(caret) # for model building (6.0-90)
require(gbm) # for machine learning (2.1.8)
require(xgboost) # for machine learning (1.5.0.1)
require(randomForest) # for machine learning (4.6-14)
require(PRROC) # for precision-recall curves (1.3.1)
require(MLmetrics) # for multi-class learning (1.1.1)

## Write ML function
mlMulticlass <- function(metaData = metadataSamplesMergedRT_PT_HMS,
                         countData = countMergedRT_PT_HMS,
                         sampleType = "Primary Tumor",
                         seqCenter = "Harvard Medical School",
                         modelType = "xgbTree",
                         numResampleIter = 1,
                         numKFold = 10){
  
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- gsub('([[:punct:]])|\\s+','',seqCenter)
  numSeed <- 42
  xgbGrid <- data.frame(nrounds = 10,
                        max_depth = 4,
                        eta = .1,
                        gamma = 0,
                        colsample_bytree = .7,
                        min_child_weight = 1,
                        subsample = .8)
  
  metaDataFilt <- metaData %>% column_to_rownames("sampleid")
  
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  countData <- countData[rownames(metaDataFilt),] # re-align in case
  
  mlDataY <- metaDataFilt
  mlDataX <- countData
  
  # Use cross-validation:
  trainX <- mlDataX
  trainY <- mlDataY[,"predY"]
  print(table(trainY))
  
  set.seed(numSeed) # have to restate seed again, 
  # as ctrl defines the cross validation sampling during training
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       sampling = "up",
                       summaryFunction = multiClassSummary,
                       classProbs = TRUE,
                       verboseIter = FALSE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  print("Now training model...")
  set.seed(numSeed)
  mlModel <- train(x = trainX,
                   y = trainY,
                   method = modelType,
                   # preProcess = c("zv"), # remove zero variance features
                   nthread = 1,
                   trControl = trainControl(method = "repeatedcv",
                                            number = numKFold,
                                            repeats = numResampleIter,
                                            sampling = "up",
                                            summaryFunction = multiClassSummary,
                                            classProbs = TRUE,
                                            verboseIter = FALSE,
                                            savePredictions = TRUE,
                                            allowParallel=TRUE),
                   # metric = "ROC",
                   tuneGrid = xgbGrid,
                   verbose = FALSE)
  
  resPredFun <- function(mlModel){
    resPred <- mlModel$pred
    
    ## Split folds and calculate perf on each fold
    resPredSplit <- split(resPred, resPred$Resample)
    repX_perf <- list()
    for(zz in seq_along(resPredSplit)){
      resPredSingleRep <- resPredSplit[[zz]]
      predProbs <- resPredSingleRep
      multiClass <- resPredSingleRep
      rep_perfTmp <- multiClassSummary(multiClass, lev = levels(multiClass$obs))
      repX_perf[[zz]] <- data.frame(as.list(rep_perfTmp))
    }
    
    # SUMMARIZE MODEL PERFORMANCES
    rep_perf <- do.call(rbind, repX_perf)
    # print(rep_perf)
    return(rep_perf)
  }
  
  print("Obtaining performance values...")
  resPredAll <- mlModel$pred
  perfCombinedAll <-  multiClassSummary(resPredAll, lev = levels(resPredAll$obs))
  repPerfCombinedAll <- resPredFun(mlModel)
  
  # Print AUC
  print(perfCombinedAll[2])
  
  ## Save predictions and perf
  # Preds
  baseFilename <- paste0("multiclassCV_RT_",st,"_",sc,"_k",numKFold,"_modelType_",modelType)
  write.csv(resPredAll, file = paste0("ML_results/pred",baseFilename,".csv"))
  save(resPredAll, file = paste0("ML_results/pred",baseFilename,".RData"))
  # Overall perf
  write.csv(perfCombinedAll, file = paste0("ML_results/perfCombinedAll",baseFilename,".csv"))
  save(perfCombinedAll, file = paste0("ML_results/perfCombinedAll",baseFilename,".RData"))
  # Rep perf
  write.csv(resPredAll, file = paste0("ML_results/repPerfCombinedAll",baseFilename,".csv"))
  save(resPredAll, file = paste0("ML_results/repPerfCombinedAll",baseFilename,".RData"))
  
  ## Variable importances
  xgb_imp <- as.data.frame(xgb.importance(feature_names = mlModel$finalModel$feature_names,
                                          model = mlModel$finalModel))
  write.csv(xgb_imp, file = paste0("ML_results/featureImp",baseFilename,".csv"))
  # Print top 5 features
  print("Top feature importances:")
  print(head(xgb_imp, 50))
  
  res <- list(resPredAll=resPredAll,
              mlModel=mlModel,
              varImp=xgb_imp,
              perfCombinedAll=perfCombinedAll,
              repPerfCombinedAll=repPerfCombinedAll)
  return(res)
  
  rm(mlModel)
}

# Call multiclass function
hmsMulticlassMLRes_PT <- mlMulticlass()

# Examine overall performance:
hmsMulticlassMLRes_PT$perfCombinedAll
# logLoss                    AUC                  prAUC 
# 0.5636223              0.9634837              0.8342653 
# Accuracy                  Kappa                Mean_F1 
# 0.8838710              0.8030496              0.8508288 
# Mean_Sensitivity       Mean_Specificity    Mean_Pos_Pred_Value 
# 0.8786378              0.9436705              0.8333686 
# Mean_Neg_Pred_Value         Mean_Precision            Mean_Recall 
# 0.9350384              0.8333686              0.8786378 
# Mean_Detection_Rate Mean_Balanced_Accuracy 
# 0.2946237              0.9111542 

# Examine per-CV performance:
hmsMulticlassMLRes_PT$repPerfCombinedAll

# Create confusion matrix
caret::confusionMatrix(hmsMulticlassMLRes_PT$resPredAll$pred, 
                       hmsMulticlassMLRes_PT$resPredAll$obs)

# Confusion Matrix and Statistics
# 
# Reference
# Prediction BLCA BRCA HNSC
# BLCA   42    2    5
# BRCA    6   17    2
# HNSC    3    0   78
# 
# Overall Statistics
# 
# Accuracy : 0.8839          
# 95% CI : (0.8227, 0.9297)
# No Information Rate : 0.5484          
# P-Value [Acc > NIR] : <2e-16          
# 
# Kappa : 0.803           
# 
# Mcnemar's Test P-Value : 0.2123          
# 
# Statistics by Class:
# 
#                      Class: BLCA Class: BRCA Class: HNSC
# Sensitivity               0.8235      0.8947      0.9176
# Specificity               0.9327      0.9412      0.9571
# Pos Pred Value            0.8571      0.6800      0.9630
# Neg Pred Value            0.9151      0.9846      0.9054
# Prevalence                0.3290      0.1226      0.5484
# Detection Rate            0.2710      0.1097      0.5032
# Detection Prevalence      0.3161      0.1613      0.5226
# Balanced Accuracy         0.8781      0.9180      0.9374

# #----------------------------------------------------------#
# # HMS PT breakdown 2 class ML (for later plotting)
# #----------------------------------------------------------#
# 
# hmsMulticlassMLRes_PT_ <- mlMulticlass()



#----------------------------------------------------------#
# BDN 2-class machine learning
#
# - Since HMS is the only seq-center subset that had >1 cancer types,
# it will be the only multi-class ML instance.
# - Although not shown here, this could also be applied to the 
# PT 2-class subsets above
#----------------------------------------------------------#

# Call multiclass function
hmsMulticlassMLRes_BDN <- mlMulticlass(metaData = metadataSamplesMergedRT_BDN_HMS,
                                   countData = countMergedRT_BDN_HMS,
                                   sampleType = "Blood Derived Normal",
                                   seqCenter = "Harvard Medical School",
                                   modelType = "xgbTree",
                                   numResampleIter = 1,
                                   numKFold = 10)

# Examine overall performance:
hmsMulticlassMLRes_BDN$perfCombinedAll
# logLoss               AUC             prAUC          Accuracy 
# 0.4221401         0.9002770         0.8374327         0.8947368 
# Kappa                F1       Sensitivity       Specificity 
# 0.6835443         0.7500000         0.7894737         0.9210526 
# Pos_Pred_Value    Neg_Pred_Value         Precision            Recall 
# 0.7142857         0.9459459         0.7142857         0.7894737 
# Detection_Rate Balanced_Accuracy 
# 0.1578947         0.8552632 

# Examine per-CV performance:
hmsMulticlassMLRes_BDN$repPerfCombinedAll

# Create confusion matrix
caret::confusionMatrix(hmsMulticlassMLRes_BDN$resPredAll$pred, 
                       hmsMulticlassMLRes_BDN$resPredAll$obs)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction BRCA HNSC
# BRCA   15    6
# HNSC    4   70
# 
# Accuracy : 0.8947          
# 95% CI : (0.8149, 0.9484)
# No Information Rate : 0.8             
# P-Value [Acc > NIR] : 0.01053         
# 
# Kappa : 0.6835          
# 
# Mcnemar's Test P-Value : 0.75183         
#                                           
#             Sensitivity : 0.7895          
#             Specificity : 0.9211          
#          Pos Pred Value : 0.7143          
#          Neg Pred Value : 0.9459          
#              Prevalence : 0.2000          
#          Detection Rate : 0.1579          
#    Detection Prevalence : 0.2211          
#       Balanced Accuracy : 0.8553          
#                                           
#        'Positive' Class : BRCA


#----------------------------------------------------------#
# MDA-based PT 2-class machine learning
#----------------------------------------------------------#

# Call multiclass function
mdaMulticlassMLRes_PT <- mlMulticlass(metaData = metadataSamplesMergedRT_PT_MDA,
                                       countData = countMergedRT_PT_MDA,
                                       sampleType = "Primary Tumor",
                                       seqCenter = "MD Anderson - Institute for Applied Cancer Science",
                                       modelType = "xgbTree",
                                       numResampleIter = 1,
                                       numKFold = 10)

# Examine overall performance:
mdaMulticlassMLRes_PT$perfCombinedAll
# logLoss               AUC             prAUC          Accuracy 
# 0.3214088         0.9807988         0.9201309         0.9397590 
# Kappa                F1       Sensitivity       Specificity 
# 0.8380804         0.9600000         0.9677419         0.8571429 
# Pos_Pred_Value    Neg_Pred_Value         Precision            Recall 
# 0.9523810         0.9000000         0.9523810         0.9677419 
# Detection_Rate Balanced_Accuracy 
# 0.7228916         0.9124424 

# Examine per-CV performance:
mdaMulticlassMLRes_PT$repPerfCombinedAll

# Create confusion matrix
caret::confusionMatrix(mdaMulticlassMLRes_PT$resPredAll$pred, 
                       mdaMulticlassMLRes_PT$resPredAll$obs)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction BLCA HNSC
# BLCA   60    3
# HNSC    2   18
# 
# Accuracy : 0.9398         
# 95% CI : (0.865, 0.9802)
# No Information Rate : 0.747          
# P-Value [Acc > NIR] : 4.815e-06      
# 
# Kappa : 0.8381         
# 
# Mcnemar's Test P-Value : 1              
#                                          
#             Sensitivity : 0.9677         
#             Specificity : 0.8571         
#          Pos Pred Value : 0.9524         
#          Neg Pred Value : 0.9000         
#              Prevalence : 0.7470         
#          Detection Rate : 0.7229         
#    Detection Prevalence : 0.7590         
#       Balanced Accuracy : 0.9124         
#                                          
#        'Positive' Class : BLCA 

#----------------------------------------------------------#
# Broad-based PT 2-class machine learning
#----------------------------------------------------------#

# Call multiclass function
broadMulticlassMLRes_PT <- mlMulticlass(metaData = metadataSamplesMergedRT_PT_Broad,
                                      countData = countMergedRT_PT_Broad,
                                      sampleType = "Primary Tumor",
                                      seqCenter = "Broad Institute of MIT and Harvard",
                                      modelType = "xgbTree",
                                      numResampleIter = 1,
                                      numKFold = 10)

# Examine overall performance:
broadMulticlassMLRes_PT$perfCombinedAll
# logLoss               AUC             prAUC          Accuracy 
# 0.4070960         0.9294355         0.8871614         0.8510638 
# Kappa                F1       Sensitivity       Specificity 
# 0.6827387         0.8000000         0.8750000         0.8387097 
# Pos_Pred_Value    Neg_Pred_Value         Precision            Recall 
# 0.7368421         0.9285714         0.7368421         0.8750000 
# Detection_Rate Balanced_Accuracy 
# 0.2978723         0.8568548

# Examine per-CV performance:
broadMulticlassMLRes_PT$repPerfCombinedAll

# Create confusion matrix
caret::confusionMatrix(broadMulticlassMLRes_PT$resPredAll$pred, 
                       broadMulticlassMLRes_PT$resPredAll$obs)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction BLCA HNSC
# BLCA   14    5
# HNSC    2   26
# 
# Accuracy : 0.8511         
# 95% CI : (0.7169, 0.938)
# No Information Rate : 0.6596         
# P-Value [Acc > NIR] : 0.002843       
# 
# Kappa : 0.6827         
# 
# Mcnemar's Test P-Value : 0.449692       
#                                          
#             Sensitivity : 0.8750         
#             Specificity : 0.8387         
#          Pos Pred Value : 0.7368         
#          Neg Pred Value : 0.9286         
#              Prevalence : 0.3404         
#          Detection Rate : 0.2979         
#    Detection Prevalence : 0.4043         
#       Balanced Accuracy : 0.8569         
#                                          
#        'Positive' Class : BLCA

#----------------------------------------------------------#
# Merge ROC curves for 2 class-examples
#----------------------------------------------------------#

hmsMulticlassMLRes_PT$resPredAll
hmsMulticlassMLRes_BDN
mdaMulticlassMLRes_PT
broadMulticlassMLRes_PT

## Separate out HMS multiclass into 2 class to plot ROCs
# HMS BLCA vs all others
hmsMulticlassMLRes_PT_BLCA <- hmsMulticlassMLRes_PT$resPredAll %>%
  mutate(obs = ifelse(obs=="BLCA","BLCA","Other"))

hmsMulticlassMLRes_PT_BLCA_ROC <- roc(hmsMulticlassMLRes_PT_BLCA$obs ~ 
                                        hmsMulticlassMLRes_PT_BLCA$BLCA, 
                                      levels = c("Other","BLCA"),
                                  direction = "<")
# HMS BRCA vs all others
hmsMulticlassMLRes_PT_BRCA <- hmsMulticlassMLRes_PT$resPredAll %>%
  mutate(obs = ifelse(obs=="BRCA","BRCA","Other"))

hmsMulticlassMLRes_PT_BRCA_ROC <- roc(hmsMulticlassMLRes_PT_BRCA$obs ~ 
                                        hmsMulticlassMLRes_PT_BRCA$BRCA, 
                                      levels = c("Other","BRCA"),
                                      direction = "<")
# HMS HNSC vs all others
hmsMulticlassMLRes_PT_HNSC <- hmsMulticlassMLRes_PT$resPredAll %>%
  mutate(obs = ifelse(obs=="HNSC","HNSC","Other"))

hmsMulticlassMLRes_PT_HNSC_ROC <- roc(hmsMulticlassMLRes_PT_HNSC$obs ~ 
                                        hmsMulticlassMLRes_PT_HNSC$HNSC, 
                                      levels = c("Other","HNSC"),
                                      direction = "<")

## HMS BDN
hmsMulticlassMLRes_BDN_ROC <- roc(hmsMulticlassMLRes_BDN$resPredAll$obs ~ 
                                    hmsMulticlassMLRes_BDN$resPredAll$BRCA, 
                                      levels = c("HNSC","BRCA"),
                                      direction = "<")

## MDA PT
mdaMulticlassMLRes_PT_ROC <- roc(mdaMulticlassMLRes_PT$resPredAll$obs ~ 
                                    mdaMulticlassMLRes_PT$resPredAll$BLCA, 
                                  levels = c("HNSC","BLCA"),
                                  direction = "<")

## Broad PT
broadMulticlassMLRes_PT_ROC <- roc(broadMulticlassMLRes_PT$resPredAll$obs ~ 
                                   broadMulticlassMLRes_PT$resPredAll$BLCA, 
                                 levels = c("HNSC","BLCA"),
                                 direction = "<")

## PT combined ROC plot
rocList_PT <- list(HMS_BLCA = hmsMulticlassMLRes_PT_BLCA_ROC,
                   HMS_BRCA = hmsMulticlassMLRes_PT_BRCA_ROC,
                   HMS_HNSC = hmsMulticlassMLRes_PT_HNSC_ROC,
                   MDA = mdaMulticlassMLRes_PT_ROC,
                   Broad = broadMulticlassMLRes_PT_ROC)
legendTopYPos <- 0.25
ggroc(rocList_PT, size = 1, legacy.axes = TRUE) + 
  theme_pubr() + 
  labs(x = "False Positive Rate", y = "True Positive Rate") + 
  coord_fixed(ratio = 1) +
  scale_color_nejm() +
  theme(legend.position = "none", text=element_text(size=16)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", linetype="dashed") +
  # HMS BLCA
  annotate("segment", x = 0.2, xend = 0.3, y = legendTopYPos, yend = legendTopYPos, color = "#BC3C29FF") +
  annotate("text", x = 0.32, y = legendTopYPos, color = "#BC3C29FF", 
           label = paste0("HMS BLCA: ",100*round(hmsMulticlassMLRes_PT_BLCA_ROC$auc,4),"%"), hjust = 0, size=5) +
  # HMS BRCA
  annotate("segment", x = 0.2, xend = 0.3, y = legendTopYPos-0.07*1, yend = legendTopYPos-0.07*1, color = "#0072B5FF") +
  annotate("text", x = 0.32, y = legendTopYPos-0.07*1, color = "#0072B5FF", 
           label = paste0("HMS BRCA: ",100*round(hmsMulticlassMLRes_PT_BRCA_ROC$auc,4),"%"), hjust = 0, size=5)  +
  # HMS HNSC
  annotate("segment", x = 0.2, xend = 0.3, y = legendTopYPos-0.07*2, yend = legendTopYPos-0.07*2, color = "#E18727FF") +
  annotate("text", x = 0.32, y = legendTopYPos-0.07*2, color = "#E18727FF", 
           label = paste0("HMS HNSC: ",100*round(hmsMulticlassMLRes_PT_HNSC_ROC$auc,4),"%"), hjust = 0, size=5) +
  # MDA
  annotate("segment", x = 0.2, xend = 0.3, y = legendTopYPos-0.07*3, yend = legendTopYPos-0.07*3, color = "#20854EFF") +
  annotate("text", x = 0.32, y = legendTopYPos-0.07*3, color = "#20854EFF", 
           label = paste0("MDA: ",100*round(mdaMulticlassMLRes_PT_ROC$auc,4),"%"), hjust = 0, size=5) +
  # Broad
  annotate("segment", x = 0.2, xend = 0.3, y = legendTopYPos-0.07*4, yend = legendTopYPos-0.07*4, color = "#7876B1FF") +
  annotate("text", x = 0.32, y = legendTopYPos-0.07*4, color = "#7876B1FF", 
           label = paste0("Broad: ",100*round(broadMulticlassMLRes_PT_ROC$auc,4),"%"), hjust = 0, size=5)
ggsave(filename = "Figures/roc_PT_5Aug23.jpeg",
       dpi = "retina", units = "in", width = 4, height = 4)

## BDN ROC plot
rocList_BDN <- list(HMS = hmsMulticlassMLRes_BDN_ROC)
legendTopYPos <- 0.15
ggroc(rocList_BDN, size = 1, legacy.axes = TRUE) + 
  theme_pubr() + 
  labs(x = "False Positive Rate", y = "True Positive Rate") + 
  coord_fixed(ratio = 1) +
  scale_color_nejm() +
  theme(legend.position = "none", text=element_text(size=16)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", linetype="dashed") +
  # HMS
  annotate("segment", x = 0.2, xend = 0.3, y = legendTopYPos, yend = legendTopYPos, color = "#BC3C29FF") +
  annotate("text", x = 0.32, y = legendTopYPos, color = "#BC3C29FF", 
           label = paste0("HMS BDN: ",100*round(hmsMulticlassMLRes_BDN_ROC$auc,4),"%"), hjust = 0, size=5)
ggsave(filename = "Figures/roc_BDN_5Aug23.jpeg",
       dpi = "retina", units = "in", width = 4, height = 4)

#----------------------------------------------------------#
# Save data for Qiime PT data
#----------------------------------------------------------#

## Create mock taxa file for Qiime containing all genera
taxaFile <- data.frame(`Feature ID` = colnames(countMergedRT),
                       Taxon = colnames(countMergedRT),
                       check.names = FALSE)
write.table(taxaFile,
            file = "./Qiime_RT/Qiime_input_data/taxa.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

## Save metadata subsets by seqcenter
write.table(metadataSamplesMergedRT_PT_HMS,
            file = "./Qiime_RT/Qiime_input_data/metadataSamplesMergedRT_PT_HMS.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metadataSamplesMergedRT_PT_Broad,
            file = "./Qiime_RT/Qiime_input_data/metadataSamplesMergedRT_PT_Broad.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metadataSamplesMergedRT_PT_MDA,
            file = "./Qiime_RT/Qiime_input_data/metadataSamplesMergedRT_PT_MDA.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metadataSamplesMergedRT_BDN_HMS,
            file = "./Qiime_RT/Qiime_input_data/metadataSamplesMergedRT_BDN_HMS.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

## Save raw count data subsets by seqcenter
# HMS PT
countMergedRT_PT_HMS_BIOM <- make_biom(t(countMergedRT_PT_HMS))
write_biom(countMergedRT_PT_HMS_BIOM, 
           biom_file = "./Qiime_RT/Qiime_input_data/countMergedRT_PT_HMS.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMergedRT_PT_HMS))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     283     597   26994    8517  545212

# Broad PT
countMergedRT_PT_Broad_BIOM <- make_biom(t(countMergedRT_PT_Broad))
write_biom(countMergedRT_PT_Broad_BIOM, 
           biom_file = "./Qiime_RT/Qiime_input_data/countMergedRT_PT_Broad.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMergedRT_PT_Broad))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 22.0    90.0   172.0  1793.6   370.5 70423.0

# MDA PT
countMergedRT_PT_MDA_BIOM <- make_biom(t(countMergedRT_PT_MDA))
write_biom(countMergedRT_PT_MDA_BIOM, 
           biom_file = "./Qiime_RT/Qiime_input_data/countMergedRT_PT_MDA.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMergedRT_PT_MDA))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 114.0    244.5    336.0  16005.4    816.0 544568.0

# HMS BDN
countMergedRT_BDN_HMS_BIOM <- make_biom(t(countMergedRT_BDN_HMS))
write_biom(countMergedRT_BDN_HMS_BIOM, 
           biom_file = "./Qiime_RT/Qiime_input_data/countMergedRT_BDN_HMS.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMergedRT_BDN_HMS))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 55     204     328   13930     653  829677

##--> Qiime commands were run using ./Qiime_RT/qiime2_tcga_analyses_31July23.ipynb
## using the Qiime conda env version qiime2-2022.2
##--> Can view outputted .qzv files on https://view.qiime2.org/ 

#----------------------------------------------------------#
# ANCOM-BC
#----------------------------------------------------------#





