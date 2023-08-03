# tcga_gihawi_rebuttal_WIS_subset_31July23.R
# Author: Greg Poore
# Date: Aug 3, 2023
# Purpose: To explore cancer type differences in data released by 
# Gihawi et al. 2023 bioRxiv: https://www.biorxiv.org/content/10.1101/2023.07.28.550993v1
# only using bacterial and fungal genera found in tumors by the Weizmann Institute of Science (WIS)

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
# Load Weizmann (WIS) data summarized at genus level and subset Kraken features
# NOTE: WIS genera come from highly decontaminated 16S and ITS data, as published 
# in Nejman et al. 2020 Science and Narunsky-Haziza et al. 2022 Cell
#------------------------------------------------------#
## Load WIS genera from 'hit list'
load("Supporting_files/wis-bacteria-fungi-genera-species-bio-24July22.RData", verbose = TRUE)

wisGenera <- data.frame(tax_table(psWzBacteriaAndFungi_genus_Bio))
wisGeneraKnown <- wisGenera %>% 
  filter(!grepl("Unknown",genus)) %>%
  filter(!grepl("other",genus)) %>%
  droplevels()

wisGeneraKnownUnique <- unique(wisGeneraKnown$genus)
length(unique(wisGeneraKnownUnique)) # 437

## Subset tables
sum(colnames(blcaKraken) %in% wisGeneraKnownUnique) # 151
sum(colnames(hnscKraken) %in% wisGeneraKnownUnique) # 161
sum(colnames(brcaKraken) %in% wisGeneraKnownUnique) # 159

blcaKrakenWIS <- blcaKraken[,colnames(blcaKraken) %in% wisGeneraKnownUnique]
hnscKrakenWIS <- hnscKraken[,colnames(hnscKraken) %in% wisGeneraKnownUnique]
brcaKrakenWIS <- brcaKraken[,colnames(brcaKraken) %in% wisGeneraKnownUnique]

#------------------------------------------------------#
# Merge tables
# - Selected conservative method to retain only overlapping features
# (although that reduce some effect size since the files are
# separated by cancer type)
# - Note that WIS did not contain "Homo" genus calls, so they are inherently excluded
#------------------------------------------------------#
sharedFeat <- Reduce(intersect, list(colnames(blcaKrakenWIS),
                                     colnames(hnscKrakenWIS),
                                     colnames(brcaKrakenWIS)))

countMergedWIS <- smartbind(blcaKrakenWIS[,sharedFeat],
                            hnscKrakenWIS[,sharedFeat],
                            brcaKrakenWIS[,sharedFeat])
# Missing values after the merge should be converted to 0s
countMergedWIS[is.na(countMergedWIS)] <- 0 
rownames(countMergedWIS) <- c(rownames(blcaKrakenWIS),
                           rownames(hnscKrakenWIS),
                           rownames(brcaKrakenWIS))
dim(countMergedWIS) # 728 149

# Subset metadata to samples reflected by Gihawi et al. 2023
metadataSamplesMergedWIS <- metadataSamplesAll[rownames(countMergedWIS),]

# Subset to primary tumor ("PT") and blod derived normal ("BDN") samples
metadataSamplesMergedWIS_PT <- metadataSamplesMergedWIS %>% 
  filter(sample_type == "Primary Tumor") %>% droplevels()
metadataSamplesMergedWIS_BDN <- metadataSamplesMergedWIS %>% 
  filter(sample_type == "Blood Derived Normal") %>% droplevels()

#------------------------------------------------------#
# Examine PT sequencing center subsets
#
# NOTE: The goal is to have individual seq-center subsets
# to avoid needing to batch correction. This allows us to use
# the raw data for comparisons.
#------------------------------------------------------#

# Count PT subsets
metadataSamplesMergedWIS_PT %>% count(data_submitting_center_label, investigation)
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

metadataSamplesMergedWIS_PT_HMS <- metadataSamplesMergedWIS_PT %>%
  filter(data_submitting_center_label == "Harvard Medical School") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  droplevels()
countMergedWIS_PT_HMS <- countMergedWIS[metadataSamplesMergedWIS_PT_HMS$sampleid,]

metadataSamplesMergedWIS_PT_Broad <- metadataSamplesMergedWIS_PT %>%
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  droplevels()
countMergedWIS_PT_Broad <- countMergedWIS[metadataSamplesMergedWIS_PT_Broad$sampleid,]

metadataSamplesMergedWIS_PT_MDA <- metadataSamplesMergedWIS_PT %>%
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  filter(investigation != "TCGA-BRCA") %>% # Drop the 2 samples
  droplevels()
countMergedWIS_PT_MDA <- countMergedWIS[metadataSamplesMergedWIS_PT_MDA$sampleid,]

#------------------------------------------------------#
# Examine BDN sequencing center subsets
#
# NOTE: The goal is to have individual seq-center subsets
# to avoid needing batch correction. This allows us to use
# the raw data for comparisons.
#------------------------------------------------------#
metadataSamplesMergedWIS_BDN %>% count(data_submitting_center_label, investigation)
# data_submitting_center_label investigation  n
# 1                         Baylor College of Medicine     TCGA-HNSC 28
# 2                 Broad Institute of MIT and Harvard     TCGA-HNSC 14
# 3                             Harvard Medical School     TCGA-BRCA 19
# 4                             Harvard Medical School     TCGA-HNSC 76
# 5 MD Anderson - Institute for Applied Cancer Science     TCGA-HNSC 22
# 6           Washington University School of Medicine     TCGA-BRCA 87

##--> Harvard Med School (HMS) has 2 cancer types
##--> Everywhere else cannot be used without batch correction, so they will be excluded

metadataSamplesMergedWIS_BDN_HMS <- metadataSamplesMergedWIS_BDN %>%
  filter(data_submitting_center_label == "Harvard Medical School") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  droplevels()
countMergedWIS_BDN_HMS <- countMergedWIS[metadataSamplesMergedWIS_BDN_HMS$sampleid,]

#------------------------------------------------------#
# For completeness - sequencing platform
#
# - NOTE: 1 sample was processed on a MiSeq, but it was at WashU,
# which is not included in the PT or BDN subsets (so ignore)
#------------------------------------------------------#
metadataSamplesMergedWIS %>% count(platform, data_submitting_center_label)
# platform                       data_submitting_center_label   n
# 1 Illumina HiSeq                         Baylor College of Medicine  65
# 2 Illumina HiSeq                 Broad Institute of MIT and Harvard  67
# 3 Illumina HiSeq                             Harvard Medical School 284
# 4 Illumina HiSeq MD Anderson - Institute for Applied Cancer Science 113
# 5 Illumina HiSeq      MD Anderson - RPPA Core Facility (Proteomics)   1
# 6 Illumina HiSeq           Washington University School of Medicine 197
# 7 Illumina MiSeq           Washington University School of Medicine   1

#----------------------------------------------------------#
# PT multi-class machine learning
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
mlMulticlass <- function(metaData = metadataSamplesMergedWIS_PT_HMS,
                         countData = countMergedWIS_PT_HMS,
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
  
  ## Save predictions and perf
  # Preds
  baseFilename <- paste0("multiclassCV_WIS_",st,"_",sc,"_k",numKFold,"_modelType_",modelType)
  write.csv(resPredAll, file = paste0("ML_results/pred",baseFilename,".csv"))
  save(resPredAll, file = paste0("ML_results/pred",baseFilename,".RData"))
  # Overall perf
  write.csv(perfCombinedAll, file = paste0("ML_results/perfCombinedAll",baseFilename,".csv"))
  save(perfCombinedAll, file = paste0("ML_results/perfCombinedAll",baseFilename,".RData"))
  # Rep perf
  write.csv(resPredAll, file = paste0("ML_results/repPerfCombinedAll",baseFilename,".csv"))
  save(resPredAll, file = paste0("ML_results/repPerfCombinedAll",baseFilename,".RData"))
  
  rm(mlModel)
  
  res <- list(resPredAll=resPredAll,
              perfCombinedAll=perfCombinedAll,
              repPerfCombinedAll=repPerfCombinedAll)
  return(res)
}


# Call multiclass function
hmsMulticlassMLRes_PT <- mlMulticlass()

# Examine overall performance:
hmsMulticlassMLRes_PT$perfCombinedAll
# logLoss                    AUC                  prAUC 
# 0.4944489              0.9919219              0.9560662 
# Accuracy                  Kappa                Mean_F1 
# 0.9290323              0.8771437              0.9217252 
# Mean_Sensitivity       Mean_Specificity    Mean_Pos_Pred_Value 
# 0.9107671              0.9616300              0.9400454 
# Mean_Neg_Pred_Value         Mean_Precision            Mean_Recall 
# 0.9589059              0.9400454              0.9107671 
# Mean_Detection_Rate Mean_Balanced_Accuracy 
# 0.3096774              0.9361986 

# Examine per-CV performance:
hmsMulticlassMLRes_PT$repPerfCombinedAll

# Create confusion matrix
caret::confusionMatrix(hmsMulticlassMLRes_PT$resPredAll$pred, 
                       hmsMulticlassMLRes_PT$resPredAll$obs)

# Confusion Matrix and Statistics
# 
# Reference
# Prediction BLCA BRCA HNSC
# BLCA   48    3    5
# BRCA    0   15    0
# HNSC    3    1   80
# 
# Overall Statistics
# 
# Accuracy : 0.9226          
# 95% CI : (0.8687, 0.9594)
# No Information Rate : 0.5484          
# P-Value [Acc > NIR] : <2e-16          
# 
# Kappa : 0.8647          
# 
# Mcnemar's Test P-Value : 0.2123          
# 
# Statistics by Class:
# 
#                      Class: BLCA Class: BRCA Class: HNSC
# Sensitivity               0.9412     0.78947      0.9412
# Specificity               0.9231     1.00000      0.9429
# Pos Pred Value            0.8571     1.00000      0.9524
# Neg Pred Value            0.9697     0.97143      0.9296
# Prevalence                0.3290     0.12258      0.5484
# Detection Rate            0.3097     0.09677      0.5161
# Detection Prevalence      0.3613     0.09677      0.5419
# Balanced Accuracy         0.9321     0.89474      0.9420

#----------------------------------------------------------#
# BDN 2-class machine learning
#
# - Since HMS is the only seq-center subset that had >1 cancer types,
# it will be the only multi-class ML instance.
# - Although not shown here, this could also be applied to the 
# PT 2-class subsets above
#----------------------------------------------------------#

# Call multiclass function
hmsMulticlassMLRes_BDN <- mlMulticlass(metaData = metadataSamplesMergedWIS_BDN_HMS,
                                   countData = countMergedWIS_BDN_HMS,
                                   sampleType = "Blood Derived Normal",
                                   seqCenter = "Harvard Medical School",
                                   modelType = "xgbTree",
                                   numResampleIter = 1,
                                   numKFold = 10)

# Examine overall performance:
hmsMulticlassMLRes_BDN$perfCombinedAll
# logLoss               AUC             prAUC          Accuracy 
# 0.3366396         0.9799169         0.9111785         0.9473684 
# Kappa                F1       Sensitivity       Specificity 
# 0.8251748         0.8571429         0.7894737         0.9868421 
# Pos_Pred_Value    Neg_Pred_Value         Precision            Recall 
# 0.9375000         0.9493671         0.9375000         0.7894737 
# Detection_Rate Balanced_Accuracy 
# 0.1578947         0.8881579 

# Examine per-CV performance:
hmsMulticlassMLRes_BDN$repPerfCombinedAll

# Create confusion matrix
caret::confusionMatrix(hmsMulticlassMLRes_BDN$resPredAll$pred, 
                       hmsMulticlassMLRes_BDN$resPredAll$obs)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction BRCA HNSC
# BRCA   15    1
# HNSC    4   75
# 
# Accuracy : 0.9474          
# 95% CI : (0.8814, 0.9827)
# No Information Rate : 0.8             
# P-Value [Acc > NIR] : 4.444e-05       
# 
# Kappa : 0.8252          
# 
# Mcnemar's Test P-Value : 0.3711          
#                                           
#             Sensitivity : 0.7895          
#             Specificity : 0.9868          
#          Pos Pred Value : 0.9375          
#          Neg Pred Value : 0.9494          
#              Prevalence : 0.2000          
#          Detection Rate : 0.1579          
#    Detection Prevalence : 0.1684          
#       Balanced Accuracy : 0.8882          
#                                           
#        'Positive' Class : BRCA



