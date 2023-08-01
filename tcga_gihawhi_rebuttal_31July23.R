# tcga_gihawi_rebuttal_31July23.R
# Author: Greg Poore
# Date: July 31, 2023
# Purpose: To explore cancer type differences in data released by 
# Gihawi et al. 2023 bioRxiv: https://www.biorxiv.org/content/10.1101/2023.07.28.550993v1

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
# Merge tables
# - Selected conservative method to retain only overlapping features
# (although that reduce some effect size since the files are
# separated by cancer type)
# - Also removed "Homo" genus calls to ensure cancer type differences
# were coming from non-human organisms in the database
#------------------------------------------------------#
sharedFeat <- Reduce(intersect, list(colnames(blcaKraken),
                                     colnames(hnscKraken),
                                     colnames(brcaKraken)))
sharedFeatNoHuman <- sharedFeat[-which(sharedFeat=="Homo")] # Remove human calls

# Merge matrices using smartbind
countMerged <- smartbind(blcaKraken[,sharedFeatNoHuman],
                         hnscKraken[,sharedFeatNoHuman],
                         brcaKraken[,sharedFeatNoHuman])
# Missing values after the merge should be converted to 0s
countMerged[is.na(countMerged)] <- 0 
rownames(countMerged) <- c(rownames(blcaKraken),
                           rownames(hnscKraken),
                           rownames(brcaKraken))
dim(countMerged) # 728 847

# Subset metadata to samples reflected by Gihawi et al. 2023
metadataSamplesMerged <- metadataSamplesAll[rownames(countMerged),]

# Subset to primary tumor ("PT") and blod derived normal ("BDN") samples
metadataSamplesMerged_PT <- metadataSamplesMerged %>% 
  filter(sample_type == "Primary Tumor") %>% droplevels()
metadataSamplesMerged_BDN <- metadataSamplesMerged %>% 
  filter(sample_type == "Blood Derived Normal") %>% droplevels()

#------------------------------------------------------#
# Examine PT sequencing center subsets
#
# NOTE: The goal is to have individual seq-center subsets
# to avoid needing to batch correction. This allows us to use
# the raw data for comparisons.
#------------------------------------------------------#

# Count PT subsets
metadataSamplesMerged_PT %>% count(data_submitting_center_label, investigation)
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

metadataSamplesMerged_PT_HMS <- metadataSamplesMerged_PT %>%
  filter(data_submitting_center_label == "Harvard Medical School") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  droplevels()
countMerged_PT_HMS <- countMerged[metadataSamplesMerged_PT_HMS$sampleid,]

metadataSamplesMerged_PT_Broad <- metadataSamplesMerged_PT %>%
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  droplevels()
countMerged_PT_Broad <- countMerged[metadataSamplesMerged_PT_Broad$sampleid,]

metadataSamplesMerged_PT_MDA <- metadataSamplesMerged_PT %>%
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  filter(investigation != "TCGA-BRCA") %>% # Drop the 2 samples
  droplevels()
countMerged_PT_MDA <- countMerged[metadataSamplesMerged_PT_MDA$sampleid,]

#------------------------------------------------------#
# Examine BDN sequencing center subsets
#
# NOTE: The goal is to have individual seq-center subsets
# to avoid needing batch correction. This allows us to use
# the raw data for comparisons.
#------------------------------------------------------#
metadataSamplesMerged_BDN %>% count(data_submitting_center_label, investigation)
# data_submitting_center_label investigation  n
# 1                         Baylor College of Medicine     TCGA-HNSC 28
# 2                 Broad Institute of MIT and Harvard     TCGA-HNSC 14
# 3                             Harvard Medical School     TCGA-BRCA 19
# 4                             Harvard Medical School     TCGA-HNSC 76
# 5 MD Anderson - Institute for Applied Cancer Science     TCGA-HNSC 22
# 6           Washington University School of Medicine     TCGA-BRCA 87

##--> Harvard Med School (HMS) has 2 cancer types
##--> Everywhere else cannot be used without batch correction, so they will be excluded

metadataSamplesMerged_BDN_HMS <- metadataSamplesMerged_BDN %>%
  filter(data_submitting_center_label == "Harvard Medical School") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>%
  droplevels()
countMerged_BDN_HMS <- countMerged[metadataSamplesMerged_BDN_HMS$sampleid,]

#------------------------------------------------------#
# For completeness - sequencing platform
#
# - NOTE: 1 sample was processed on a MiSeq, but it was at WashU,
# which is not included in the PT or BDN subsets (so ignore)
#------------------------------------------------------#
metadataSamplesMerged %>% count(platform, data_submitting_center_label)
# platform                       data_submitting_center_label   n
# 1 Illumina HiSeq                         Baylor College of Medicine  65
# 2 Illumina HiSeq                 Broad Institute of MIT and Harvard  67
# 3 Illumina HiSeq                             Harvard Medical School 284
# 4 Illumina HiSeq MD Anderson - Institute for Applied Cancer Science 113
# 5 Illumina HiSeq      MD Anderson - RPPA Core Facility (Proteomics)   1
# 6 Illumina HiSeq           Washington University School of Medicine 197
# 7 Illumina MiSeq           Washington University School of Medicine   1

#----------------------------------------------------------#
# Save data for Qiime PT data
#----------------------------------------------------------#

## Create mock taxa file for Qiime containing all genera
taxaFile <- data.frame(`Feature ID` = colnames(countMerged),
                       Taxon = colnames(countMerged),
                       check.names = FALSE)
write.table(taxaFile,
            file = "./Qiime/Qiime_input_data/taxa.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

## Save metadata subsets by seqcenter
write.table(metadataSamplesMerged_PT_HMS,
            file = "./Qiime/Qiime_input_data/metadataSamplesMerged_PT_HMS.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metadataSamplesMerged_PT_Broad,
            file = "./Qiime/Qiime_input_data/metadataSamplesMerged_PT_Broad.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metadataSamplesMerged_PT_MDA,
            file = "./Qiime/Qiime_input_data/metadataSamplesMerged_PT_MDA.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metadataSamplesMerged_BDN_HMS,
            file = "./Qiime/Qiime_input_data/metadataSamplesMerged_BDN_HMS.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

## Save raw count data subsets by seqcenter
# HMS PT
countMerged_PT_HMS_BIOM <- make_biom(t(countMerged_PT_HMS))
write_biom(countMerged_PT_HMS_BIOM, 
           biom_file = "./Qiime/Qiime_input_data/countMerged_PT_HMS.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMerged_PT_HMS))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 1.0     618.5    1799.0   51852.0   17551.5 1863068.0 

# Broad PT
countMerged_PT_Broad_BIOM <- make_biom(t(countMerged_PT_Broad))
write_biom(countMerged_PT_Broad_BIOM, 
           biom_file = "./Qiime/Qiime_input_data/countMerged_PT_Broad.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMerged_PT_Broad))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 77     396     733    8981    4636   81314 

# MDA PT
countMerged_PT_MDA_BIOM <- make_biom(t(countMerged_PT_MDA))
write_biom(countMerged_PT_MDA_BIOM, 
           biom_file = "./Qiime/Qiime_input_data/countMerged_PT_MDA.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMerged_PT_MDA))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 157.0    337.5    492.0  27181.9   1731.5 758959.0

# HMS BDN
countMerged_BDN_HMS_BIOM <- make_biom(t(countMerged_BDN_HMS))
write_biom(countMerged_BDN_HMS_BIOM, 
           biom_file = "./Qiime/Qiime_input_data/countMerged_BDN_HMS.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMerged_BDN_HMS))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 164.0     600.5    1607.0   22539.8    4843.0 1154723.0

##--> Qiime commands were run using ./Qiime/qiime2_tcga_analyses_31July23.ipynb
## using the Qiime conda env version qiime2-2022.2
##--> Can view outputted .qzv files on https://view.qiime2.org/ 

#----------------------------------------------------------#
# Plot alpha diversities for the HMS PT and BDN subsets
#----------------------------------------------------------#
# Instructions: 
# - Upload the following 2 files to https://view.qiime2.org/
# (1) ./Qiime/core_metrics_pt_hms/shannon_vector_significance.qzv
# (2) ./Qiime/core_metrics_bdn_hms/shannon_vector_significance.qzv
# - Then click "Download raw data as TSV" (just below and left of the plot)
# - Save in the respective folder (eg ./Qiime/core_metrics_pt_hms/) as a CSV file


##------------PT------------##
# HMS
alphaDivShannon_HMS_PT <- read.csv("Qiime/core_metrics_pt_hms/shannon_vector_significance_raw_values.csv",
                                   stringsAsFactors = FALSE, comment.char = "#")
alphaDivShannon_HMS_PT$investigation <- gsub("^TCGA\\-","",alphaDivShannon_HMS_PT$investigation)

# Plot
alphaDivShannon_HMS_PT %>%
  ggplot(aes(reorder(investigation, shannon_entropy, FUN=median), shannon_entropy, fill=investigation)) +
  geom_boxplot() + xlab("Primary tumors (HMS)") + ylab("Shannon Entropy") + 
  geom_quasirandom(method = "smiley", size=0.8) +
  theme_pubr(legend = "none") +
  rotate_x_text(0) +
  scale_fill_nejm() +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 1) + 
  stat_n_text()
ggsave(filename = "Figures/alphaDivShannonPlot_HMS_PT.jpeg",
       dpi = "retina", units = "in", height = 4, width = 2.5)

##------------BDN------------##
# HMS
alphaDivShannon_HMS_BDN <- read.csv("Qiime/core_metrics_bdn_hms/shannon_vector_significance_raw_values.csv",
                                   stringsAsFactors = FALSE, comment.char = "#")
alphaDivShannon_HMS_BDN$investigation <- gsub("^TCGA\\-","",alphaDivShannon_HMS_BDN$investigation)

# Plot
alphaDivShannon_HMS_BDN %>%
  ggplot(aes(reorder(investigation, shannon_entropy, FUN=median), shannon_entropy, fill=investigation)) +
  geom_boxplot() + xlab("Blood samples (HMS)") + ylab("Shannon Entropy") + 
  geom_quasirandom(method = "smiley", size=0.8) +
  theme_pubr(legend = "none") +
  rotate_x_text(0) +
  scale_fill_manual(values = c("#0072B5FF","#E18727FF")) +
  stat_compare_means(label.x.npc = 0, label.y.npc = 1) + 
  stat_n_text()
ggsave(filename = "Figures/alphaDivShannonPlot_HMS_BDN.jpeg",
       dpi = "retina", units = "in", height = 4, width = 2.5)

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
mlMulticlass <- function(metaData = metadataSamplesMerged_PT_HMS,
                         countData = countMerged_PT_HMS,
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
  baseFilename <- paste0("multiclassCV_",st,"_",sc,"_k",numKFold,"_modelType_",modelType)
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
# 0.4464495              0.9970843              0.9547654 
# Accuracy                  Kappa                Mean_F1 
# 0.9677419              0.9438569              0.9737518 
# Mean_Sensitivity       Mean_Specificity    Mean_Pos_Pred_Value 
# 0.9725490              0.9793040              0.9750388 
# Mean_Neg_Pred_Value         Mean_Precision            Mean_Recall 
# 0.9808144              0.9750388              0.9725490 
# Mean_Detection_Rate Mean_Balanced_Accuracy 
# 0.3225806              0.9759265 

# Examine per-CV performance:
hmsMulticlassMLRes_PT$repPerfCombinedAll

# Create confusion matrix
caret::confusionMatrix(hmsMulticlassMLRes_PT$resPredAll$pred, 
                       hmsMulticlassMLRes_PT$resPredAll$obs)

# Confusion Matrix and Statistics
# 
# Reference
# Prediction BLCA BRCA HNSC
# BLCA   48    0    2
# BRCA    0   19    0
# HNSC    3    0   83
# 
# Overall Statistics
# 
# Accuracy : 0.9677          
# 95% CI : (0.9263, 0.9894)
# No Information Rate : 0.5484          
# P-Value [Acc > NIR] : < 2.2e-16       
# 
# Kappa : 0.9439          
# 
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: BLCA Class: BRCA Class: HNSC
# Sensitivity               0.9412      1.0000      0.9765
# Specificity               0.9808      1.0000      0.9571
# Pos Pred Value            0.9600      1.0000      0.9651
# Neg Pred Value            0.9714      1.0000      0.9710
# Prevalence                0.3290      0.1226      0.5484
# Detection Rate            0.3097      0.1226      0.5355
# Detection Prevalence      0.3226      0.1226      0.5548
# Balanced Accuracy         0.9610      1.0000      0.9668

#----------------------------------------------------------#
# BDN 2-class machine learning
#
# - Since HMS is the only seq-center subset that had >1 cancer types,
# it will be the only multi-class ML instance.
# - Although not shown here, this could also be applied to the 
# PT 2-class subsets above
#----------------------------------------------------------#

# Call multiclass function
hmsMulticlassMLRes_BDN <- mlMulticlass(metaData = metadataSamplesMerged_BDN_HMS,
                                   countData = countMerged_BDN_HMS,
                                   sampleType = "Blood Derived Normal",
                                   seqCenter = "Harvard Medical School",
                                   modelType = "xgbTree",
                                   numResampleIter = 1,
                                   numKFold = 10)

# Examine overall performance:
hmsMulticlassMLRes_BDN$perfCombinedAll
# logLoss               AUC             prAUC          Accuracy 
# 0.2762710         0.9916898         0.9192374         0.9473684 
# Kappa                F1       Sensitivity       Specificity 
# 0.8387097         0.8717949         0.8947368         0.9605263 
# Pos_Pred_Value    Neg_Pred_Value         Precision            Recall 
# 0.8500000         0.9733333         0.8500000         0.8947368 
# Detection_Rate Balanced_Accuracy 
# 0.1789474         0.9276316 

# Examine per-CV performance:
hmsMulticlassMLRes_BDN$repPerfCombinedAll

# Create confusion matrix
caret::confusionMatrix(hmsMulticlassMLRes_BDN$resPredAll$pred, 
                       hmsMulticlassMLRes_BDN$resPredAll$obs)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction BRCA HNSC
# BRCA   17    3
# HNSC    2   73
# 
# Accuracy : 0.9474          
# 95% CI : (0.8814, 0.9827)
# No Information Rate : 0.8             
# P-Value [Acc > NIR] : 4.444e-05       
# 
# Kappa : 0.8387          
# 
# Mcnemar's Test P-Value : 1               
#                                           
#             Sensitivity : 0.8947          
#             Specificity : 0.9605          
#          Pos Pred Value : 0.8500          
#          Neg Pred Value : 0.9733          
#              Prevalence : 0.2000          
#          Detection Rate : 0.1789          
#    Detection Prevalence : 0.2105          
#       Balanced Accuracy : 0.9276          
#                                           
#        'Positive' Class : BRCA   



