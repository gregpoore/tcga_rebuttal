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
# - Since the tables are already intersected with the WIS features, they 
# are simply joined here. Any missing values after joining the tables are filled 
# in with 0s.
# - Note that WIS did not contain "Homo" genus calls, so they are inherently excluded
#------------------------------------------------------#
countMergedWIS <- smartbind(blcaKrakenWIS,
                            hnscKrakenWIS,
                            brcaKrakenWIS)
# Missing values after the merge should be converted to 0s
countMergedWIS[is.na(countMergedWIS)] <- 0 
rownames(countMergedWIS) <- c(rownames(blcaKrakenWIS),
                           rownames(hnscKrakenWIS),
                           rownames(brcaKrakenWIS))
dim(countMergedWIS) # 728 149 | 728 161

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
# Save data for Qiime PT data
#----------------------------------------------------------#

## Create mock taxa file for Qiime containing all genera
taxaFile <- data.frame(`Feature ID` = colnames(countMergedWIS),
                       Taxon = colnames(countMergedWIS),
                       check.names = FALSE)
write.table(taxaFile,
            file = "./Qiime_WIS/Qiime_input_data/taxa.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

## Save metadata subsets by seqcenter
write.table(metadataSamplesMergedWIS_PT_HMS,
            file = "./Qiime_WIS/Qiime_input_data/metadataSamplesMergedWIS_PT_HMS.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metadataSamplesMergedWIS_PT_Broad,
            file = "./Qiime_WIS/Qiime_input_data/metadataSamplesMergedWIS_PT_Broad.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metadataSamplesMergedWIS_PT_MDA,
            file = "./Qiime_WIS/Qiime_input_data/metadataSamplesMergedWIS_PT_MDA.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metadataSamplesMergedWIS_BDN_HMS,
            file = "./Qiime_WIS/Qiime_input_data/metadataSamplesMergedWIS_BDN_HMS.txt",
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

## Save raw count data subsets by seqcenter
# HMS PT
countMergedWIS_PT_HMS_BIOM <- make_biom(t(countMergedWIS_PT_HMS))
write_biom(countMergedWIS_PT_HMS_BIOM, 
           biom_file = "./Qiime_WIS/Qiime_input_data/countMergedWIS_PT_HMS.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMergedWIS_PT_HMS))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0     359.5     802.0   50367.0   17343.5 1858680.0

# Broad PT
countMergedWIS_PT_Broad_BIOM <- make_biom(t(countMergedWIS_PT_Broad))
write_biom(countMergedWIS_PT_Broad_BIOM, 
           biom_file = "./Qiime_WIS/Qiime_input_data/countMergedWIS_PT_Broad.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMergedWIS_PT_Broad))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.0   191.5   417.0  2942.3  1491.5 78342.0

# MDA PT
countMergedWIS_PT_MDA_BIOM <- make_biom(t(countMergedWIS_PT_MDA))
write_biom(countMergedWIS_PT_MDA_BIOM, 
           biom_file = "./Qiime_WIS/Qiime_input_data/countMergedWIS_PT_MDA.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMergedWIS_PT_MDA))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 136.0    270.5    374.0  26416.7   1561.0 740872.0 

# HMS BDN
countMergedWIS_BDN_HMS_BIOM <- make_biom(t(countMergedWIS_BDN_HMS))
write_biom(countMergedWIS_BDN_HMS_BIOM, 
           biom_file = "./Qiime_WIS/Qiime_input_data/countMergedWIS_BDN_HMS.biom")
# Check distribution of sample read counts (use 1st quartile for alpha div)
summary(rowSums(countMergedWIS_BDN_HMS))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 118     284     456   18790    1027 1083150

##--> Qiime commands were run using ./Qiime_WIS/qiime2_tcga_analyses_31July23.ipynb
## using the Qiime conda env version qiime2-2022.2
##--> Can view outputted .qzv files on https://view.qiime2.org/ 

#----------------------------------------------------------#
# Plot alpha diversities for the HMS PT and BDN subsets
#----------------------------------------------------------#
# Instructions: 
# - Upload the following 2 files to https://view.qiime2.org/
# (1) ./Qiime_WIS/core_metrics_pt_hms/shannon_vector_significance.qzv
# (2) ./Qiime_WIS/core_metrics_bdn_hms/shannon_vector_significance.qzv
# - Then click "Download raw data as TSV" (just below and left of the plot)
# - Save in the respective folder (eg ./Qiime_WIS/core_metrics_pt_hms/) as a CSV file


##------------PT------------##
# HMS
alphaDivShannon_HMS_PT <- read.csv("Qiime_WIS/core_metrics_pt_hms/shannon_vector_significance_raw_values.csv",
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
ggsave(filename = "Figures/alphaDivShannonPlotWIS_HMS_PT.jpeg",
       dpi = "retina", units = "in", height = 4, width = 2.5)

##------------BDN------------##
# HMS
alphaDivShannon_HMS_BDN <- read.csv("Qiime_WIS/core_metrics_bdn_hms/shannon_vector_significance_raw_values.csv",
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
ggsave(filename = "Figures/alphaDivShannonPlotWIS_HMS_BDN.jpeg",
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
# 0.4998314              0.9907398              0.9499644 
# Accuracy                  Kappa                Mean_F1 
# 0.9225806              0.8646682              0.9087649 
# Mean_Sensitivity       Mean_Specificity    Mean_Pos_Pred_Value 
# 0.8906089              0.9553114              0.9365079 
# Mean_Neg_Pred_Value         Mean_Precision            Mean_Recall 
# 0.9569010              0.9365079              0.8906089 
# Mean_Detection_Rate Mean_Balanced_Accuracy 
# 0.3075269              0.9229601 

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
# 0.3289573         0.9750693         0.9329673         0.9473684 
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



