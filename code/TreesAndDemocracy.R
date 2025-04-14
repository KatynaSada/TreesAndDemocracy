mac <- "/Users/katyna/Library/CloudStorage/OneDrive----/"
angel <- "C:/Users/arubio/OneDrive - Tecnun/"
computer <- mac
setwd(computer) # choose directory

# Loading Libraries and Packages --------------------------------------------
#The first part is to load all the libraries and external functions that are required for the analysis. # nolint # nolint

library(readxl)
library(RColorBrewer)
library(matrixStats)
library(partykit)
library(glmnet)

library(impute)
library(knitr)
library(reticulate)
library(readr)
library(limma)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(ggsci)
library(matrixStats)
# library(maftools)
# library(qvalue)
library(tibble)
library(Rtsne)
library(writexl)
library(matrixStats)
library(pheatmap)
library(DT)
library("dplyr")
library(tictoc)
library(ODT)

# add option to load to macos fonts
if (Sys.info()[1] == "Roboto") {
  loadfonts(device = "macOS")
}

source("TreesAndDemocracy/XAIfunctions.R")

#windowsFonts("Roboto" = windowsFont("Roboto"))

# Folder to save results
folder_dir <- "TreesAndDemocracy/output"

# Loading Data -------------------------------------------- 
# Load clinical data 
clinical <- read_excel("GithubProjects/interpretable-ai-review/data/input/BeatAML2/beataml_wv1to4_clinical.xlsx",
                       sheet = "summary"
)
clinical$dbgap_subject_id <- as.character(clinical$dbgap_subject_id)
# keep only first specimen collected of each subject
clinical <- clinical[!clinical$timeOfSampleCollectionRelativeToInclusion>0,]
# keep subjects that have dna and rna
clinical <- clinical[complete.cases(clinical[, c("dbgap_dnaseq_sample", "dbgap_rnaseq_sample")]),]
# make sure no samples are duplicated: keep only one rna and dna sample per subject (some have more)
sum(duplicated(clinical$dbgap_subject_id)) # no samples are duplicated
# clinical <- clinical[!duplicated(clinical$dbgap_subject_id),]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$dbgap_subject_id

#Then, we will load all the data for the analysis.
#Data selected was the BeatAML2 cohort for which we had mutational data, drug sensitivity data and gene expression data.

# Load gene expression data
expression <- read_excel("GithubProjects/interpretable-ai-review/data/input/BeatAML2/beataml_waves1to4_norm_exp_dbgap.xlsx")
expression <- as.data.frame(expression)
# set symbols as row names 
rownames(expression) <- expression$display_label 
# set dbgap_subject_id as colnames
rna_ids <- clinical$dbgap_rnaseq_sample[clinical$dbgap_rnaseq_sample  %in% colnames(expression)]
expression <- expression[,rna_ids]
subject_ids <-na.omit(clinical$dbgap_subject_id[match(colnames(expression),clinical$dbgap_rnaseq_sample)])
colnames(expression) <- subject_ids

# Load mutations data
gene_variants <- read_excel("GithubProjects/interpretable-ai-review/data/input/BeatAML2/beataml_wes_wv1to4_mutations_dbgap.xlsx")
colnames(gene_variants)[1]="dbgap_dnaseq_sample"
# add dbgap_subject_id and keep only samples in clinical (have also dna)
gene_variants <-merge(gene_variants,clinical[,c("dbgap_subject_id","dbgap_dnaseq_sample")], by="dbgap_dnaseq_sample")

# Load drug sensitivity data
drug_response <- read_excel("GithubProjects/interpretable-ai-review/data/input/BeatAML2/beataml_probit_curve_fits_v4_dbgap.xlsx")
drug_response$dbgap_subject_id <- as.character(drug_response$dbgap_subject_id)

# Pre-processing drug and mutational information  --------------------------------------------

# Now it is time to pre-process the data in order to obtain two different matrices, a matrix containing all patients and their mutations available and a second matrix containing all patients and their sensitivity to the different drugs.
# For doing so, we will impute the missing values in the drug matrix.

# DRUG RESPONSE MATRIX 
drug_names <- unique(drug_response$inhibitor)
subjects_ids <- clinical$dbgap_subject_id

# Build Drug Matrix 
drug_matrix <- matrix( # empty matrix
  data = "NA", nrow = length(drug_names),
  ncol = length(subjects_ids)
) 

colnames(drug_matrix) <- subjects_ids # columns are subjects
rownames(drug_matrix) <- drug_names # rows are drugs

for (subject in subjects_ids) { # add the ic50
  inds_subject <- which(drug_response$dbgap_subject_id == subject) # get indices of subject
  drugs <- as.character(drug_response$inhibitor[inds_subject]) # get drug names of subject
  drug_matrix[drugs, subject] <- drug_response$ic50[inds_subject] # add ic50 of those drugs
}

# Identifying missing values
is.na(drug_matrix) <- drug_matrix == "NA"

# Change IC50 format to numeric
for (i in 2:(dim(drug_matrix)[2] - 1)) {
  drug_matrix[, i] <- as.numeric(drug_matrix[, i])
}

# Delete patients with few drug response and input missing drug response 
drug_matrix_2<-as.matrix(t(drug_matrix)) # transpose drug 
if (sum(is.na(drug_matrix_2)) > 0) { 
  drug_matrix_3 <- drug_matrix_2[
    # NAs percentage has to be less than 20% in rows and less than 30% in columns
    -which((rowSums(is.na(drug_matrix_2)) / dim(drug_matrix_2)[2]) > 0.8), # delete those tested in a few drugs, 
    -which((colSums(is.na(drug_matrix_2)) / dim(drug_matrix_2)[1]) > 0.7)] 
  out <- impute.knn(as.matrix(drug_matrix_3), # impute missing values
                    k = 10, rowmax = 0.8,
                    colmax = 0.75, maxp = 1500, rng.seed = 362436069
  )
  drug_matrix_4 <- out$data[, -which(colnames(out$data) %in% c("Elesclomol", "JNJ-7706621"))] # Los quitamos por ser toxicos, no los queremos, cuando hicimos la normalizacion se deberian de haber ido solitos
} else {
  drug_matrix_4 <- drug_matrix_2[, -which(colnames(drug_matrix_2) %in% c("Elesclomol", "JNJ-7706621"))]
} 

# MUTATIONS MATRIX 
# Build Mutations matrix
mutations <- matrix(0, # empty matrix
                    nrow = length(subjects_ids),
                    ncol = length(unique(gene_variants$symbol))
)

rownames(mutations) <- subjects_ids
colnames(mutations) <-unique(gene_variants$symbol)

for (i in as.character(unique(gene_variants$dbgap_subject_id))) {
  jj <- unique(gene_variants$symbol[which(gene_variants$dbgap_subject_id == i)])
  mutations[i, jj] <- 1 # add 1 if there is a mutation
}

# Add Translocations 
# look for translocations:  filtering the rows to only those that contain "AML with inv" or "AML with t" 
translocations_clinical<-clinical[,c("dbgap_subject_id","specificDxAtAcquisition")][grep("AML with inv", clinical$specificDxAtAcquisition),]
translocations_clinical<-rbind(translocations_clinical,clinical[,c("dbgap_subject_id","specificDxAtAcquisition")][grep("AML with t", clinical$specificDxAtAcquisition),])

trans<-as.character(unique(translocations_clinical$specificDxAtAcquisition))
trans_matrix<-matrix(0, nrow=nrow(mutations), ncol=length(trans)) # add columns for translocations
colnames(trans_matrix)<-trans
rownames(trans_matrix)<-rownames(mutations)

for(j in colnames(trans_matrix)){
  p<-as.character(translocations_clinical$dbgap_subject_id[which(translocations_clinical$specificDxAtAcquisition ==j)])
  p<-p[which(p %in% rownames(trans_matrix))]
  trans_matrix[p,j]<-1 # add 1 if there is a translocation
}

# Change column name and also add t: extract the second element of each translocation name (i.e., the specific type of translocation) and use it as the new column name in add. 
trans<-sapply(trans, function(X){unlist(strsplit(X, split = "; "))[2]})
colnames(trans_matrix)<-trans
mutations<-cbind(mutations,trans_matrix)
mutations<-as.matrix(mutations)

# Keep only samples in the 4 matrices (gene_variants,expression,drug_response,clinical)
drug_subjects <- rownames(drug_matrix_4)
mutation_subjects <- rownames(mutations)
expression_subjects <- colnames(expression)
clinical_subjects <- clinical$dbgap_subject_id

keep_subjects <- Reduce(intersect, list(drug_subjects,mutation_subjects,expression_subjects,clinical_subjects))
drug_matrix_4 <- drug_matrix_4[keep_subjects,]
drug_matrix_2 <- drug_matrix_2[keep_subjects,colnames(drug_matrix_4)] # same drugs and patients
mutations <- mutations[keep_subjects,]
# Filter to keep only the 5000 genes with the highest variance among patients
gene_variances <- rowVars(as.matrix(expression))
top_genes <- order(gene_variances, decreasing = TRUE)[1:5000]
expression <- expression[top_genes, keep_subjects]
expression <- t(expression) # transpose expression matrix
clinical <- clinical[keep_subjects,]

table(clinical$cohort)

# Filter: keep only mutations present in at least 1% of patients 
mutations<-mutations[, colSums(mutations)>0.01*nrow(mutations)] 

# IC50* Calculation 
drug<-log10(drug_matrix_4)-1
drug <- t(t(drug) - colMeans(drug, na.rm=T))  # unify the drugs by columns, we compensate with the dosage
drug_response_matrix <- drug - min(drug) 


gamma <- 0.7
pseudop <- exp(-gamma*drug_response_matrix) # Convert weights into pseudoprobabilities for voting scheme
pseudop <- pseudop /rowSums(pseudop) # Probabilities sum up one

# Separate variables into train and test (Waves 1+2 or Both and Waves 3+4) REVISAR
# get ids of subjects
subject_id_w12 <- as.character(clinical$dbgap_subject_id[which(clinical$cohort=="Waves1+2" | clinical$cohort=="Both")])
subject_id_w34 <- as.character(clinical$dbgap_subject_id[which(clinical$cohort=="Waves3+4")])

drug_response_w12 <- drug_response_matrix[subject_id_w12,]
drug_response_w34 <- drug_response_matrix[subject_id_w34,]

pseudop_w12 <- pseudop[subject_id_w12,]
pseudop_w34 <- pseudop[subject_id_w34,]

mutations_w12 <- mutations[subject_id_w12,]
mutations_w34 <- mutations[subject_id_w34,]

expression_w12 <- expression[subject_id_w12,]
expression_w34 <- expression[subject_id_w34,]

# Tidying up the variables (delete extra variables)
rm(list = c("drug","drug_matrix_2", "drug_matrix_3", "drug_matrix_4", "drug_response", "gene_variants"))
rm(list = c("p","j","out", "translocations_clinical", "trans_matrix", "drug_matrix","i","jj")) 

# Save to use in GDSC testing
#save(drug_response_w12,pseudop_w12,expression_w12,mutations_w12, file = paste(folder_dir,"/training_data_5000.RData",sep=""))
# load(file = paste(folder_dir,"Rdata/training_data.RData",sep=""))

#############################################################
# BeatAML Waves1+2 models
#############################################################

# Comparison of methods with Cross-Validation --------------------------------------------
"
We performed a 5-fold cross-validation using the BeatAML dataset.
We trained all models with genetic variants data from 319 patients, dividing the cohort between the training samples 4-folds and testing samples the selected 1-fold.
Each of the folds were tested, and the predicted IC50* for the 5-fold testing was compared for all the methods and compared against the Oracle -the drug with the optimum IC50.
We calculated the Oracle as the minimum IC50* value for each patient.
"

set.seed(2022)
folds <- 5
groups <- sample(1:folds, nrow(drug_response_w12), replace = T)
treatment_RF_mut <- treatment_XGBoost_mut <- treatment_ODT_mut <- rep(NA, length(groups))
treatment_RF_exp <- treatment_XGBoost_exp <- treatment_ODT_exp <- rep(NA, length(groups))

treatment_RF_mut_v3 <- treatment_RF_mut_v5 <- treatment_RF_mut_v10 <- treatment_RF_mut_v20 <- rep(NA, length(groups))
treatment_XGBoost_mut_v3 <- treatment_XGBoost_mut_v5 <- treatment_XGBoost_mut_v10 <- treatment_XGBoost_mut_v20 <- rep(NA, length(groups))

treatment_RF_exp_v3 <- treatment_RF_exp_v5 <- treatment_RF_exp_v10 <- treatment_RF_exp_v20 <- rep(NA, length(groups))
treatment_XGBoost_exp_v3 <- treatment_XGBoost_exp_v5 <- treatment_XGBoost_exp_v10 <- treatment_XGBoost_exp_v20 <- rep(NA, length(groups))

# ODT parameters
minbucket_mut <- 6
minbucket_exp <- 17

# RF parameters
num_trees <- 500

# XGBoost parameters
num_rounds <- 20

# Predict functions require the drug_response_matrix just to extract the names of the drugs
# Train and test
for (group in 1:folds) {
  cat("Fold: ", group, "\n")
  dejar <- groups != group
  quitar <- !dejar
  
  # ODTMut
  cat("ODTMut\n")
  ODTMut <- trainTree(PatientData = mutations_w12[dejar,], PatientSensitivity = drug_response_w12[dejar,], minbucket = minbucket_mut)
  treatment_ODT_mut[quitar] <- predictTree(tree = ODTMut, PatientSensitivityTrain = drug_response_w12[quitar,], PatientData = mutations_w12[quitar,])
  
  # RFMut
  cat("RFMut\n")
  RFMut <- train_rf_hard(mutations_w12[dejar,], drug_response_w12[dejar,],num_trees=num_trees)
  treatment_RF_mut[quitar] <- predict_rf_hard(RFMut,mutations_w12[quitar,])
  
  # RFMut vote=3
  RFMut_v3 <- train_rf_votes(mutations_w12[dejar,], pseudop_w12[dejar,],num_trees=num_trees, num_votes = 3)
  treatment_RF_mut_v3[quitar] <- predict_rf(RFMut_v3,mutations_w12[quitar,])
  
  # RFMut vote=5
  RFMut_v5 <- train_rf_votes(mutations_w12[dejar,], pseudop_w12[dejar,],num_trees=num_trees, num_votes = 5)
  treatment_RF_mut_v5[quitar] <- predict_rf(RFMut_v5,mutations_w12[quitar,])
  
  # RFMut vote=10
  RFMut_v10 <- train_rf_votes(mutations_w12[dejar,], pseudop_w12[dejar,],num_trees=num_trees, num_votes = 10)
  treatment_RF_mut_v10[quitar] <- predict_rf(RFMut_v10,mutations_w12[quitar,])
  
  # RFMut vote=20
  RFMut_v20 <- train_rf_votes(mutations_w12[dejar,], pseudop_w12[dejar,],num_trees=num_trees, num_votes = 20)
  treatment_RF_mut_v20[quitar] <- predict_rf(RFMut_v20,mutations_w12[quitar,])
  
  # XGBoostMut
  cat("XGBoostMut\n")
  XGBoostMut <- train_xgb_hard(mutations_w12[dejar,], drug_response_w12[dejar,],num_rounds=num_rounds)
  treatment_XGBoost_mut[quitar] <- predict_xgb(XGBoostMut,mutations_w12[quitar,])

  # XGBoostMut vote=3
  cat("XGBoostMut vote=3\n")
  XGBoostMut_v3 <- train_xgb_votes(mutations_w12[dejar,], pseudop_w12[dejar,], num_rounds=num_rounds, num_votes = 3)
  treatment_XGBoost_mut_v3[quitar] <- predict_xgb(XGBoostMut_v3, mutations_w12[quitar,])
  
  # XGBoostMut vote=5
  cat("XGBoostMut vote=5\n")
  XGBoostMut_v5 <- train_xgb_votes(mutations_w12[dejar,], pseudop_w12[dejar,], num_rounds=num_rounds, num_votes = 5)
  treatment_XGBoost_mut_v5[quitar] <- predict_xgb(XGBoostMut_v5, mutations_w12[quitar,])
  
  # XGBoostMut vote=10
  cat("XGBoostMut vote=10\n")
  XGBoostMut_v10 <- train_xgb_votes(mutations_w12[dejar,], pseudop_w12[dejar,], num_rounds=num_rounds, num_votes = 10)
  treatment_XGBoost_mut_v10[quitar] <- predict_xgb(XGBoostMut_v10, mutations_w12[quitar,])
  
  # XGBoostMut vote=20
  cat("XGBoostMut vote=20\n")
  XGBoostMut_v20 <- train_xgb_votes(mutations_w12[dejar,], pseudop_w12[dejar,], num_rounds=num_rounds, num_votes = 20)
  treatment_XGBoost_mut_v20[quitar] <- predict_xgb(XGBoostMut_v20, mutations_w12[quitar,])

  # ODTExp
  cat("ODTExp\n")
  ODTExp <- trainTree(PatientData = expression_w12[dejar,], PatientSensitivity = drug_response_w12[dejar,], minbucket = minbucket_exp)
  treatment_ODT_exp[quitar] <- predictTree(tree = ODTExp, PatientSensitivityTrain = drug_response_w12[quitar,], PatientData = expression_w12[quitar,])
  
  # RFExp
  cat("RFExp\n")
  RFExp <- train_rf_hard(expression_w12[dejar,], drug_response_w12[dejar,],num_trees=num_trees)
  treatment_RF_exp[quitar] <- predict_rf_hard(RFExp,expression_w12[quitar,])

  # RFExp vote=3
  RFExp_v3 <- train_rf_votes(expression_w12[dejar,], pseudop_w12[dejar,], num_trees=num_trees, num_votes = 3)
  treatment_RF_exp_v3[quitar] <- predict_rf(RFExp_v3, expression_w12[quitar,])
  
  # RFExp vote=5
  RFExp_v5 <- train_rf_votes(expression_w12[dejar,], pseudop_w12[dejar,], num_trees=num_trees, num_votes = 5)
  treatment_RF_exp_v5[quitar] <- predict_rf(RFExp_v5, expression_w12[quitar,])
  
  # RFExp vote=10
  RFExp_v10 <- train_rf_votes(expression_w12[dejar,], pseudop_w12[dejar,], num_trees=num_trees, num_votes = 10)
  treatment_RF_exp_v10[quitar] <- predict_rf(RFExp_v10, expression_w12[quitar,])
  
  # RFExp vote=20
  RFExp_v20 <- train_rf_votes(expression_w12[dejar,], pseudop_w12[dejar,], num_trees=num_trees, num_votes = 20)
  treatment_RF_exp_v20[quitar] <- predict_rf(RFExp_v20, expression_w12[quitar,])
  
  # XGBoostExp
  cat("XGBoostExp\n")
  XGBoostExp <- train_xgb_hard(expression_w12[dejar,], drug_response_w12[dejar,],num_rounds=num_rounds)
  treatment_XGBoost_exp[quitar] <- predict_xgb(XGBoostExp,expression_w12[quitar,])

  # XGBoostExp vote=3
  cat("XGBoostExp vote=3\n")
  XGBoostExp_v3 <- train_xgb_votes(expression_w12[dejar,], pseudop_w12[dejar,], num_rounds=num_rounds, num_votes = 3)
  treatment_XGBoost_exp_v3[quitar] <- predict_xgb(XGBoostExp_v3, expression_w12[quitar,])
  
  # XGBoostExp vote=5
  cat("XGBoostExp vote=5\n")
  XGBoostExp_v5 <- train_xgb_votes(expression_w12[dejar,], pseudop_w12[dejar,], num_rounds=num_rounds, num_votes = 5)
  treatment_XGBoost_exp_v5[quitar] <- predict_xgb(XGBoostExp_v5, expression_w12[quitar,])
  
  # XGBoostExp vote=10
  cat("XGBoostExp vote=10\n")
  XGBoostExp_v10 <- train_xgb_votes(expression_w12[dejar,], pseudop_w12[dejar,], num_rounds=num_rounds, num_votes = 10)
  treatment_XGBoost_exp_v10[quitar] <- predict_xgb(XGBoostExp_v10, expression_w12[quitar,])
  
  # XGBoostExp vote=20
  cat("XGBoostExp vote=20\n")
  XGBoostExp_v20 <- train_xgb_votes(expression_w12[dejar,], pseudop_w12[dejar,], num_rounds=num_rounds, num_votes = 20)
  treatment_XGBoost_exp_v20[quitar] <- predict_xgb(XGBoostExp_v20, expression_w12[quitar,])
  
}

treatmentOracle <- apply(drug_response_w12, 1, which.min)
### Plotting BeatAML Cross-validation Results and comparing means difference
treatment_plot <- data.frame(Method = "ORACLE", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatmentOracle)])
treatment_plot <- rbind(treatment_plot, data.frame(Method = "ODTMut", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_ODT_mut)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_mut)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_mut)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "ODTExp", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_ODT_exp)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_exp)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_exp)]))

# Adding plots for RFMut with votes
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut_v3", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_mut_v3)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut_v5", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_mut_v5)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut_v10", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_mut_v10)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut_v20", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_mut_v20)]))

# Adding plots for XGBoostMut with votes
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut_v3", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_mut_v3)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut_v5", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_mut_v5)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut_v10", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_mut_v10)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut_v20", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_mut_v20)]))

# Adding plots for RFExp with votes
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp_v3", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_exp_v3)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp_v5", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_exp_v5)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp_v10", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_exp_v10)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp_v20", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_RF_exp_v20)]))

# Adding plots for XGBoostExp with votes
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp_v3", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_exp_v3)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp_v5", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_exp_v5)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp_v10", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_exp_v10)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp_v20", IC50 = drug_response_w12[cbind(1:nrow(drug_response_w12), treatment_XGBoost_exp_v20)]))

treatment_plot$Method <- factor(treatment_plot$Method, levels = c(
  "ORACLE", "ODTMut", "RFMut", "RFMut_v3", "RFMut_v5", "RFMut_v10", "RFMut_v20",
  "XGBoostMut", "XGBoostMut_v3", "XGBoostMut_v5", "XGBoostMut_v10", "XGBoostMut_v20",
  "ODTExp", "RFExp", "RFExp_v3", "RFExp_v5", "RFExp_v10", "RFExp_v20",
  "XGBoostExp", "XGBoostExp_v3", "XGBoostExp_v5", "XGBoostExp_v10", "XGBoostExp_v20"
))


colors <- c("ORACLE" = "#E36192", "ODT" = "#39A7AE", "RF" = "#96D6B6", "XGBoost" = "#FDE5B0")
# Assign colors and patterns based on model type
treatment_plot$ColorGroup <- ifelse(grepl("ORACLE", treatment_plot$Method), "ORACLE",
               ifelse(grepl("ODT", treatment_plot$Method), "ODT",
               ifelse(grepl("RF", treatment_plot$Method), "RF", "XGBoost")))

treatment_plot$PatternGroup <- ifelse(grepl("Mut", treatment_plot$Method), "Mut",
                 ifelse(grepl("Exp", treatment_plot$Method), "Exp", "ORACLE"))

treatment_plot_gg <- ggplot(treatment_plot, aes(x = Method, y = IC50, fill = ColorGroup, pattern = PatternGroup)) +
  geom_boxplot_pattern(
  lwd = 0.8, position = position_dodge(0.8),
  pattern_fill = "black", pattern_density = 0.2, pattern_spacing = 0.02
  ) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  scale_pattern_manual(values = c("ORACLE" = "none", "Mut" = "circle", "Exp" = "stripe")) +
  ylab("IC50*") +
  ggtitle("5-fold CV in BeatAML (Waves 1+2)") +
  theme(
  text = element_text(family = "Roboto"),
  axis.text.x = element_text(size = 24, angle = 80, vjust = 0.6),
  axis.text.y = element_text(size = 20),
  axis.title = element_text(size = 27),
  plot.title = element_text(size = 32, face = "bold"),
  plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
  panel.background = element_rect(fill = "transparent"),
  legend.position = "right"
  ) +
  ylim(0, 5) +
  guides(
  fill = guide_legend(title = "Model Type"),
  pattern = guide_legend(title = "Pattern Type")
  )

#treatment_plot_gg
ggsave(paste(folder_dir, "/images/boxplot_validate_5000_0_7.png", sep = ""), treatment_plot_gg, width = 12, height = 7, dpi = 500)


# Calculate mean and median IC50 for each method
summary_stats <- treatment_plot %>%
  group_by(Method) %>%
  summarise(
    Mean_IC50 = mean(IC50, na.rm = TRUE),
    Median_IC50 = median(IC50, na.rm = TRUE)
  ) %>%
  arrange(Mean_IC50)
write_xlsx(summary_stats, paste(folder_dir, "/ic50_results_validation_5000_0_7.xlsx", sep = ""))
# Training models in complete BeatAML Waves1+2 data, testing in Waves3+4 and timing
"
We trained the models in BeatAML cohort and measured the training time for each model.
Technical implementation refers to the computational burden and software that the method requires.Despite it could be considered less important, some of the algorithms require hours of computing time for the BeatAML of the subset of AML samples in GDSC -that be considered to be small/medium size.
" 
# Train with ALL the data and time the performance
timing_results <- data.frame(Method = character(), Time = character(), Variables = integer(), stringsAsFactors = FALSE)

# ODTMut
cat("ODTMut\n")
start_time <- Sys.time()
ODTMut_All <- trainTree(PatientData = mutations_w12, PatientSensitivity = drug_response_w12, minbucket = minbucket_mut)
treatmentODTMut_w34 <- predictTree(tree = ODTMut_All, PatientSensitivityTrain = drug_response_w12, PatientData = mutations_w34)
end_time <- Sys.time()
time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
timing_results <- rbind(timing_results, data.frame(Method = "ODTMut", Time = time_taken, Variables = 0))

# RFMut
cat("RFMut\n")
start_time <- Sys.time()
RFMut_All <- train_rf_hard(mutations_w12, drug_response_w12)
treatmentRFMut_w34 <- predict_rf_hard(RFMut_All, mutations_w34)
end_time <- Sys.time()
time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
num_variables <- length(RFMut_All$forest$independent.variable.names)
timing_results <- rbind(timing_results, data.frame(Method = "RFMut", Time = time_taken, Variables = num_variables))

# XGBoostMut
cat("XGBoostMut\n")
start_time <- Sys.time()
XGBoostMut_All <- train_xgb_hard(mutations_w12, drug_response_w12, num_rounds = num_rounds)
treatmentXGBoostMut_w34 <- predict_xgb(XGBoostMut_All, mutations_w34)
end_time <- Sys.time()
time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
num_variables <- length(XGBoostMut_All$feature_names)
timing_results <- rbind(timing_results, data.frame(Method = "XGBoostMut", Time = time_taken, Variables = num_variables))

# RFMut with votes
for (votes in c(3, 5, 10, 20)) {
  cat(paste("RFMut vote=", votes, "\n", sep = ""))
  start_time <- Sys.time()
  RFMut_All <- train_rf_votes(mutations_w12, pseudop_w12, num_trees = num_trees, num_votes = votes)
  assign(paste0("RFMut_All_v", votes), RFMut_All)  # Assign model to a unique variable name
  assign(paste0("treatmentRFMut_w34_v", votes), predict_rf(RFMut_All, mutations_w34))
  end_time <- Sys.time()
  time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
  num_variables <- length(RFMut_All$forest$independent.variable.names)
  timing_results <- rbind(timing_results, data.frame(Method = paste0("RFMut_v", votes), Time = time_taken, Variables = num_variables))
}

# XGBoostMut with votes
for (votes in c(3, 5, 10, 20)) {
  cat(paste("XGBoostMut vote=", votes, "\n", sep = ""))
  start_time <- Sys.time()
  XGBoostMut_All <- train_xgb_votes(mutations_w12, pseudop_w12, num_rounds = num_rounds, num_votes = votes)
  assign(paste0("XGBoostMut_All_v", votes), XGBoostMut_All)  # Assign model to a unique variable name
  assign(paste0("treatmentXGBoostMut_w34_v", votes), predict_xgb(XGBoostMut_All, mutations_w34))
  end_time <- Sys.time()
  time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
  num_variables <- length(XGBoostMut_All$feature_names)
  timing_results <- rbind(timing_results, data.frame(Method = paste0("XGBoostMut_v", votes), Time = time_taken, Variables = num_variables))
}

# ODTExp
cat("ODTExp\n")
start_time <- Sys.time()
ODTExp_All <- trainTree(PatientData = expression_w12, PatientSensitivity = drug_response_w12, minbucket = minbucket_exp)
treatmentODTExp_w34 <- predictTree(tree = ODTExp_All, PatientSensitivityTrain = drug_response_w12, PatientData = expression_w34)
end_time <- Sys.time()
time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
timing_results <- rbind(timing_results, data.frame(Method = "ODTExp", Time = time_taken, Variables = 0))

# RFExp
cat("RFExp\n")
start_time <- Sys.time()
RFExp_All <- train_rf_hard(expression_w12, drug_response_w12)
treatmentRFExp_w34 <- predict_rf_hard(RFExp_All, expression_w34)
end_time <- Sys.time()
time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
num_variables <- length(RFExp_All$forest$independent.variable.names)
timing_results <- rbind(timing_results, data.frame(Method = "RFExp", Time = time_taken, Variables = num_variables))

# RFExp with votes
for (votes in c(3, 5, 10, 20)) {
  cat(paste("RFExp vote=", votes, "\n", sep = ""))
  start_time <- Sys.time()
  RFExp_All <- train_rf_votes(expression_w12, pseudop_w12, num_trees = num_trees, num_votes = votes)
  assign(paste0("RFExp_All_v", votes), RFExp_All)  # Assign model to a unique variable name
  assign(paste0("treatmentRFExp_w34_v", votes), predict_rf(RFExp_All, expression_w34))
  end_time <- Sys.time()
  time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
  num_variables <- length(RFExp_All$forest$independent.variable.names)
  timing_results <- rbind(timing_results, data.frame(Method = paste0("RFExp_v", votes), Time = time_taken, Variables = num_variables))
}

# XGBoostExp
cat("XGBoostExp\n")
start_time <- Sys.time()
XGBoostExp_All <- train_xgb_hard(expression_w12, drug_response_w12, num_rounds = num_rounds)
treatmentXGBoostExp_w34 <- predict_xgb(XGBoostExp_All, expression_w34)
end_time <- Sys.time()
time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
num_variables <- length(XGBoostExp_All$feature_names)
timing_results <- rbind(timing_results, data.frame(Method = "XGBoostExp", Time = time_taken, Variables = num_variables))

# XGBoostExp with votes
for (votes in c(3, 5, 10, 20)) {
  cat(paste("XGBoostExp vote=", votes, "\n", sep = ""))
  start_time <- Sys.time()
  XGBoostExp_All <- train_xgb_votes(expression_w12, pseudop_w12, num_rounds = num_rounds, num_votes = votes)
  assign(paste0("XGBoostExp_All_v", votes), XGBoostExp_All)  # Assign model to a unique variable name
  assign(paste0("treatmentXGBoostExp_w34_v", votes), predict_xgb(XGBoostExp_All, expression_w34))
  end_time <- Sys.time()
  time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
  num_variables <- length(XGBoostExp_All$feature_names)
  timing_results <- rbind(timing_results, data.frame(Method = paste0("XGBoostExp_v", votes), Time = time_taken, Variables = num_variables))
}


# Save timing results to an Excel file
write_xlsx(timing_results, paste(folder_dir, "/timing_results_5000.xlsx", sep = ""))

treatmentOracle_w34 <- apply(drug_response_w34, 1, which.min)

### Plotting BeatAML Cross-validation Results and comparing means difference
treatment_plot_test<-data.frame(Method="ORACLE", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34),treatmentOracle_w34)])
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="ODTMut", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTMut_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFMut", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFMut_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFMut_v3", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFMut_w34_v3)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFMut_v5", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFMut_w34_v5)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFMut_v10", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFMut_w34_v10)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFMut_v20", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFMut_w34_v20)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostMut", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostMut_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostMut_v3", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostMut_w34_v3)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostMut_v5", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostMut_w34_v5)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostMut_v10", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostMut_w34_v10)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostMut_v20", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostMut_w34_v20)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="ODTExp", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentODTExp_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFExp", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFExp_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFExp_v3", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFExp_w34_v3)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFExp_v5", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFExp_w34_v5)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFExp_v10", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFExp_w34_v10)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="RFExp_v20", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentRFExp_w34_v20)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostExp", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostExp_w34)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostExp_v3", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostExp_w34_v3)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostExp_v5", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostExp_w34_v5)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostExp_v10", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostExp_w34_v10)]))
treatment_plot_test<-rbind(treatment_plot_test,data.frame(Method="XGBoostExp_v20", IC50=drug_response_w34[cbind(1:nrow(drug_response_w34), treatmentXGBoostExp_w34_v20)]))

treatment_plot_test$Method <- factor(treatment_plot_test$Method, levels = c(
  "ORACLE", "ODTMut", "RFMut", "RFMut_v3", "RFMut_v5", "RFMut_v10", "RFMut_v20",
  "XGBoostMut", "XGBoostMut_v3", "XGBoostMut_v5", "XGBoostMut_v10", "XGBoostMut_v20",
  "ODTExp", "RFExp", "RFExp_v3", "RFExp_v5", "RFExp_v10", "RFExp_v20",
  "XGBoostExp", "XGBoostExp_v3", "XGBoostExp_v5", "XGBoostExp_v10", "XGBoostExp_v20"
))

colors <- c("ORACLE" = "#E36192", "ODT" = "#39A7AE", "RF" = "#96D6B6", "XGBoost" = "#FDE5B0")
# Assign colors and patterns based on model type
treatment_plot_test$ColorGroup <- ifelse(grepl("ORACLE", treatment_plot_test$Method), "ORACLE",
                                  ifelse(grepl("ODT", treatment_plot_test$Method), "ODT",
                                  ifelse(grepl("RF", treatment_plot_test$Method), "RF", "XGBoost")))

treatment_plot_test$PatternGroup <- ifelse(grepl("Mut", treatment_plot_test$Method), "Mut",
                                    ifelse(grepl("Exp", treatment_plot_test$Method), "Exp", "ORACLE"))

treatment_plot_test_gg <- ggplot(treatment_plot_test, aes(x = Method, y = IC50, fill = ColorGroup, pattern = PatternGroup)) +
  geom_boxplot_pattern(
    lwd = 0.8, position = position_dodge(0.8),
    pattern_fill = "black", pattern_density = 0.2, pattern_spacing = 0.02
  ) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  scale_pattern_manual(values = c("ORACLE" = "none", "Mut" = "circle", "Exp" = "stripe")) +
  ylab("IC50*") +
  ggtitle("Test in BeatAML (Waves 3+4)") +
  theme(
    text = element_text(family = "Roboto"),
    axis.text.x = element_text(size = 24, angle = 80, vjust = 0.6),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 27),
    plot.title = element_text(size = 32, face = "bold"),
    plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = "right"
  )+
  ylim(0, 5) +
  guides(
    fill = guide_legend(title = "Model Type"),
    pattern = guide_legend(title = "Pattern Type")
  )


ggsave(paste(folder_dir, "/images/boxplot_test_5000_0_7.png", sep = ""), treatment_plot_test_gg, width = 12, height = 7, dpi = 500)

# Calculate mean and median IC50 for each method
summary_stats <- treatment_plot_test %>%
  group_by(Method) %>%
  summarise(
    Mean_IC50 = mean(IC50, na.rm = TRUE),
    Median_IC50 = median(IC50, na.rm = TRUE)
  ) %>%
  arrange(Mean_IC50)
write_xlsx(summary_stats, paste(folder_dir, "/ic50_results_test_5000_0_7.xlsx", sep = ""))

# Save results for GDSC
save(ODTMut_All, RFMut_All, RFMut_All_v3, RFMut_All_v5, RFMut_All_v10, RFMut_All_v20,
  XGBoostMut_All, XGBoostMut_All_v3, XGBoostMut_All_v5, XGBoostMut_All_v10, XGBoostMut_All_v20,
  drug_response_w12,
  file = paste(folder_dir,"/Rdata/models_all_5000_0_7.RData",sep=""))

# # Save results
# save(ODTMut_All, RFMut_All, 
#      XGBoostMut_All, 
#      ODTExp_All, RFExp_All, 
#      XGBoostExp_All, drug_response_w12,
#      file = paste(folder_dir,"/Rdata/models_all_5000_0_7.RData",sep=""))

# save(treatmentODTMut_w34, treatmentRFMut_w34, treatmentRFMut_w34_v3, treatmentRFMut_w34_v5, treatmentRFMut_w34_v10, treatmentRFMut_w34_v20, 
#   treatmentXGBoostMut_w34, treatmentXGBoostMut_w34_v3, treatmentXGBoostMut_w34_v5, treatmentXGBoostMut_w34_v10, treatmentXGBoostMut_w34_v20, 
#   treatmentODTExp_w34, treatmentRFExp_w34, treatmentRFExp_w34_v3, treatmentRFExp_w34_v5, treatmentRFExp_w34_v10, treatmentRFExp_w34_v20, 
#   treatmentXGBoostExp_w34, treatmentXGBoostExp_w34_v3, treatmentXGBoostExp_w34_v5, treatmentXGBoostExp_w34_v10, treatmentXGBoostExp_w34_v20, 
#   file = paste(folder_dir,"/Rdata/results_test_5000_0_7.RData",sep=""))


# Load results
# load(file = paste(folder_dir,"Rdata/mut_models_all.RData",sep=""))
# load(file = paste(folder_dir,"Rdata/mut_results_test.RData",sep=""))  

#

# Comparison of explainability --------------------------------------------
"
Feature importance and graphic
"
# ODT
niceTree(ODTMut_All)
niceTree(ODTExp_All)

# RF - RFExp_v10
RFExp_All <- train_rf_votes(expression_w12, pseudop_w12, num_trees = num_trees, num_votes = 10)
plot(RFExp_All)
varImpPlot(RFExp_All)

# XGBoost - XGBoostExp_v5
# Extract variable importance
XGBoostExp_All <- train_xgb_votes(expression_w12, pseudop_w12, num_rounds = num_rounds, num_votes = 5)
importance_matrix <- xgb.importance(feature_names = colnames(expression_w34), model = XGBoostExp_All)
# Print the importance matrix
print(importance_matrix)
# Plot the importance
xgb.plot.importance(importance_matrix = importance_matrix)

# Plot feature importance -----

# Extract variable importance for Random Forest
var_importance_rf <- importance(RFExp_All)
var_importance_rf_df <- data.frame(
  Genes = rownames(var_importance_rf),
  Importance = var_importance_rf[, "MeanDecreaseGini"],
  Model = "RFExp_v10",
  stringsAsFactors = FALSE
)
var_importance_rf_df <- var_importance_rf_df[order(-var_importance_rf_df$Importance), ][1:20, ]
var_importance_rf_df$Genes <- factor(var_importance_rf_df$Genes, levels = var_importance_rf_df$Genes)

# Extract variable importance for XGBoost
importance_matrix_xgb <- xgb.importance(feature_names = colnames(expression_w34), model = XGBoostExp_All)
importance_df_xgb <- as.data.frame(importance_matrix_xgb)
importance_df_xgb <- importance_df_xgb[order(-importance_df_xgb$Gain), ][1:20, ]
importance_df_xgb <- importance_df_xgb %>%
  rename(Genes = Feature, Importance = Gain) %>%
  mutate(Model = "XGBoostExp_v5") %>%
  select(Genes, Importance, Model)  # Select only the necessary columns

# Combine the data frames
combined_data <- rbind(var_importance_rf_df, importance_df_xgb)

# Define custom colors for each model
model_colors <- c("RFExp_v10" = "#96D6B6", "XGBoostExp_v5" = "#FDE5B0")

# Create the combined plot with free scales and model-specific colors
combined_plot <- ggplot(combined_data, aes(x = Genes, y = Importance, fill = Model)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8) +
  facet_wrap(~ Model, scales = "free", nrow = 1) +  # Free scales for each facet
  theme_bw() +
  scale_fill_manual(values = model_colors) +  # Use custom colors for each model
  ylab("Importance") +
  ggtitle("Variable Importance (Top 20 Genes)") +
  theme(
    text = element_text(family = "Roboto"),
    axis.text.x = element_text(size = 20, angle = 80, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 35),
    plot.title = element_text(size = 35, face = "bold"),
    plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = "none",  # No legend needed for variable names
    strip.text = element_text(size = 30, face = "bold")  # Customize facet labels
  ) +
  guides(
    fill = guide_legend(title = "Model")
  )

# Display the plot
print(combined_plot)

# Save the plot
ggsave(paste(folder_dir, "/images/combined_importance_plot.png", sep = ""), combined_plot, width = 14, height = 7, dpi = 500)
# Comparison of time and variables --------------------------------------------
"
Implementability: Time and variables
"
# Load required libraries
library(ggplot2)
library(ggpattern)
library(tidyr)  # For reshaping data

# Create a data frame with the provided data
data <- data.frame(
  Method = c("ODTMut", "RFMut", "XGBoostMut", "ODTExp", "RFExp", "XGBoostExp"),
  Time = c(0.008, 0.329, 0.031, 5.002, 18.749, 3.382),
  Variables = c(11, 70, 70, 5, 5000, 5000)
)

# Reshape the data for faceting
data_long <- pivot_longer(data, cols = c(Time, Variables), names_to = "Metric", values_to = "Value")

# Define colors for each method
colors <- c("ODT" = "#39A7AE", "RF" = "#96D6B6", "XGBoost" = "#FDE5B0")

# Assign colors and patterns based on method
data_long$ColorGroup <- ifelse(grepl("ODT", data_long$Method), "ODT",
                               ifelse(grepl("RF", data_long$Method), "RF", "XGBoost"))

data_long$PatternGroup <- ifelse(grepl("Mut", data_long$Method), "Mut",
                                 ifelse(grepl("Exp", data_long$Method), "Exp", "none"))

# Create custom labels for facets
custom_labels <- c("Time" = "Training time (sec)", "Variables" = "Number of Features")

# Sort the data within each Metric group and convert Method to factor
data_long <- data_long %>%
  group_by(Metric) %>%
  arrange(desc(Value)) %>%
  mutate(Method = factor(Method, levels = unique(Method))) %>%
  ungroup()

# Createthe combined plot with custom facet labels and horizontal bars
combined_plot <- ggplot(data_long, aes(x = Value, y = Method, fill = ColorGroup, pattern = PatternGroup)) +
  geom_col_pattern(
    position = position_dodge2(preserve = "single"),
    pattern_fill = "black", pattern_density = 0.1, pattern_spacing = 0.05,
    width = 0.8
  ) +
  geom_text(aes(label = round(Value, 2)), hjust = 0, size=9, color = "black",face = "bold") +  # Add values to the right of each bar
  facet_wrap(~ Metric, scales = "free_x", ncol = 1, labeller = as_labeller(custom_labels)) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  scale_pattern_manual(values = c("Mut" = "circle", "Exp" = "stripe", "none" = "none")) +
  xlab("Value") +
  ylab("Method") +
  ggtitle("Comparison of Time and Number of Features") +
  theme(
    text = element_text(family = "Roboto"),
    axis.text.y = element_text(size = 32),
    axis.text.x = element_text(size = 31),
    axis.title = element_text(size = 34),
    plot.title = element_text(size = 36, face = "bold"),
    plot.background = element_rect(fill = "transparent", size = 0, color = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    legend.position = "right",
    strip.text = element_text(size = 32, face = "bold"),  # Customize facet labels
    legend.text = element_text(size = 30),  # Increase legend text size
    legend.title = element_text(size = 30),  # Increase legend title size
    legend.key.size = unit(30, "points")  # Increase legend key size
  ) +
  guides(
    fill = guide_legend(title = "Model Type"),
    #pattern = guide_legend(title = "Data Type")
  )
# Display the plot
#print(combined_plot)

# Save the plot
ggsave(paste(folder_dir, "/images/combined_implementability_plot.png", sep = ""), combined_plot, width = 25, height = 8, dpi = 500)

