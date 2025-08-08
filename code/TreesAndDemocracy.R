################################################################################
# Trees & Democracy: Drug Response Prediction Analysis
# 
# This script analyzes drug response prediction using different machine learning
# methods on the BeatAML dataset.
################################################################################

#===============================================================================
# 1. Setup and Configuration

# Directory Setup
mac <- "/Users/katyna/Library/CloudStorage/OneDrive----/"
angel <- "C:/Users/arubio/OneDrive - Tecnun/"
computer <- mac
setwd(computer)

# Output Directory
folder_dir <- "TreesAndDemocracy/output"

#===============================================================================
# 2. Load Required Libraries
library(readxl)
library(matrixStats)
library(partykit)
library(impute)
library(ggplot2)
library(ggpattern)
library(writexl)
library(pheatmap)
library(dplyr)
library(ODT)

# add option to load to macos fonts
if (Sys.info()[1] == "Roboto") {
  loadfonts(device = "macOS")
}
#windowsFonts("Roboto" = windowsFont("Roboto"))

source("./GithubProjects/Trees&Democracy/code/XAIfunctions.R")
# Best parameters from HP tunning

# Import parameters from parameters.R
source("./GithubProjects/Trees&Democracy/code/parameters.R")  # Adjust the path as necessary

#===============================================================================
# 3. Data Loading and Pre-processing 

## 3.1 Clinical Data Loading ---
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

# Fill ageAtDiagnosis with ageAtSpecimenAcquisition where missing
clinical$ageAtDiagnosis[is.na(clinical$ageAtDiagnosis)] <- 
  clinical$ageAtSpecimenAcquisition[is.na(clinical$ageAtDiagnosis)]

# Calculate percentage of NAs in each column
na_percentages <- colMeans(is.na(clinical)) * 100

# Create summary dataframe of columns with less than 10% NAs
complete_cols <- data.frame(
  Column = names(na_percentages),
  NA_Percentage = na_percentages,
  Unique_Values = sapply(clinical, function(x) length(unique(na.omit(x))))
) %>%
  filter(NA_Percentage < 10) %>%
  arrange(NA_Percentage)

# Print summary
print(complete_cols, row.names = FALSE)

clinical_wanted <- clinical[,c("dbgap_subject_id","ageAtDiagnosis","consensus_sex")]
clinical_wanted <- clinical_wanted[complete.cases(clinical_wanted),]

## 3.2 Gene Expression Data Loading --- 
expression <- read_excel("GithubProjects/interpretable-ai-review/data/input/BeatAML2/beataml_waves1to4_norm_exp_dbgap.xlsx")
expression <- as.data.frame(expression)
# set symbols as row names 
rownames(expression) <- expression$display_label 
# set dbgap_subject_id as colnames
rna_ids <- clinical$dbgap_rnaseq_sample[clinical$dbgap_rnaseq_sample  %in% colnames(expression)]
expression <- expression[,rna_ids]
subject_ids <-na.omit(clinical$dbgap_subject_id[match(colnames(expression),clinical$dbgap_rnaseq_sample)])
colnames(expression) <- subject_ids
# Check distribution of expression values for a few samples
boxplot(t(expression)[1:10,1:100], main="Boxplot of Expression Values", ylab="Expression", col="lightblue")

# expression <- log2(expression + 1) # Add 1 to avoid log(0) issues - Already transformed
expression_scaled <- t(scale(t(expression), center = TRUE, scale = TRUE)) # Center and scale the data by genes
boxplot(t(expression_scaled)[1:10,1:100], main="Boxplot of Expression Values", ylab="Expression", col="lightblue")

## 3.3 Mutations Data Loading ---
gene_variants <- read_excel("GithubProjects/interpretable-ai-review/data/input/BeatAML2/beataml_wes_wv1to4_mutations_dbgap.xlsx")
colnames(gene_variants)[1]="dbgap_dnaseq_sample"
# add dbgap_subject_id and keep only samples in clinical (have also dna)
gene_variants <-merge(gene_variants,clinical[,c("dbgap_subject_id","dbgap_dnaseq_sample")], by="dbgap_dnaseq_sample")

## 3.4 Drug Response Data Loading ---
drug_response <- read_excel("GithubProjects/interpretable-ai-review/data/input/BeatAML2/beataml_probit_curve_fits_v4_dbgap.xlsx")
drug_response$dbgap_subject_id <- as.character(drug_response$dbgap_subject_id)

#===============================================================================
# 4. Data Matrix Creation and Processing

# Now it is time to pre-process the data in order to obtain two different matrices, a matrix containing all patients and their mutations available and a second matrix containing all patients and their sensitivity to the different drugs.
# For doing so, we will impute the missing values in the drug matrix.

## 4.1 Drug Response Matrix Creation and Imputation --- 
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
  drug_matrix[drugs, subject] <- drug_response$ic50[inds_subject] # add ic50 of those drugs # THEY ALSO HAVE AUC
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

## 4.2 Mutations Matrix Creation ---
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

## 4.4 Data Filtering and Final Processing ---

# Keep only samples in the 4 matrices (gene_variants,expression,drug_response,clinical)
drug_subjects <- rownames(drug_matrix_4)
mutation_subjects <- rownames(mutations)
expression_subjects <- colnames(expression)
clinical_subjects <- clinical$dbgap_subject_id
clinical_wanted_subjects <- clinical_wanted$dbgap_subject_id

keep_subjects <- Reduce(intersect, list(drug_subjects,mutation_subjects,expression_subjects,clinical_subjects))
keep_subjects2 <- Reduce(intersect, list(drug_subjects,mutation_subjects,expression_subjects,clinical_subjects,clinical_wanted_subjects))

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

# Compute IC50*
drug<-log10(drug_matrix_4)-1
drug <- t(t(drug) - colMeans(drug, na.rm=T))  # unify the drugs by columns, we compensate with the dosage
drug_response_matrix <- drug - min(drug) 

# Compute pseudoprobabilities
gamma <- 1.6 # 0.7 1 5 10 20 try more 
pseudop <- exp(-gamma*drug_response_matrix) # Convert weights into pseudoprobabilities for SEATS
pseudop <- pseudop /rowSums(pseudop) # Probabilities sum up one

gamma <- 10 # 0.7 1 5 10 20 try more 
pseudop_g10 <- exp(-gamma*drug_response_matrix) # Convert weights into pseudoprobabilities for SEATS
pseudop_g10 <- pseudop_g10 /rowSums(pseudop_g10) # Probabilities sum up one


gamma <- 20 # 0.7 1 5 10 20 try more 
pseudop_g20 <- exp(-gamma*drug_response_matrix) # Convert weights into pseudoprobabilities for SEATS
pseudop_g20 <- pseudop_g20 /rowSums(pseudop_g20) # Probabilities sum up one


# Separate variables into train and test (Waves 1+2 or Both and Waves 3+4) 
# get ids of subjects
subject_id_w12 <- as.character(clinical$dbgap_subject_id[which(clinical$cohort=="Waves1+2" | clinical$cohort=="Both")])
subject_id_w34 <- as.character(clinical$dbgap_subject_id[which(clinical$cohort=="Waves3+4")])

# We are only going to compare with pseudoprobabilities, so we do not need the drug response matrix
drug_response_w12 <- drug_response_matrix[subject_id_w12,]
drug_response_w34 <- drug_response_matrix[subject_id_w34,]

pseudop_w12 <- pseudop[subject_id_w12,]
pseudop_w34 <- pseudop[subject_id_w34,]

pseudop_w12_g10 <- pseudop_g10[subject_id_w12,]
pseudop_w34_g10 <- pseudop_g10[subject_id_w34,]

pseudop_w12_g20 <- pseudop_g20[subject_id_w12,]
pseudop_w34_g20 <- pseudop_g20[subject_id_w34,]

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

#===============================================================================
# 5. Cross-Validation Analysis

"
BeatAML Waves1+2 models

We performed a 5-fold cross-validation using the BeatAML dataset.
We trained all models with genetic variants data from 257 patients, dividing the cohort between the training samples 4-folds and testing samples the selected 1-fold.
Each of the folds were tested, and the predicted IC50* for the 5-fold testing was compared for all the methods and compared against the Oracle -the drug with the optimum IC50.
We calculated the Oracle as the minimum IC50* value for each patient.
"

## 5.1 Cross-Validation Setup --- 
set.seed(2022)
folds <- 5

# Cell-blind splitting - keep the same patients in the same fold
groups <- sample(1:folds, nrow(drug_response_w12), replace = T)

treatment_RF_mut <- treatment_XGBoost_mut <- treatment_ODT_mut <- rep(NA, length(groups))
treatment_RF_exp <- treatment_XGBoost_exp <- treatment_ODT_exp <- rep(NA, length(groups))

treatment_RF_mut_s1 <- treatment_RF_mut_s6 <- treatment_RF_mut_s16 <- treatment_RF_mut_s19 <- treatment_RF_mut_s22 <- rep(NA, length(groups))
treatment_XGBoost_mut_s1 <- treatment_XGBoost_mut_s6 <- treatment_XGBoost_mut_s16 <- treatment_XGBoost_mut_s19 <- treatment_XGBoost_mut_s22 <- rep(NA, length(groups))

treatment_RF_exp_s1 <- treatment_RF_exp_s6 <- treatment_RF_exp_s16 <- treatment_RF_exp_s19 <- treatment_RF_exp_s22 <- rep(NA, length(groups))
treatment_XGBoost_exp_s1 <- treatment_XGBoost_exp_s6 <- treatment_XGBoost_exp_s16 <- treatment_XGBoost_exp_s19 <- treatment_XGBoost_exp_s22 <- rep(NA, length(groups))

# ODT parameters
minbucket_mut <- 6
minbucket_exp <- 17

# # RF parameters
# num_trees <- 500
# 
# # XGBoost parameters
# nrounds <- 100

## 5.2 Model Training and Testing in CV ---

# Predict functions require the drug_response_matrix just to extract the names of the drugs
# Train and test
for (group in 1:folds) { #folds
  cat("Fold: ", group, "\n")
  keep_patients <- groups != group
  remove_patients <- !keep_patients
  
  # ODTMut
  cat("ODTMut\n")
  ODTMut <- trainTree(PatientData = mutations_w12[keep_patients,], PatientSensitivity = drug_response_w12[keep_patients,], minbucket = minbucket_mut)
  treatment_ODT_mut[remove_patients] <- predictTree(tree = ODTMut, PatientSensitivityTrain = drug_response_w12[keep_patients,], PatientData = mutations_w12[remove_patients,])
  
  # # RFMut
  # cat("RFMut\n")
  # RFMut <- train_rf_hard(mutations_w12[keep_patients,], drug_response_w12[keep_patients,],num_trees=num_trees)
  # treatment_RF_mut[remove_patients] <- predict_rf(RFMut,mutations_w12[remove_patients,])
  
  
  
  # RFMut seats=1
  RFMut_s1 <- train_rf_seats(mutations_w12[keep_patients,],  pseudop_w12[keep_patients,], num_trees = params_err_rfmut_seats$num_trees, mtry = params_err_rfmut_seats$mtry, min_node_size = params_err_rfmut_seats$min_node_size, max_depth = params_err_rfmut_seats$max_depth, sample_fraction = params_err_rfmut_seats$sample_fraction, num_seats = 1)
  treatment_RF_mut_s1[remove_patients] <- predict_rf(RFMut_s1, mutations_w12[remove_patients,])
  
  # RFMut seats=6
  RFMut_s6 <- train_rf_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], num_trees = params_err_rfmut_seats$num_trees, mtry = params_err_rfmut_seats$mtry, min_node_size = params_err_rfmut_seats$min_node_size, max_depth = params_err_rfmut_seats$max_depth, sample_fraction = params_err_rfmut_seats$sample_fraction, num_seats = 6)
  treatment_RF_mut_s6[remove_patients] <- predict_rf(RFMut_s6, mutations_w12[remove_patients,])
  
  # RFMut seats=16
  RFMut_s16 <- train_rf_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], num_trees = params_err_rfmut_seats$num_trees, mtry = params_err_rfmut_seats$mtry, min_node_size = params_err_rfmut_seats$min_node_size, max_depth = params_err_rfmut_seats$max_depth, sample_fraction = params_err_rfmut_seats$sample_fraction, num_seats = 16)
  treatment_RF_mut_s16[remove_patients] <- predict_rf(RFMut_s16, mutations_w12[remove_patients,])
  
  # RFMut seats=19
  RFMut_s19 <- train_rf_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], num_trees = params_err_rfmut_seats$num_trees, mtry = params_err_rfmut_seats$mtry, min_node_size = params_err_rfmut_seats$min_node_size, max_depth = params_err_rfmut_seats$max_depth, sample_fraction = params_err_rfmut_seats$sample_fraction, num_seats = 19)
  treatment_RF_mut_s19[remove_patients] <- predict_rf(RFMut_s19, mutations_w12[remove_patients,])
  
  # RFMut seats=22
  RFMut_s22 <- train_rf_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], num_trees = params_err_rfmut_seats$num_trees, mtry = params_err_rfmut_seats$mtry, min_node_size = params_err_rfmut_seats$min_node_size, max_depth = params_err_rfmut_seats$max_depth, sample_fraction = params_err_rfmut_seats$sample_fraction, num_seats = 22)
  treatment_RF_mut_s22[remove_patients] <- predict_rf(RFMut_s22, mutations_w12[remove_patients,])
  
  
  # XGBoostMut
  # cat("XGBoostMut\n")
  # XGBoostMut <- train_xgb_hard(mutations_w12[keep_patients,], drug_response_w12[keep_patients,],nrounds=nrounds)
  # treatment_XGBoost_mut[remove_patients] <- predict_xgb(XGBoostMut,mutations_w12[remove_patients,])
  
  train_xgb_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], max_depth = params_err_xgboostmut$max_depth, eta = params_err_xgboostmut$eta, gamma = params_err_xgboostmut$gamma, min_child_weight = params_err_xgboostmut$min_child_weight, subsample = params_err_xgboostmut$subsample, colsample_bytree = params_err_xgboostmut$colsample_bytree, nrounds = params_err_xgboostmut$nrounds, early_stopping_rounds = 10, num_seats = 1)
  
  # XGBoostMut seats=1
  cat("XGBoostMut seats=1\n")
  XGBoostMut_s1 <- train_xgb_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], max_depth = params_err_xgboostmut$max_depth, eta = params_err_xgboostmut$eta, gamma = params_err_xgboostmut$gamma, min_child_weight = params_err_xgboostmut$min_child_weight, subsample = params_err_xgboostmut$subsample, colsample_bytree = params_err_xgboostmut$colsample_bytree, nrounds = params_err_xgboostmut$nrounds, early_stopping_rounds = 10, num_seats = 1)
  treatment_XGBoost_mut_s1[remove_patients] <- predict_xgb(XGBoostMut_s1, mutations_w12[remove_patients,])
  
  # XGBoostMut seats=6
  cat("XGBoostMut seats=6\n")
  XGBoostMut_s6 <- train_xgb_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], max_depth = params_err_xgboostmut$max_depth, eta = params_err_xgboostmut$eta, gamma = params_err_xgboostmut$gamma, min_child_weight = params_err_xgboostmut$min_child_weight, subsample = params_err_xgboostmut$subsample, colsample_bytree = params_err_xgboostmut$colsample_bytree, nrounds = params_err_xgboostmut$nrounds, early_stopping_rounds = 10, num_seats = 6)
  treatment_XGBoost_mut_s6[remove_patients] <- predict_xgb(XGBoostMut_s6, mutations_w12[remove_patients,])
  
  # XGBoostMut seats=16
  cat("XGBoostMut seats=16\n")
  XGBoostMut_s16 <- train_xgb_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], max_depth = params_err_xgboostmut$max_depth, eta = params_err_xgboostmut$eta, gamma = params_err_xgboostmut$gamma, min_child_weight = params_err_xgboostmut$min_child_weight, subsample = params_err_xgboostmut$subsample, colsample_bytree = params_err_xgboostmut$colsample_bytree, nrounds = params_err_xgboostmut$nrounds, early_stopping_rounds = 10, num_seats = 16)
  treatment_XGBoost_mut_s16[remove_patients] <- predict_xgb(XGBoostMut_s16, mutations_w12[remove_patients,])
  
  # XGBoostMut seats=19
  cat("XGBoostMut seats=19\n")
  XGBoostMut_s19 <- train_xgb_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], max_depth = params_err_xgboostmut$max_depth, eta = params_err_xgboostmut$eta, gamma = params_err_xgboostmut$gamma, min_child_weight = params_err_xgboostmut$min_child_weight, subsample = params_err_xgboostmut$subsample, colsample_bytree = params_err_xgboostmut$colsample_bytree, nrounds = params_err_xgboostmut$nrounds, early_stopping_rounds = 10, num_seats = 19)
  treatment_XGBoost_mut_s19[remove_patients] <- predict_xgb(XGBoostMut_s19, mutations_w12[remove_patients,])
  
  # XGBoostMut seats=22
  cat("XGBoostMut seats=22\n")
  XGBoostMut_s22 <- train_xgb_seats(mutations_w12[keep_patients,], pseudop_w12[keep_patients,], max_depth = params_err_xgboostmut$max_depth, eta = params_err_xgboostmut$eta, gamma = params_err_xgboostmut$gamma, min_child_weight = params_err_xgboostmut$min_child_weight, subsample = params_err_xgboostmut$subsample, colsample_bytree = params_err_xgboostmut$colsample_bytree, nrounds = params_err_xgboostmut$nrounds, early_stopping_rounds = 10, num_seats = 22)
  treatment_XGBoost_mut_s22[remove_patients] <- predict_xgb(XGBoostMut_s22, mutations_w12[remove_patients,])
  
  # ODTExp NEED TO REDO THIS
  cat("ODTExp\n")
  ODTExp <- trainTree(PatientData = expression_w12[keep_patients,], PatientSensitivity = drug_response_w12[keep_patients,], minbucket = minbucket_exp)
  treatment_ODT_exp[remove_patients] <- predictTree(tree = ODTExp, PatientSensitivityTrain = drug_response_w12[keep_patients,], PatientData = expression_w12[remove_patients,])
  
  # RFExp
  # cat("RFExp\n")
  # RFExp <- train_rf_hard(expression_w12[keep_patients,], drug_response_w12[keep_patients,],num_trees=num_trees)
  # treatment_RF_exp[remove_patients] <- predict_rf(RFExp,expression_w12[remove_patients,])
  
  # RFExp seats=1
  RFExp_s1 <- train_rf_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,], 
                             num_trees = params_err_rfexp_seats$num_trees, 
                             mtry = params_err_rfexp_seats$mtry, 
                             min_node_size = params_err_rfexp_seats$min_node_size, 
                             max_depth = params_err_rfexp_seats$max_depth, 
                             sample_fraction = params_err_rfexp_seats$sample_fraction, 
                             num_seats = 1)
  treatment_RF_exp_s1[remove_patients] <- predict_rf(RFExp_s1, expression_w12[remove_patients,])
  
  # RFExp seats=6
  RFExp_s6 <- train_rf_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,], 
                             num_trees = params_err_rfexp_seats$num_trees, 
                             mtry = params_err_rfexp_seats$mtry, 
                             min_node_size = params_err_rfexp_seats$min_node_size, 
                             max_depth = params_err_rfexp_seats$max_depth, 
                             sample_fraction = params_err_rfexp_seats$sample_fraction, 
                             num_seats = 6)
  treatment_RF_exp_s6[remove_patients] <- predict_rf(RFExp_s6, expression_w12[remove_patients,])
  
  # RFExp seats=16
  RFExp_s16 <- train_rf_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,], 
                              num_trees = params_err_rfexp_seats$num_trees, 
                              mtry = params_err_rfexp_seats$mtry, 
                              min_node_size = params_err_rfexp_seats$min_node_size, 
                              max_depth = params_err_rfexp_seats$max_depth, 
                              sample_fraction = params_err_rfexp_seats$sample_fraction, 
                              num_seats = 16)
  treatment_RF_exp_s16[remove_patients] <- predict_rf(RFExp_s16, expression_w12[remove_patients,])
  
  # RFExp seats=19
  RFExp_s19 <- train_rf_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,], 
                              num_trees = params_err_rfexp_seats$num_trees, 
                              mtry = params_err_rfexp_seats$mtry, 
                              min_node_size = params_err_rfexp_seats$min_node_size, 
                              max_depth = params_err_rfexp_seats$max_depth, 
                              sample_fraction = params_err_rfexp_seats$sample_fraction, 
                              num_seats = 19)
  treatment_RF_exp_s19[remove_patients] <- predict_rf(RFExp_s19, expression_w12[remove_patients,])
  
  # RFExp seats=22
  RFExp_s22 <- train_rf_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,], 
                              num_trees = params_err_rfexp_seats$num_trees, 
                              mtry = params_err_rfexp_seats$mtry, 
                              min_node_size = params_err_rfexp_seats$min_node_size, 
                              max_depth = params_err_rfexp_seats$max_depth, 
                              sample_fraction = params_err_rfexp_seats$sample_fraction, 
                              num_seats = 22)
  treatment_RF_exp_s22[remove_patients] <- predict_rf(RFExp_s22, expression_w12[remove_patients,])
  
  # XGBoostExp
  # cat("XGBoostExp\n")
  # XGBoostExp <- train_xgb_hard(expression_w12[keep_patients,], drug_response_w12[keep_patients,],nrounds=nrounds)
  # treatment_XGBoost_exp[remove_patients] <- predict_xgb(XGBoostExp,expression_w12[remove_patients,])
  
  # XGBoostExp seats=1
  cat("XGBoostExp seats=1\n")
  XGBoostExp_s1 <- train_xgb_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,], 
                                   max_depth = params_err_xgboostexp_seats$max_depth, 
                                   eta = params_err_xgboostexp_seats$eta, 
                                   gamma = params_err_xgboostexp_seats$gamma, 
                                   min_child_weight = params_err_xgboostexp_seats$min_child_weight, 
                                   subsample = params_err_xgboostexp_seats$subsample, 
                                   colsample_bytree = params_err_xgboostexp_seats$colsample_bytree, 
                                   nrounds = params_err_xgboostexp_seats$nrounds, 
                                   early_stopping_rounds = 10, 
                                   num_seats = 1)
  treatment_XGBoost_exp_s1[remove_patients] <- predict_xgb(XGBoostExp_s1, expression_w12[remove_patients,])
  
  # XGBoostExp seats=6
  cat("XGBoostExp seats=6\n")
  XGBoostExp_s6 <- train_xgb_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,], 
                                   max_depth = params_err_xgboostexp_seats$max_depth, 
                                   eta = params_err_xgboostexp_seats$eta, 
                                   gamma = params_err_xgboostexp_seats$gamma, 
                                   min_child_weight = params_err_xgboostexp_seats$min_child_weight, 
                                   subsample = params_err_xgboostexp_seats$subsample, 
                                   colsample_bytree = params_err_xgboostexp_seats$colsample_bytree, 
                                   nrounds = params_err_xgboostexp_seats$nrounds, 
                                   early_stopping_rounds = 10, 
                                   num_seats = 6)
  treatment_XGBoost_exp_s6[remove_patients] <- predict_xgb(XGBoostExp_s6, expression_w12[remove_patients,])
  
  # XGBoostExp seats=16
  cat("XGBoostExp seats=16\n")
  XGBoostExp_s16 <- train_xgb_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,], 
                                    max_depth = params_err_xgboostexp_seats$max_depth, 
                                    eta = params_err_xgboostexp_seats$eta, 
                                    gamma = params_err_xgboostexp_seats$gamma, 
                                    min_child_weight = params_err_xgboostexp_seats$min_child_weight, 
                                    subsample = params_err_xgboostexp_seats$subsample, 
                                    colsample_bytree = params_err_xgboostexp_seats$colsample_bytree, 
                                    nrounds = params_err_xgboostexp_seats$nrounds, 
                                    early_stopping_rounds = 10, 
                                    num_seats = 16)
  treatment_XGBoost_exp_s16[remove_patients] <- predict_xgb(XGBoostExp_s16, expression_w12[remove_patients,])
  
  # XGBoostExp seats=19
  cat("XGBoostExp seats=19\n")
  XGBoostExp_s19 <- train_xgb_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,],
                                    max_depth = params_err_xgboostexp_seats$max_depth,
                                    eta = params_err_xgboostexp_seats$eta,
                                    gamma = params_err_xgboostexp_seats$gamma,
                                    min_child_weight = params_err_xgboostexp_seats$min_child_weight,
                                    subsample = params_err_xgboostexp_seats$subsample,
                                    colsample_bytree = params_err_xgboostexp_seats$colsample_bytree,
                                    nrounds = params_err_xgboostexp_seats$nrounds,
                                    early_stopping_rounds = 10,
                                    num_seats = 19)
  treatment_XGBoost_exp_s19[remove_patients] <- predict_xgb(XGBoostExp_s19, expression_w12[remove_patients,])

  # XGBoostExp seats=22
  cat("XGBoostExp seats=22\n")
  XGBoostExp_s22 <- train_xgb_seats(expression_w12[keep_patients,], pseudop_w12[keep_patients,],
                                    max_depth = params_err_xgboostexp_seats$max_depth,
                                    eta = params_err_xgboostexp_seats$eta,
                                    gamma = params_err_xgboostexp_seats$gamma,
                                    min_child_weight = params_err_xgboostexp_seats$min_child_weight,
                                    subsample = params_err_xgboostexp_seats$subsample,
                                    colsample_bytree = params_err_xgboostexp_seats$colsample_bytree,
                                    nrounds = params_err_xgboostexp_seats$nrounds,
                                    early_stopping_rounds = 10,
                                    num_seats = 22)
  treatment_XGBoost_exp_s22[remove_patients] <- predict_xgb(XGBoostExp_s22, expression_w12[remove_patients,])
}

# Change to drug names 
treatment_ODT_mut <- colnames(drug_response_w12)[as.numeric(treatment_ODT_mut)]
treatment_ODT_exp <- colnames(drug_response_w12)[as.numeric(treatment_ODT_exp)]

# Obtain oracle
treatmentOracle <- colnames(drug_response_w12)[apply(drug_response_w12, 1, which.min)]

## 5.3 Visualize Results ---
### Plotting BeatAML Cross-validation Results and comparing means difference
treatment_plot <- data.frame(Method = "ORACLE", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatmentOracle)])
treatment_plot <- rbind(treatment_plot, data.frame(Method = "ODTMut", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_ODT_mut)]))
# treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_mut)]))
# treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_mut)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "ODTExp", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_ODT_exp)]))
# treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_exp)]))
# treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_exp)]))

# Adding plots for RFMut with num_seats
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut_s1", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_mut_s1)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut_s6", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_mut_s6)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut_s16", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_mut_s16)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut_s19", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_mut_s19)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFMut_s22", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_mut_s22)]))

# Adding plots for XGBoostMut with num_seats
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut_s1", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_mut_s1)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut_s6", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_mut_s6)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut_s16", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_mut_s16)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut_s19", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_mut_s19)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostMut_s22", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_mut_s22)]))

# Adding plots for RFExp with num_seats
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp_s1", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_exp_s1)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp_s6", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_exp_s6)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp_s16", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_exp_s16)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp_s19", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_exp_s19)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "RFExp_s22", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_RF_exp_s22)]))

# Adding plots for XGBoostExp with num_seats
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp_s1", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_exp_s1)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp_s6", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_exp_s6)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp_s16", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_exp_s16)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp_s19", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_exp_s19)]))
treatment_plot <- rbind(treatment_plot, data.frame(Method = "XGBoostExp_s22", IC50 = drug_response_w12[cbind(rownames(drug_response_w12), treatment_XGBoost_exp_s22)]))

treatment_plot$Method <- factor(treatment_plot$Method, levels = c(
  "ORACLE", "ODTMut", "RFMut", "RFMut_s1", "RFMut_s6", "RFMut_s16", "RFMut_s19", "RFMut_s22",
  "XGBoostMut", "XGBoostMut_s1", "XGBoostMut_s6", "XGBoostMut_s16", "XGBoostMut_s19", "XGBoostMut_s22",
  "ODTExp", "RFExp", "RFExp_s1", "RFExp_s6", "RFExp_s16", "RFExp_s19", "RFExp_s22",
  "XGBoostExp", "XGBoostExp_s1", "XGBoostExp_s6", "XGBoostExp_s16", "XGBoostExp_s19", "XGBoostExp_s22"
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
ggsave(paste(folder_dir, "/images/boxplot_validate_final.png", sep = ""), treatment_plot_gg, width = 12, height = 7, dpi = 500)

#ggsave("/Users/katyna/Desktop/boxplot_validate_final.png", treatment_plot_gg, width = 12, height = 7, dpi = 500)


# Calculate mean and median IC50 for each method
summary_stats_val <- treatment_plot %>%
  group_by(Method) %>%
  summarise(
    Mean_IC50 = mean(IC50, na.rm = TRUE),
    Median_IC50 = median(IC50, na.rm = TRUE)
  ) %>%
  arrange(Mean_IC50)
write_xlsx(summary_stats_val, paste(folder_dir, "/ic50_results_validation_final.xlsx", sep = ""))

#===============================================================================
# 6. Model Training on Complete Dataset

# Training models in complete BeatAML Waves1+2 data, testing in Waves3+4 and timing
"
We trained the models in BeatAML cohort and measured the training time for each model.
Technical implementation refers to the computational burden and software that the method requires.Despite it could be considered less important, some of the algorithms require hours of computing time for the BeatAML of the subset of AML samples in GDSC -that be considered to be small/medium size.
" 
## 6.1 Training Models

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
# cat("RFMut\n")
# start_time <- Sys.time()
# #RFMut_All <- train_rf_hard(mutations_w12, drug_response_w12,num_trees=num_trees)
# RFMut_All <- train_rf_hard(mutations_w12, drug_response_w12, num_trees = 100, mtry = 4, min_node_size = 10, max_depth = 100, sample_fraction = 1)
# treatmentRFMut_w34 <- predict_rf(RFMut_All, mutations_w34)
# end_time <- Sys.time()
# time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
# num_variables <- length(RFMut_All$forest$independent.variable.names)
# timing_results <- rbind(timing_results, data.frame(Method = "RFMut", Time = time_taken, Variables = num_variables))

# XGBoostMut
# cat("XGBoostMut\n")
# start_time <- Sys.time()
# XGBoostMut_All <- train_xgb_hard(mutations_w12, drug_response_w12, max_depth = 3, eta = 0.3, gamma = 0, min_child_weight = 1, subsample = 1, colsample_bytree = 1, nrounds = 100, nthread = 8, early_stopping_rounds = 20)
# treatmentXGBoostMut_w34 <- predict_xgb(XGBoostMut_All, mutations_w34)
# end_time <- Sys.time()
# time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
# num_variables <- length(XGBoostMut_All$feature_names)
# timing_results <- rbind(timing_results, data.frame(Method = "XGBoostMut", Time = time_taken, Variables = num_variables))

# RFMut with num_seats
for (num_seats in c(1, 6, 16, 19, 22)) {
  cat(paste("RFMut seats=", num_seats, "\n", sep = ""))
  start_time <- Sys.time()
  RFMut_All <- train_rf_seats(mutations_w12,  pseudop_w12, num_trees = params_err_rfmut_seats$num_trees, mtry = params_err_rfmut_seats$mtry, min_node_size = params_err_rfmut_seats$min_node_size, max_depth = params_err_rfmut_seats$max_depth, sample_fraction = params_err_rfmut_seats$sample_fraction, num_seats = num_seats)
  assign(paste0("RFMut_All_s", num_seats), RFMut_All)  # Assign model to a unique variable name
  assign(paste0("treatmentRFMut_w34_s", num_seats), predict_rf(RFMut_All, mutations_w34))
  end_time <- Sys.time()
  time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
  num_variables <- length(RFMut_All$forest$independent.variable.names)
  timing_results <- rbind(timing_results, data.frame(Method = paste0("RFMut_s", num_seats), Time = time_taken, Variables = num_variables))
}

# XGBoostMut with num_seats
for (num_seats in c(1, 6, 16, 19, 22)) {
  cat(paste("XGBoostMut vote=", num_seats, "\n", sep = ""))
  start_time <- Sys.time()
  XGBoostMut_All <- train_xgb_seats(mutations_w12, pseudop_w12, max_depth = params_err_xgboostmut$max_depth, eta = params_err_xgboostmut$eta, gamma = params_err_xgboostmut$gamma, min_child_weight = params_err_xgboostmut$min_child_weight, subsample = params_err_xgboostmut$subsample, colsample_bytree = params_err_xgboostmut$colsample_bytree, nrounds = params_err_xgboostmut$nrounds, early_stopping_rounds = 10, num_seats = num_seats)
  assign(paste0("XGBoostMut_All_s", num_seats), XGBoostMut_All)  # Assign model to a unique variable name
  assign(paste0("treatmentXGBoostMut_w34_s", num_seats), predict_xgb(XGBoostMut_All, mutations_w34))
  end_time <- Sys.time()
  time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
  num_variables <- length(XGBoostMut_All$feature_names)
  timing_results <- rbind(timing_results, data.frame(Method = paste0("XGBoostMut_s", num_seats), Time = time_taken, Variables = num_variables))
}

# ODTExp
cat("ODTExp\n")
start_time <- Sys.time()
ODTExp_All <- trainTree(PatientData = expression_w12, PatientSensitivity = drug_response_w12, minbucket = minbucket_exp)
treatmentODTExp_w34 <- predictTree(tree = ODTExp_All, PatientSensitivityTrain = drug_response_w12, PatientData = expression_w34)
end_time <- Sys.time()
time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
timing_results <- rbind(timing_results, data.frame(Method = "ODTExp", Time = time_taken, Variables = 0))

# # RFExp
# cat("RFExp\n")
# start_time <- Sys.time()
# RFExp_All <- train_rf_hard(expression_w12, drug_response_w12, num_trees = 100, mtry = 40, min_node_size = 5, max_depth = 100, sample_fraction = 0.8)

# treatmentRFExp_w34 <- predict_rf(RFExp_All, expression_w34)
# end_time <- Sys.time()
# time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
# num_variables <- length(RFExp_All$forest$independent.variable.names)
# timing_results <- rbind(timing_results, data.frame(Method = "RFExp", Time = time_taken, Variables = num_variables))

# RFExp with num_seats
for (num_seats in c(1, 6, 16, 19, 22)) {
  cat(paste("RFExp seats=", num_seats, "\n", sep = ""))
  start_time <- Sys.time()
  RFExp_All <- train_rf_seats(expression_w12, pseudop_w12, 
                              num_trees = params_err_rfexp_seats$num_trees, 
                              mtry = params_err_rfexp_seats$mtry, 
                              min_node_size = params_err_rfexp_seats$min_node_size, 
                              max_depth = params_err_rfexp_seats$max_depth, 
                              sample_fraction = params_err_rfexp_seats$sample_fraction, 
                              num_seats = num_seats)
  assign(paste0("RFExp_All_s", num_seats), RFExp_All)  # Assign model to a unique variable name
  assign(paste0("treatmentRFExp_w34_s", num_seats), predict_rf(RFExp_All, expression_w34))
  end_time <- Sys.time()
  time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
  num_variables <- length(RFExp_All$forest$independent.variable.names)
  timing_results <- rbind(timing_results, data.frame(Method = paste0("RFExp_s", num_seats), Time = time_taken, Variables = num_variables))
}

# # XGBoostExp
# cat("XGBoostExp\n")
# start_time <- Sys.time()
# XGBoostExp_All <- train_xgb_hard(expression_w12, drug_response_w12, max_depth = 3, eta = 0.3, gamma = 0, min_child_weight = 1, subsample = 1, colsample_bytree = 1, nrounds = 100, nthread = 8, early_stopping_rounds = 20)

# treatmentXGBoostExp_w34 <- predict_xgb(XGBoostExp_All, expression_w34)
# end_time <- Sys.time()
# time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
# num_variables <- length(XGBoostExp_All$feature_names)
# timing_results <- rbind(timing_results, data.frame(Method = "XGBoostExp", Time = time_taken, Variables = num_variables))

# XGBoostExp with num_seats
for (num_seats in c(1, 6, 16, 19, 22)) {
  cat(paste("XGBoostExp seats=", num_seats, "\n", sep = ""))
  start_time <- Sys.time()
  XGBoostExp_All <- train_xgb_seats(expression_w12, pseudop_w12, 
                  max_depth = params_err_xgboostexp_seats$max_depth, 
                  eta = params_err_xgboostexp_seats$eta, 
                  gamma = params_err_xgboostexp_seats$gamma, 
                  min_child_weight = params_err_xgboostexp_seats$min_child_weight, 
                  subsample = params_err_xgboostexp_seats$subsample, 
                  colsample_bytree = params_err_xgboostexp_seats$colsample_bytree, 
                  nrounds = params_err_xgboostexp_seats$nrounds, 
                  early_stopping_rounds = 10, 
                  num_seats = num_seats)
  assign(paste0("XGBoostExp_All_s", num_seats), XGBoostExp_All)  # Assign model to a unique variable name
  assign(paste0("treatmentXGBoostExp_w34_s", num_seats), predict_xgb(XGBoostExp_All, expression_w34))
  end_time <- Sys.time()
  time_taken <- sprintf("%.3f sec", as.numeric(difftime(end_time, start_time, units = "secs")))
  num_variables <- length(XGBoostExp_All$feature_names)
  timing_results <- rbind(timing_results, data.frame(Method = paste0("XGBoostExp_s", num_seats), Time = time_taken, Variables = num_variables))
}

## 6.2 Visualize Results 

# Save timing results to an Excel file
write_xlsx(timing_results, paste(folder_dir, "/timing_results_final.xlsx", sep = ""))

# Change to drug names 
treatmentODTMut_w34 <- colnames(drug_response_w12)[treatmentODTMut_w34]
treatmentODTExp_w34 <-colnames(drug_response_w12)[treatmentODTExp_w34]

# Obtain oracle
treatmentOracle_w34 <- colnames(drug_response_w34)[apply(drug_response_w34, 1, which.min)]

### Plotting BeatAML Cross-validation Results and comparing means difference
treatment_plot_test <- data.frame(Method = "ORACLE", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentOracle_w34)])
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "ODTMut", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentODTMut_w34)]))
#treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFMut", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFMut_w34)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFMut_s1", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFMut_w34_s1)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFMut_s6", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFMut_w34_s6)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFMut_s16", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFMut_w34_s16)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFMut_s19", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFMut_w34_s19)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFMut_s22", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFMut_w34_s22)]))

#treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostMut", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostMut_w34)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostMut_s1", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostMut_w34_s1)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostMut_s6", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostMut_w34_s6)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostMut_s16", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostMut_w34_s16)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostMut_s19", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostMut_w34_s19)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostMut_s22", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostMut_w34_s22)]))

treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "ODTExp", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentODTExp_w34)]))
#treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFExp", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFExp_w34)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFExp_s1", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFExp_w34_s1)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFExp_s6", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFExp_w34_s6)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFExp_s16", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFExp_w34_s16)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFExp_s19", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFExp_w34_s19)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "RFExp_s22", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentRFExp_w34_s22)]))

#treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostExp", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostExp_w34)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostExp_s1", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostExp_w34_s1)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostExp_s6", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostExp_w34_s6)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostExp_s16", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostExp_w34_s16)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostExp_s19", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostExp_w34_s19)]))
treatment_plot_test <- rbind(treatment_plot_test, data.frame(Method = "XGBoostExp_s22", IC50 = drug_response_w34[cbind(rownames(drug_response_w34), treatmentXGBoostExp_w34_s22)]))

treatment_plot_test$Method <- factor(treatment_plot_test$Method, levels = c(
  "ORACLE", "ODTMut", "RFMut_s1", "RFMut_s6", "RFMut_s16", "RFMut_s19", "RFMut_s22",
  "XGBoostMut_s1", "XGBoostMut_s6", "XGBoostMut_s16", "XGBoostMut_s19", "XGBoostMut_s22",
  "ODTExp", "RFExp_s1", "RFExp_s6", "RFExp_s16", "RFExp_s19", "RFExp_s22",
  "XGBoostExp_s1", "XGBoostExp_s6", "XGBoostExp_s16", "XGBoostExp_s19", "XGBoostExp_s22"
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


ggsave(paste(folder_dir, "/images/boxplot_test_final.png", sep = ""), treatment_plot_test_gg, width = 12, height = 7, dpi = 500)
ggsave("/Users/katyna/Desktop/w34.png", treatment_plot_test_gg, width = 12, height = 7, dpi = 500)


# Calculate mean and median IC50 for each method
summary_stats_test <- treatment_plot_test %>%
  group_by(Method) %>%
  summarise(
    Mean_IC50 = mean(IC50, na.rm = TRUE),
    Median_IC50 = median(IC50, na.rm = TRUE)
  ) %>%
  arrange(Mean_IC50)
write_xlsx(summary_stats_test, paste(folder_dir, "/ic50_results_test_final.xlsx", sep = ""))

# Save results for GDSC
save(ODTMut_All, RFMut_All_s1, RFMut_All_s6, RFMut_All_s16, RFMut_All_s19, RFMut_All_s22,
     XGBoostMut_All_s1, XGBoostMut_All_s6, XGBoostMut_All_s16, XGBoostMut_All_s19, XGBoostMut_All_s22,
     pseudop_w12, drug_response_w12,
     file = paste(folder_dir,"/Rdata/models_all_final.RData",sep=""))

# # Save results
# save(ODTMut_All, RFMut_All, 
#      XGBoostMut_All, 
#      ODTExp_All, RFExp_All, 
#      XGBoostExp_All, drug_response_w12,
#      file = paste(folder_dir,"/Rdata/models_all_5000_0_7.RData",sep=""))

# save(treatmentODTMut_w34, treatmentRFMut_w34, treatmentRFMut_w34_s3, treatmentRFMut_w34_s5, treatmentRFMut_w34_s10, treatmentRFMut_w34_s20, 
#   treatmentXGBoostMut_w34, treatmentXGBoostMut_w34_s3, treatmentXGBoostMut_w34_s5, treatmentXGBoostMut_w34_s10, treatmentXGBoostMut_w34_s20, 
#   treatmentODTExp_w34, treatmentRFExp_w34, treatmentRFExp_w34_s3, treatmentRFExp_w34_s5, treatmentRFExp_w34_s10, treatmentRFExp_w34_s20, 
#   treatmentXGBoostExp_w34, treatmentXGBoostExp_w34_s3, treatmentXGBoostExp_w34_s5, treatmentXGBoostExp_w34_s10, treatmentXGBoostExp_w34_s20, 
#   file = paste(folder_dir,"/Rdata/results_test_5000_0_7.RData",sep=""))


# Load results
# load(file = paste(folder_dir,"Rdata/mut_models_all.RData",sep=""))
# load(file = paste(folder_dir,"Rdata/mut_results_test.RData",sep=""))  

#

# Compute statistical differences ---------------------------
# Load required libraries
library(gtools)

# Assuming validation_plot is your data frame
# Calculate means by method
means <- treatment_plot_test %>%
  group_by(Method) %>%
  summarise(Mean_IC50 = mean(IC50, na.rm = TRUE)) %>%
  ungroup()

# Function to perform Wilcoxon test and return p-value
wilcoxon_test_pair <- function(x, y) {
  test <- try(wilcox.test(x, y, exact = FALSE), silent = TRUE)
  if (inherits(test, "try-error")) {
    return(NA)  # Return NA if the test fails
  } else {
    return(test$p.value)
  }
}

# List of all unique methods
methods <- unique(treatment_plot_test$Method)

# Generate all possible pairs of methods
method_pairs <- combn(methods, 2, simplify = FALSE)

# Apply the Wilcoxon test to each pair and store results in a data frame
pairwise_tests <- data.frame(Method1 = character(), Method2 = character(), p_value = numeric())
for (i in seq_along(method_pairs)) {
  pair <- method_pairs[[i]]
  pairwise_tests <- rbind(pairwise_tests, data.frame(
    Method1 = pair[1],
    Method2 = pair[2],
    p_value = wilcoxon_test_pair(treatment_plot_test$IC50[treatment_plot_test$Method == pair[1]],
                                 treatment_plot_test$IC50[treatment_plot_test$Method == pair[2]])
  ))
}

# Print the results data frame
View(pairwise_tests)


#===============================================================================
# 7. Model Explainability Plots and Perfomance Metrics

"
Feature importance and graphic
"
# ODT
niceTree(ODTMut_All)
niceTree(ODTExp_All)

# RF - RFExp_s22
RFExp_All <- train_rf_seats(expression_w12, pseudop_w12, 
                            num_trees = params_err_rfexp_seats$num_trees, 
                            mtry = params_err_rfexp_seats$mtry, 
                            min_node_size = params_err_rfexp_seats$min_node_size, 
                            max_depth = params_err_rfexp_seats$max_depth, 
                            sample_fraction = params_err_rfexp_seats$sample_fraction, 
                            num_seats = 6)

plot(RFExp_All)
varImpPlot(RFExp_All)

# XGBoost - XGBoostExp_s6
# Extract variable importance
XGBoostExp_All <- train_xgb_seats(expression_w12, pseudop_w12, 
                                  max_depth = params_err_xgboostexp_seats$max_depth, 
                                  eta = params_err_xgboostexp_seats$eta, 
                                  gamma = params_err_xgboostexp_seats$gamma, 
                                  min_child_weight = params_err_xgboostexp_seats$min_child_weight, 
                                  subsample = params_err_xgboostexp_seats$subsample, 
                                  colsample_bytree = params_err_xgboostexp_seats$colsample_bytree, 
                                  nrounds = params_err_xgboostexp_seats$nrounds, 
                                  early_stopping_rounds = 10, 
                                  num_seats = 6)

importance_matrix <- xgb.importance(feature_names = colnames(expression_w12), model = XGBoostExp_All)
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
  Model = "RFExp_s10",
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
  mutate(Model = "XGBoostExp_s5") %>%
  select(Genes, Importance, Model)  # Select only the necessary columns

# Combine the data frames
combined_data <- rbind(var_importance_rf_df, importance_df_xgb)

# Define custom colors for each model
model_colors <- c("RFExp_s10" = "#96D6B6", "XGBoostExp_s5" = "#FDE5B0")

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

## 7.2 Performance Visualization --- 

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


# MATRICES OF DRUG ASSIGNMENT --------------------------------------------
# 
# fakeBeat <- matrix(rnorm(250 * 119), 250, 119) # patient x drugs (measured)
# fakeBeat2 <- 1*fakeBeat + 1*matrix(rnorm(250 * 119), 250, 119) # predicted
# indication <- apply(fakeBeat2,1, which.max) # suggested indication
# rankfakeBeat <- t(apply(-fakeBeat, 1, function(x) rank(x, ties.method = "first")))
# table(factor(rankfakeBeat[cbind(1:length(indication),indication)],1:119) )

# rank the real/measured drug responses to see which is the top drug for a patient, then the second best drug, until the worst 
rank_drug_response_w34 <- t(apply(drug_response_w34, 1, function(x) rank(x, ties.method = "first")))

plot(table(factor(rank_drug_response_w34[cbind(rownames(drug_response_w34),treatmentOracle_w34)],1:119) ))

# Assuming drug_response_w34 is your matrix of drug responses (patients x drugs)
rank_drug_response_w34 <- t(apply(drug_response_w34, 1, function(x) rank(x, ties.method = "first")))

# List of all treatment vectors
treatment_vectors <- list(
  Oracle = treatmentOracle_w34,
  ODTMut = treatmentODTMut_w34,
  RFMut_s1 = treatmentRFMut_w34_s1,
  RFMut_s6 = treatmentRFMut_w34_s6,
  RFMut_s16 = treatmentRFMut_w34_s16,
  RFMut_s19 = treatmentRFMut_w34_s19,
  RFMut_s22 = treatmentRFMut_w34_s22,
  XGBoostMut_s1 = treatmentXGBoostMut_w34_s1,
  XGBoostMut_s6 = treatmentXGBoostMut_w34_s6,
  XGBoostMut_s16 = treatmentXGBoostMut_w34_s16,
  XGBoostMut_s19 = treatmentXGBoostMut_w34_s19,
  XGBoostMut_s22 = treatmentXGBoostMut_w34_s22,
  ODTExp = treatmentODTExp_w34,
  RFExp_s1 = treatmentRFExp_w34_s1,
  RFExp_s6 = treatmentRFExp_w34_s6,
  RFExp_s16 = treatmentRFExp_w34_s16,
  RFExp_s19 = treatmentRFExp_w34_s19,
  RFExp_s22 = treatmentRFExp_w34_s22,
  XGBoostExp_s1 = treatmentXGBoostExp_w34_s1,
  XGBoostExp_s6 = treatmentXGBoostExp_w34_s6,
  XGBoostExp_s16 = treatmentXGBoostExp_w34_s16,
  XGBoostExp_s19 = treatmentXGBoostExp_w34_s19,
  XGBoostExp_s22 = treatmentXGBoostExp_w34_s22
)

# Initialize an empty data frame with ranks as columns
rank_counts <- data.frame(matrix(0, nrow = length(treatment_vectors), ncol = 119))
colnames(rank_counts) <- 1:119
rownames(rank_counts) <- names(treatment_vectors)

oracle_labels <-  treatment_vectors[["Oracle"]]

# Iterate over each treatment method
for (method in names(treatment_vectors)) {
  indication <- treatment_vectors[[method]]
  ranks <- rank_drug_response_w34[cbind(rownames(drug_response_w34), indication)]
  
  # Count the ranks
  rank_table <- table(factor(ranks, levels = 1:119))
  
  # Update the rank_counts data frame
  rank_counts[method, ] <- rank_table
}

# Print the results
print(rank_counts)

# Heatmap of rank counts
# Load necessary library
# Load necessary libraries
library(pheatmap)
library(grid)
library(gridExtra)

# Assuming rank_counts_matrix is already defined and is a matrix
# Exclude the first row and select the first 20 columns for plotting
data_to_plot <- rank_counts_matrix[-1, 1:20]

# Create the heatmap with titles and axis labels
p <- pheatmap(
  mat = data_to_plot,
  cluster_cols = FALSE,       # Disable column clustering
  cluster_rows = FALSE,       # Disable row clustering
  main = "Drug Ranking Heatmap (Waves 3+4)", # Title of the heatmap
  fontsize = 24,              # Font size for text
  fontsize_row = 24,           # Font size for row labels
  fontsize_col = 24,           # Font size for column labels
  silent = TRUE               # Suppress drawing
)

# Extract the gtable from the pheatmap object
heatmap_gtable <- p$gtable

# Create custom x and y axis labels
x_label_grob <- textGrob("Drug Ranking", gp=gpar(fontsize=24), hjust=1 ,vjust=0)
y_label_grob <- textGrob("Model", gp=gpar(fontsize=24), rot=90, vjust=1)

# Combine the heatmap with labels using grid.arrange
grid.newpage()
combined_gtable <- arrangeGrob(
  arrangeGrob(y_label_grob, heatmap_gtable, ncol=2, widths=c(1, 8)),
  x_label_grob,
  nrow=2,
  heights=c(8,1)
)
grid.draw(combined_gtable)

# Save the heatmap with a transparent background
ggsave(paste(folder_dir, "/images/drug_ranking_heatmap.png", sep = ""), plot = combined_gtable, 
       bg = "transparent", width = 10, height = 8, dpi = 300)


