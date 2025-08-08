# Trees & Democracy: SEATS and ODT for Personalized Drug Assignment in AML

## Abstract
Personalized drug selection is crucial for treating complex diseases like Acute Myeloid Leukemia, where maximizing therapeutic efficacy is essential. Although precision medicine aims to tailor treatments to individual molecular profiles, existing machine learning models often fall short in selecting the best drug from multiple candidates. We present SEATS (Systematic Efficacy Assignment with Treatment Seats), which adapts conventional models like Random Forest and XGBoost for multiclass drug assignment by allocating probabilistic “treatment seats” to drugs based on efficacy. This approach helps models learn clinically relevant distinctions.

Additionally, we assess an interpretable Optimal Decision Tree (ODT) model designed specifically for drug assignment. Trained on the BeatAML2 cohort and validated on the GDSC AML cell line dataset, integrating SEATS with Random Forest and XGBoost improved prediction accuracy and consistency. The ODT model offered competitive performance with clear, interpretable decision paths and minimal feature requirements, facilitating clinical use. SEATS reorients standard models towards personalized drug selection, and combined with the ODT framework, provides effective, interpretable strategies for precision oncology, underscoring the potential of tailored machine learning solutions in supporting real-world treatment decisions.

## Repository overview
This repository contains R code to reproduce the analysis workflow for learning patient-specific drug assignments using:
- Random Forest and XGBoost with SEATs
- Optimal Decision Trees

Main scripts:
- code/TreesAndDemocracy.R — end-to-end analysis with accuracy-oriented hyperparameters (cross-validation on BeatAML Waves 1+2, testing on Waves 3+4, plots, timing, importance)
- code/TreesAndDemocracy_final_params.R — same workflow with a different parameter set reference
- code/TreesAndDemocracy_acc.R — accuracy-oriented configuration (referenced in the main script)
- code/parameters.R — hyperparameters for RF/XGBoost SEATS and others (sourced by the main scripts)
- code/XAIfunctions.R — helper functions for model training, prediction, plotting (sourced by the main scripts)

Outputs are written to TreesAndDemocracy/output (images, xlsx summaries, RData models).

## Data expected by the code
Paths are referenced relative to a OneDrive folder and a sibling project folder. The scripts expect the following Excel files (BeatAML2) to exist under:
GithubProjects/interpretable-ai-review/data/input/BeatAML2/
- beataml_wv1to4_clinical.xlsx (sheet: summary)
- beataml_waves1to4_norm_exp_dbgap.xlsx
- beataml_wes_wv1to4_mutations_dbgap.xlsx
- beataml_probit_curve_fits_v4_dbgap.xlsx

Note: The scripts remove toxic drugs Elesclomol and JNJ-7706621 from the drug matrix.

## Methods implemented (as per code)
- Data processing:
  - Build patient-by-drug IC50 matrix, log-transform and center/scale by drug, derive IC50*; kNN imputation for missing values
  - Construct binary mutation matrix (genes + clinical translocations)
  - Filter expression to top 5,000 most variable genes; keep mutations present in ≥1% of patients
  - Split BeatAML into Waves 1+2 (train/validation) and Waves 3+4 (test)
- SEATS pseudoprobabilities: pseudop = softmax(-gamma · IC50*); gamma tried at 1.6, 10, 20 (pseudop is used by RF/XGB SEATS trainers)
- Models:
  - ODT: trainTree/predictTree on mutations and on expression
  - RF SEATS: train_rf_seats/predict_rf with num_seats ∈ {1, 6, 16, 19, 22}
  - XGBoost SEATS: train_xgb_seats/predict_xgb with num_seats ∈ {1, 6, 16, 19, 22}
- Evaluation:
  - 5-fold CV on Waves 1+2; oracle = argmin(IC50*) per patient
  - Hold-out test on Waves 3+4
  - Plots saved to output/images: CV boxplots, test boxplots, feature importance, implementability (time/variables), and a drug ranking heatmap
  - Summary statistics written to Excel (ic50_results_validation_final.xlsx, ic50_results_test_final.xlsx, timing_results_final.xlsx)
  - Trained models saved to output/Rdata/models_all_final.RData

## Software requirements
R (≥ 4.0 recommended) and the following packages used by the scripts:
- readxl, writexl, matrixStats, partykit, impute, ggplot2, ggpattern, pheatmap, dplyr, ODT, gtools, tidyr
- For some figures: system fonts support for “Roboto” (optional in code)
- XGBoost and Random Forest backends are invoked through helper functions defined in code/XAIfunctions.R

Install packages in R:
```r
install.packages(c("readxl","writexl","matrixStats","partykit","impute","ggplot2","ggpattern","pheatmap","dplyr","gtools","tidyr"))
# ODT and xgboost may require:
install.packages("ODT")
install.packages("xgboost")
```

## Usage
1) Verify working directory and paths
   - The scripts set the working directory to a OneDrive root via variables mac/angel and source helper files using relative paths like ./GithubProjects/Trees&Democracy/...
   - If your local paths differ, update the variables at the top of code/TreesAndDemocracy.R (and the other scripts):
     mac <- "/path/to/OneDrive-root/"
     computer <- mac
     setwd(computer)
   - Ensure the BeatAML2 Excel files are present under GithubProjects/interpretable-ai-review/data/input/BeatAML2/ as listed above.

2) Run the main analysis
   - From R:
     source("./GithubProjects/Trees&Democracy/code/TreesAndDemocracy.R")
   - This will:
     - Load data and build matrices
     - Perform 5-fold CV on Waves 1+2 for ODT, RF/XGB with SEATS
     - Train on Waves 1+2 and test on Waves 3+4
     - Save figures and Excel summaries under TreesAndDemocracy/output

3) Alternative parameterization
   - You can run code/TreesAndDemocracy_final_params.R to execute the same workflow while sourcing parameters from code/parameters.R.

4) Outputs
   - Images: TreesAndDemocracy/output/images/*.png
   - Tables: TreesAndDemocracy/output/*.xlsx
   - Models: TreesAndDemocracy/output/Rdata/models_all_final.RData

## Notes and limitations
- Some “hard” RF/XGBoost baselines are present but commented out in the scripts; the active comparisons primarily use SEATS variants.
- The helper functions (train_rf_seats, train_xgb_seats, predict_rf, predict_xgb, niceTree, etc.) are sourced from code/XAIfunctions.R and must be available.
- The repository assumes access to BeatAML2 data files; these are not included here.
