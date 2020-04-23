# Requirements 
Code was tested on R version 3.4.3

R libraries:
* MatrixEQTL
* lattice
* biomaRt
* caret
* randomForest
* dplyr
* Boruta

# Data colllections
Large data files are gzipped and should be unzipped before usage. on github only two datasets out of 4 analysed in the paper are available (genomics data partly).
## Covariates
/data/covariates/
* Covariates_UHPOS.txt
* Covariates_NMR_Urine.txt
* Covariates_EigenMS_format.txt
* Covariates_diagnosis.txt

## Genomics data
/data/dna_matrices/

On github only first parts of dna matrices are available (500,000 SNP records) due to size limitations, full versions can be obtained from the synapse or from the corresponding author by request.

* dna_matrix_UHPOS.tsv.gz
* dna_matrix_NMR_Urine.tsv.gz

## Original metabolomics data
/data/metabolomics_matrices/
* NMR_Urine_original.txt.gz
* UHPOS_original.csv

## Normalised metabolomics data
/data/metabolomics_matrices/normalised/
* NMR_Urine.tsv.gz
* UHPOS.tsv

## Permuted by samples metabolimics data
/data/metabolomics_matrices/normalised/permuted/
* NMR_Urine_permuted.tsv.gz
* UHPOS_permuted.tsv

The full dataset availability:
Data from the AddNeuroMed study, including clinical and molecular data is available at Sage Bionetworks AD community portal https://www.synapse.org/#!Synapse:syn2790911/wiki/235388. Data from the current MS will be added to this dataset on publication.

# Metabolomics data normalisation, MS data
```
source("./scripts/normalization.R")
```
## 1) EigenMS to remove bias of unknown complexity from MS data (“Metabolomics Data Normalization with EigenMS”, Yuliya V. Karpievitch et al.) 
```
EigenMS_normalize(matrix_to_normalize, matrix_EigenMS_result, covariates_matrix_to_preserve)

EigenMS_normalize("./data/metabolomics_matrices/UHPOS_original.csv","./data/metabolomics_matrices/normalised/UHPOS_EigenMS.csv", "./data/covariates/Covariates_EigenMS_format.txt")
```
## 2) Quantile normalization to transform measurements for each metabolite into normally distributed while preserving relative rankings.
```
QN_normalize(matrix_to_normalize, matrix_QN_result)

QN_normalize("./data/metabolomics_matrices/normalised/UHPOS_EigenMS.csv","./data/metabolomics_matrices/normalised/UHPOS.tsv")
```
# Metabolomics data normalisation, NMR data
```
source("./scripts/normalization.R")
```
## 1) Quantile normalization to transform measurements for each metabolite into normally distributed while preserving relative rankings.
```
QN_normalize(matrix_to_normalize, matrix_QN_result)

QN_normalize("./data/metabolomics_matrices/NMR_Nosey_Urine_original.txt","./data/metabolomics_matrices/normalised/NMR_Nosey_Urine.tsv")
```

# Metabolic QTL analysis

We've run this script on a computational cluster using parallel computing and dna matrix separation into chunks. 

```
./scripts/run_mqtl.R "NMR_Urine"
```

Results: ./results/mqtl_NMR_Urine.txt.gz

Feature reduction using mQTL significance (FDR threshold 0.01) results with genomic annotations (gene-metabolite): ./results/result_NMR_Urine_annot.txt

Feature reduction using mQTL significance (FDR threshold 0.01) results with genomic annotations and SNPs (SNP-metabolite): ./results/result_NMR_Urine_annot_full.txt


Select metabolites based on mQTL significance:

```
source("./scripts/QTL_results_selection.R")

feature_selection_mqtl_per_assay("UHPOS")
feature_selection_mqtl_per_assay("NMR_Urine")

```

When all data are available:
```
source("./scripts/QTL_results_selection.R")
feature_selection_mqtl_per_assay("UHPOS")
feature_selection_mqtl_per_assay("URPOS")
feature_selection_mqtl_per_assay("URNEG")
feature_selection_mqtl_per_assay("NMR_Urine")

feature_selection_mqtl_snps()
feature_selection_mqtl_metabolites()
```
Intersection between GWAS SNPs from Alzheimer's pathology related traits and mQTL SNPs:

```
source("./scripts/QTL_results_selection.R")
snps_intersection("./results/gwas_snps.txt")
```
# Feature selection with Random Forests

We've tried multiple feature selection methods: linear regression, feature selection by correlation (WEKA) and Random Forests (RF) for feature selection. RF gave us the best results.

NB! RF is used to select the most important metabolomic features for further annotation, not to create the prediction model. We are aware of overfitting and do not proposing to use the final model to predict Alzheimer's disease. One more time, RF is used for the selection of the important features. See https://chrisalbon.com/machine_learning/trees_and_forests/feature_selection_using_random_forest/, https://doi.org/10.1016/j.csda.2012.09.020, etc.

```
source("./scripts/feature_selection.R")
./scripts/RF_choose_set_and_classifier.R
./scripts/RF_parameters_tuning.R
```

# Additional scripts (GWAS Catalog search, plots, etc.)

./scripts/additional/

