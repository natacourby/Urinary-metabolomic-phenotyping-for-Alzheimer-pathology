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
* gwascat (Bioconductor)

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

We've tried multiple feature selection methods: linear regression, Recursive Feature Elimination (RFE), Learning Vector Quantization (LVQ), feature selection by correlation and Random Forests (RF). RF gave us the best results.

NB! RF is used to select the most important metabolomic features for further annotation, not to create the prediction model. We are aware of overfitting and do not proposing to use the final model to predict Alzheimer's disease. One more time, RF is used for the selection of the important features. See https://chrisalbon.com/machine_learning/trees_and_forests/feature_selection_using_random_forest/, https://doi.org/10.1016/j.csda.2012.09.020, etc.

RF method was chosen for solution of the classification problem since it performs implicit feature selection and provides a good indicator of feature importance for the ranking. RF has relatively good accuracy and robustness and it is easy to apply.

The main focus of this study is on metabolites. However, there are also 6923 SNPs that metabolites are associated with and 9 SNPs found in the results of the GWAS. Finally, covariates available for the samples: age, gender and data collection centre. This extra data can help to classify samples more accurately. 

The following sets of features were considered: 
 * A - Metabolites only, 
 * B - Metabolites and SNPs,
 * C - Metabolites, SNPs and covariates,
 * D - Metabolites and covariates.
 
![Sets of features][sets]

[sets]: https://github.com/natacourby/Urinary-metabolomic-phenotyping-for-Alzheimer-pathology/blob/master/images/sets.png "Sets of features"

In addition to the feature set selection problem there was a problem of classification itself. The number of samples per class is imbalanced, especially the cMCI class. 

|                     | AD  | CTL | cMCI | sMCI |
| --------------------|:---:| ---:|-----:|-----:|
| Metabolites only    | 162 | 133 | 36   | 140  |
| Metabolites and SNPs| 121 | 119 | 23   | 80   |

Since decision trees, and as a result RF, are sensitive to class imbalance we decided in addition to the original 4 classes to test different binary classifiers:
 * I - AD/CTL/cMCI/sMCI (4 original classes), 
 * II - AD+cMCI/CTL+sMCI (2 classes when sMCI class is remapped to CTL and cMCI is remapped to AD), 
 * III - AD/CTL (2 classes when MCIs classes are removed).

The following tasks are covered in RF_feature_selection.R script:
1. Choose set and classifier
2. Tune RF parameters for chosen dataset (set and classifier)
3. Run final RF model with importance parameter to rank metabolites

```
source("./scripts/feature_selection.R")
./scripts/RF_feature_selection.R
```

# GWAS Catalog
Intersection with GWAS Catalog SNPs.
1. Search by SNPs (SNP intersection)
2. Search by genomic regions (genomic regions intersection)

```
./scripts/GWAS_Catalog.R
```

# Results for annotated metabolites
Annotated metabolites are availbale in the file: "./results/annotated_metabolites.txt".

1. RF using annotated metabolites for AD/CTL classes and then using the same model for cMCI/sMCI classes for validation since they were not used for training.
2. Logistic regression using annotated metabolites for AD/CTL classes and then using the same model for cMCI/sMCI classes for validation since they were not used for training.

```
./scripts/annotated_metabolites.R
```
