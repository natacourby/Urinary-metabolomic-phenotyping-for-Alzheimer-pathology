# Requirements 
Code was tested on R version 3.4.3
# Data colllections
Large data files are gzipped and should be unzipped before usage.
## Covariates
/data/covariates/
* Covariates_UHPOS.txt
* Covariates_NMR_Urine.txt
* Covariates_EigenMS_format.txt

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
## 1) EigenMS to remove bias of unknown complexity from MS data from “Metabolomics Data Normalization with EigenMS”, Yuliya V. Karpievitch et al. 
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
```
./scripts/run_mqtl.R "NMR_Urine"
```
# Feature selection with Random Forests
```
./scripts/feature_selection.R
```

# Additional scripts (GWAS Catalog search, plots, etc.)

./scripts/additional/

