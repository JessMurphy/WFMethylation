# Methylation Outcome Results

The table below describes the variables in each results file.

| **No.** | **Name** 	| **Description** |
|:-------:|:------------|:----------------|
| 1 | methyl		| chr.pos of the CpG |
| 2 | BMI.p		| BMI p-value (interaction model) |
| 3 | arm.p		| arm p-value (interaction model) |
| 4 | int.p		| armxBMI interaction p-value (interaction model) |
| 5 | BMI.p2		| BMI p-value (non-interaction model) |
| 6 | arm.p2		| arm p-value (non-interaction model) |
| 7 | arm2.beta1	| arm2 beta coefficient (interaction model) |
| 8 | arm3.beta1	| arm3 beta coefficient (interaction model) |
| 9 | BMI.beta1		| BMI beta coefficient (interaction model) |
| 10 | arm2xbmi.beta1	| arm2xBMI beta coefficient (interaction model) |
| 11 | arm3xbmi.beta1	| arm3xBMI beta coefficient (interaction model) |
| 12 | arm2.beta2	| arm2 beta coefficient (non-interaction model) |
| 13 | arm3.beta2	| arm3 beta coefficient (non-interaction model) |
| 14 | BMI.beta2	| BMI beta coefficient (non-interaction model) |
| 15 | NAs		| number of samples with NA methylation values at that CpG (excluding outliers) |
| 16 | zeros		| number of samples with 0% methylation at that CpG |
| 17 | non.zeros	| number of samples with nonzero methylation values at that CpG |
| 18 | ones		| number of samples with 100% methylation at that CpG |
| 19 | non.ones		| number of samples with non-one methylation values at that CpG |
| 20 | outliers		| number of samples considered an outlier at that CpG (> Q3 + 3*IQR or < Q1 - 3*IQR) |
| 21 | diff.miss	| pvalue for the calculation of differential missingness |
| 22 | ME		| metastable epiallele identifier (. = no, ME = yes) |
| 23 | Chr		| chromsome number |
| 24 | pos		| chromosome position |
| 25 | genes		| genes that directly overlap the CpG |
| 26 | genes_flanking	| genes within 5kb of the CpG |
| 27 | gene_type	| gene type (Ensembl or HAVANA) of the genes listed in the genes_flanking column |
| 28 | BMI.fdr		| FDR adjusted p-value for BMI (interaction model) |
| 29 | arm.fdr		| FDR adjusted p-value for arm (interaction model) |
| 30 | int.fdr		| FDR adjusted p-value for armxBMI interaction (interaction model) |
| 31 | BMI.fdr2		| FDR adjusted p-value for BMI (non-interaction model) |
| 32 | arm.fdr2		| FDR adjusted p-value for arm (non-interaction model) |
