# Methylation Analysis Pipeline

| **Folder**		| **Description** |
|:----------------------|:----------------|
| Annotation		| Genecode gene annotation information |
| BirthOutcomes		| Birth outcomes analysis (birth outcome ~ methylation + arm + BMI + covariates) |
| CellTypeAdjustment	| Cellular heterogeneity adjustment using ReFACTor |
| MethylOutcome		| Methylation outcome analysis (methylation ~ arm + BMI + arm x BMI + covariates) |
| QTL			| QTL analysis using QTLTools (identify/remove CpG sites associated with genetic variation) |

**Workflow:**
1) CellTypeAdjustment
2) QTL
3) Annotation
4) MethylOutcome
5) BirthOutcomes