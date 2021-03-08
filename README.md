# Woman First Trial Methylation Analysis
## Guatemala

This is the repository for the Woman First Trial DNA methylation analysis pipeline.

| **Filename**		| **Description** |
|:----------------------|:----------------|
| Annotation		| Prepare Genecode gene annotation information |
| BirthOutcomes		| Perform birth outcomes analysis (birth outcome ~ methylation + covariates) |
| CellTypeAdjustment	| Adjust for cellular heterogeneity using ReFACTor |
| MethylOutcome		| Perform methylation outcome analysis (methylation ~ arm + BMI + arm x BMI) |
| QTL			| Identify and remove CpG sites associated with genetic variation using QTLTools |

**Workflow:**
1) CellTypeAdjustment
2) QTL
3) Annotation
4) MethylOutcome
5) BirthOutcomes