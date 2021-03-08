# Woman First Trial Methylation Analysis

This is the repository for the Woman First Trial methylation analysis pipeline.

The two research questions of interest were:
1) Are study arm and maternal BMI associated with infant DNA methylation at birth?
2) Is infant DNA methylation associated with age-adjusted birth outcomes?

| **Filename**		| **Description** |
|:----------------------|:----------------|
| Annotation		| Prepare Genecode gene annotation information |
| BirthOutcomes		| Perform birth outcomes analysis (birth outcome ~ methylation + arm + BMI + covariates) |
| CellTypeAdjustment	| Adjust for cellular heterogeneity using ReFACTor |
| MethylOutcome		| Perform methylation outcome analysis (methylation ~ arm + BMI + arm x BMI + covariates) |
| QTL			| Identify and remove CpG sites associated with genetic variation using QTLTools |

**Workflow:**
1) CellTypeAdjustment
2) QTL
3) Annotation
4) MethylOutcome
5) BirthOutcomes