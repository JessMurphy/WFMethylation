# README 
*updated 3/2/2021*

| **Filename**			| **Description** |
|:------------------------------|:----------------|
| analysis_function.R 		| Performs the methylation analysis with interaction and non-interaction models |
| annotate_function.R		| Annotates the chromosome results files with overlapping and flanking genes from Gencode |
| chr_analysis.R 		| Loads the data and calls the analysis functions |
| combine_results_all.sh	| Removes the header from all of the chromosome results files and concatenates them to produce one large results file |
| combine_results_chr.sh	| Combines the results files to produce one results file per chromosome |
| create_scripts.sh		| Produces the individual chromosome R and bash scripts for the EWAS analysis |			
| data_prep.R			| Performs CpG site quality control for each chromosome (removes sites with <15% nonzero reads and >20% NAs as well as SNPs) and combines the necessary covariates with the methylation data |	
| filter_MEs.R			| Subsets ME results from overall results, adds FDR adjusted p-values, and merges results with previous ME results |
| int_qqplot.R			| Produces the QQplot for the interaction model (arm, BMI, & armxBMI interaction p-values) |
| noint_qqplot.R		| Produces the QQplot for the non-interaction model (arm & BMI p-values) |
| region_data.R			| Subsets the results to the sites within 1,000 bp of the top FDR results (used to make the region plots) |
| remove_QTLs.sh		| Removes the QTL sites from the results |
| top_data.R			| Saves the methylation data for the top FDR results (used for emmeans/emtrends) |
| top_hits.R			| Filters EWAS_results.txt (non.zeros>=15, total.NAs<=20, non.ones>=15) and adds FDR adjusted p-values. Saves the top interaction results with int.fdr <= 0.1 and the top non-interaction results with int.fdr > 0.1 and BMI.p2/arm.p2 < 1E-4. Also saves the top 1,000 results for each term: arm, BMI, armxBMI interaction, arm2, BMI2. |

*Note: all text files are tab-delimited*

**Workflow:**
1) data_prep.R
2) create_scripts.sh
3) chr#_analysis.R scripts 
4) combine_results_chr.sh
5) annotate_function.R 
6) remove_QTLs.sh
7) combine_results_all.sh
8) int_qqplot.R, noint_qqplot.R
9) filter_MEs.R
10) top_hits.R
11) region_data.R, top_data.R
