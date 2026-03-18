# Methylation Outcome Analysis
*updated 5/8/2021*

| **Filename**			| **Description** |
|:------------------------------|:----------------|
| combp				| Directory for the comb-p regional analysis code |
| plots				| Directory for the code to create the methylation outcome plots |
| analysis_function.R 		| Performs the methylation analysis with interaction and non-interaction models |
| chr_analysis.R 		| Loads the data and calls the analysis functions |
| combine_results_all.sh	| Removes the header from all of the chromosome results files and concatenates them to produce one large results file |
| combine_results_chr.sh	| Combines the results files to produce one results file per chromosome |
| create_scripts.sh		| Produces the individual chromosome R and bash scripts for the EWAS analysis |			
| data_prep.R			| Performs CpG site quality control for each chromosome (removes sites with <15% nonzero reads and >20% NAs as well as SNPs) and combines the necessary covariates with the methylation data |	
| filter_MEs.R			| Subsets ME results from overall results, adds FDR adjusted p-values, and merges results with previous ME results |
| means_function.R		| Adds the arm contrasts as well as the trends (interaction model) or marginal means (non-interaction model) to the top results |
| region_data.R			| Subsets the results to the sites within 1,000 bp of the top FDR results (used to make the region plots) |
| remove_QTLs.sh		| Removes the QTL sites from the results |
| top_data.R			| Saves the methylation data for the top FDR and comb-p results (used for emmeans/emtrends) |
| top_hits.R			| Filters EWAS_results.txt (non.zeros>=15, total.NAs<=20, non.ones>=15) and adds FDR adjusted p-values. Saves the top interaction results with int.fdr <= 0.1 and the top non-interaction results with int.fdr > 0.1 and BMI.p2/arm.p2 < 1E-4. Also saves the top 1,000 results for each term: arm, BMI, armxBMI interaction, arm2, BMI2. |

*Note: all text files are tab-delimited*

**Workflow:**
1) data_prep.R
2) create_scripts.sh
3) chr_analysis.R (one script per chromosome) 
4) combine_results_chr.sh
5) annotate_function.R (WFMethylation > Annotation directory)
6) remove_QTLs.sh
7) combine_results_all.sh
8) int_qqplot.R, noint_qqplot.R (plots subdirectory)
9) filter_MEs.R
10) top_hits.R
11) region_data.R, top_data.R
12) means_function.R
13) interaction_plots.R, region_plots.R (plots subdirectory)
14) Prep_data.R, Comb_p_script.sh (combp subdirectory)
15) combp_sites.R (combp subdirectory)
16) combp_plots.R (plots subdirectory)
