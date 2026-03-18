# Birth Outcomes Analysis
*updated 5/8/2021*

| **Filename**			| **Description** |
|:------------------------------|:----------------|
| combp				| Directory for the comb-p regional analysis code |
| plots				| Directory for the code to create the birth outcomes plots |
| analysis_function.R 		| Performs the phenotype analysis for each birth outcome (LGAZ, WGAZ, HCGAZ, WLGAZ) |
| chr_analysis.R 		| Loads the data and calls the analysis functions |
| combine_results_all.sh	| Removes the header from all of the Chr#_filtered_*.txt files (except Chr1) and concatenates them to produce one large results file for each outcome |
| combine_results_chr.sh	| Combines the results files to produce one results file per chromosome |
| create_scripts.sh		| Produces the individual chromosome R and bash scripts for the phenotype analysis |
| filter_MEs.R			| Subsets ME results from overall results, adds FDR adjusted p-values, and merges restuls with previous ME results |				
| region_data.R			| Subsets the results to the sites within 1,000 bp of the top FDR results (used to make the region plots) |
| remove_QTLs.sh		| Removes the QTL sites from the results |
| top_data.R			| Saves the methylation data for the top FDR results |
| top_hits.R			| Filters *_results.txt (non.zeros>=15, total.NAs<=20, non.ones>=15) and adds FDR adjusted p-values. Also saves the top 1,000 results for each outcome. |

*Note: all text files are tab-delimited*

**Workflow:**
1) create_scripts.sh
2) chr_analysis.R (one script per chromosome)
3) combine_results_chr.sh
4) annotate_function.R (WFMethylation > Annotation directory)
5) remove_QTLs.sh
6) combine_results_all.sh
7) qqplot.R (plots subdirectory)
8) filter_MEs.R
9) top_hits.R
10) region_data.R, plot_data.R
11) scatterplots.R, region_plots.R (plots subdirectory)
12) Prep_data.R, Comb_p_script.sh (combp subdirectory)
13) combp_sites.R (combp subdirectory)
14) combp_plots.R (plots subdirectory)