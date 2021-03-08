workflow:
1) create_scripts.sh
2) Chr#_analysis.R scripts 
3) combine_results_chr.sh
4) annotate_function.R
5) remove_QTLs.sh
6) combine_results_all.sh
7) QQplot.R
8) filter_MEs.R
9) top_hits.R
10) region_data.R, plot_data.R

README (updated 3/2/2021)

*Note: all text files are tab-delimited
 
analysis_function.R 		Performs the phenotype analysis for each birth outcome (LGAZ, WGAZ, HCGAZ, WLGAZ)

annotate_function.R		Annotates the chromosome results files with overlapping and flanking genes from Gencode

Chr_analysis.R 			Loads the data and calls the analysis functions
    
combine_results_all.sh		Removes the header from all of the Chr#_filtered_*.txt files (except Chr1) and
		    		concatenates them to produce one large results file for each outcome

combine_results_chr.sh		Combines the results files to produce one results file per chromosome 

filter_MEs.R			Subsets ME results from overall results and adds FDR adjusted p-values 
				Merges with previous ME results

QQplot.R			Produces the QQplot for the birth outcome results (LGAZ, WGAZ, HCGAZ, and WLGAZ)

create_scripts.sh		Produces the individual chromosome R and bash scripts for the phenotype analysis
				
top_hits.R			Filters *_results.txt (non.zeros>=15, total.NAs<=20, non.ones>=15) and adds 
				FDR adjusted p-values. Also saves the top 1,000 results for each outcome.
				
region_data.R			Subsets the results to the sites within 1,000 bp of the top FDR results
				(used to make the region plots) 

remove_QTLs.sh			Removes the QTL sites from the results 

top_data.R			Saves the methylation data for the top FDR results