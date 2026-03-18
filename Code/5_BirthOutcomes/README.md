# Birth Outcomes Analysis

*updated 5/8/2021*

|**Filename**|**Description**|
|-|-|
|5a\_create\_scripts.sh|Produces the individual chromosome R and bash scripts for the phenotype analysis|
|5b\_chr\_analysis.R|Loads the data and calls the analysis functions|
|5c\_combine\_results\_chr.sh|Combines the results files to produce one results file per chromosome|
|5d\_remove\_QTLs.sh|Removes the QTL sites from the results|
|5e\_combine\_results\_all.sh|Removes the header from all of the Chr#*filtered*\*.txt files (except Chr1) and concatenates them to produce one large results file for each outcome|
|5f\_filter\_MEs.R|Subsets ME results from overall results, adds FDR adjusted p-values, and merges restuls with previous ME results|
|5g\_top\_hits.R|Filters \*\_results.txt (non.zeros>=15, total.NAs<=20, non.ones>=15) and adds FDR adjusted p-values. Also saves the top 1,000 results for each outcome.|
|5h\_region\_data.R|Subsets the results to the sites within 1,000 bp of the top FDR results (used to make the region plots)|
|5i\_top\_data.R|Saves the methylation data for the top FDR results|
|analysis\_function.R|Performs the phenotype analysis for each birth outcome (LGAZ, WGAZ, HCGAZ, WLGAZ)|

*Note: all text files are tab-delimited*

