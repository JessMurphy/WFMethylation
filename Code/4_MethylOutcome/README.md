# Step 4: Methylation Outcome Analysis

*updated 3/18/2026*

|**Filename**|**Description**|
|-|-|
|4a\_data\_prep.R|Performs CpG site quality control for each chromosome (removes sites with <15% nonzero reads and >20% NAs as well as SNPs) and combines the necessary covariates with the methylation data|
|4b\_create\_scripts.sh|Produces the individual chromosome R and bash scripts for the EWAS analysis|
|4c\_chr\_analysis.R|Loads the data and calls the analysis functions|
|4d\_combine\_results\_chr.sh|Combines the results files to produce one results file per chromosome|
|4e\_remove\_QTLs.sh|Removes the QTL sites from the results|
|4f\_combine\_results\_all.sh|Removes the header from all of the chromosome results files and concatenates them to produce one large results file|
|4g\_filter\_MEs.R|Subsets ME results from overall results, adds FDR adjusted p-values, and merges results with previous ME results|
|4h\_top\_hits.R|Filters EWAS\_results.txt (non.zeros>=15, total.NAs<=20, non.ones>=15) and adds FDR adjusted p-values. Saves the top interaction results with int.fdr <= 0.1 and the top non-interaction results with int.fdr > 0.1 and BMI.p2/arm.p2 < 1E-4. Also saves the top 1,000 results for each term: arm, BMI, armxBMI interaction, arm2, BMI2.|
|4i\_region\_data.R|Subsets the results to the sites within 1,000 bp of the top FDR results (used to make the region plots)|
|4j\_top\_data.R|Saves the methylation data for the top FDR and comb-p results (used for emmeans/emtrends)|
|4k\_means\_function.R|Adds the arm contrasts as well as the trends (interaction model) or marginal means (non-interaction model) to the top results|
|analysis\_function.R|Performs the methylation analysis with interaction and non-interaction models|

*Note: all text files are tab-delimited*

