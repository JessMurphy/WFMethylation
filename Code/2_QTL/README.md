# Step2: QTL Analysis

*updated 3/18/2026*

|**Filename**|**Description**|
|-|-|
|2a\_QTL\_formatting.R|Reformats the methylation data into bed format for QTLTools. It uses R data files and map files to output a bed file for each chromosome saved in the bed\_files directory.|
|2b\_index\_bed.sh|Sorts and zips each bed file in a bed\_files directory|
|2c\_index\_vcf.sh|Copies the vcf files into a vcf\_files directory. It also sorts and zips each vcf file in the directory.|
|2d\_QTL\_scripts.sh|Generates QTL scripts for each chromosome by copying QTL.sh|
|2e\_Chr#\_QTL.sh (based off QTL.sh)|Run QTLTools for each chromosome|
|2f\_QTL\_pvalues.R|Counts the number of significant CpG sites per chromosome for several p-value thresholds based on both the nominal and adjusted p-values. It produces a chromosome by threshold table as well as the nominal and adjusted QQplots for each chromosome.|
|2g\_QTL\_significant.R|Filters the significant CpG sites (adjusted p-value < 0.01) from the permutation results|
|2h\_QTL\_sites.sh|Extracts just the CpG site names from the significant QTL files|



