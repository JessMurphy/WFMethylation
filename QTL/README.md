workflow:
1) QTL_formatting.R
2) index_bed.sh
3) index_vcf.sh
4) QTL_scripts.sh
5) Chr#_QTL.sh 
6) QTL_pvalues.R (not required)
7) QTL_significant.R
8) QTL_sites.sh

README (updated 1/25/2021)

index_bed.sh		Sorts and zips each bed file in a bed_files directory

index_vcf.sh		Copies the vcf files into a vcf_files directory. It also sorts and zips each vcf file in the directory.

QTL_formatting.R	Reformats the methylation data into bed format for QTLTools. It uses R data files and map files
		  	to output a bed file for each chromosome saved in the bed_files directory.

QTL_pvalues.R		Counts the number of significant CpG sites per chromosome for several p-value thresholds based on
	       		both the nominal and adjusted p-values. It produces a chromosome by threshold table as well as the
	       		nominal and adjusted QQplots for each chromosome.

QTL_scripts.sh		Generates QTL scripts for each chromosome

QTL_significant.R	Filters the significant CpG sites (adjusted p-value < 0.01) from the permutation results 
		   	
QTL_sites.sh		Extracts just the CpG site names from the significant QTL files