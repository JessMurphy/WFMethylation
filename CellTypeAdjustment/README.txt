workflow:
1) transpose_and_filter.R
2) merge_filtered.sh
3) remove_snps.sh
4) refactor.R

README (updated 3/2/2020)

merge_filtered.sh 		Concatonates the sample_IDs.txt file and the individual chromosome methylation matrices from the 
				data directory to produce one overall mathylation matrix (Edited to exclude the X chromosome matrix.)

refactor.R 			The refactor code downloaded from https://github.com/cozygene/refactor/releases (edited to
	    			accomodate numeric and categorical covariates). It takes a tab deliminated sites by samples
	    			methylation matrix and a tab deliminated samples by covariates text file as inputs and outputs
	    			refactor.out.components.txt and refactor.out.rankedlist.txt.

remove_snps.sh 			Determines the CpG sites that are SNPs and removes them from the methylation matrix 

transpose_and_filter.R	 	For each chromosome, it filters the samples to the 105 samples common to both the
				genetic and methylation data, performs CpG site quality control (removes sites with 
				<15% nonzero reads and >20% NA values), imputes missing values to the mean, and
				transposes the data into sites by samples methylation matrices. It takes the R data
				files (samples by sites) as input files.
