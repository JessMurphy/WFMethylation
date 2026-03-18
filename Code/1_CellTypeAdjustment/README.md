# Cell Type Adjustment

*updated 3/18/2026*

|**Filename**|**Description**|
|-|-|
|1\_transpose\_and\_filter.R|For each chromosome, it filters the samples to the 105 samples common to both the genetic and methylation data, performs CpG site quality control (removes sites with <15% nonzero reads and >20% NA values), imputes missing values to the mean, and transposes the data into sites by samples methylation matrices. It takes the R data files (samples by sites) as input files.|
|2\_merge\_filtered.sh|Concatenates the sample\_IDs.txt file and the individual chromosome methylation matrices from the data directory to produce one overall methylation matrix (Edited to exclude the X chromosome matrix.)|
|3\_remove\_snps.sh|Determines the CpG sites that are SNPs and removes them from the methylation matrix|
|4\_refactor.R|The refactor code downloaded from https://github.com/cozygene/refactor/releases (edited to accommodate numeric and categorical covariates). It takes a tab deliminited sites by samples methylation matrix and a tab delimited samples by covariates text file as inputs and outputs refactor.out.components.txt and refactor.out.rankedlist.txt.|



