############################################################
##  Women First Trial: Guatemala Methylation Analysis (Cell Type Adjustment)         
##  Written by Jessica Murphy 
##  Last edited on July 10, 2020
##  This script performs CpG site quality control for each chromosome and
##  imputes missing methylation values to the mean. It also transposes the
##  data, resulting in a CpG sites by samples methylation matrix for each 
##  chromosome. These transposed matrices will then be combined into one
##  methylation matrix (merge_filtered.sh) to be used in the ReFACTor
##  method for cell type adjustment (refactor.R).
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/Megan_August2019/GenomeWide/Data_files/RDA/")

# read in covariates
data = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/CellTypeAdjustment/refactor_metadata.txt")

# extract sample IDs
samples = as.character(data[,1]) # with ID_ prefix
samples2 = sapply(strsplit(samples, "_"), tail, 1) # without ID_ prefix

# define function for quality control
# names: CpG site names
# data: chromosome data

quality <- function(names, data) {
  
  nonzero.reads = NA.values = c()
  
  # number of nonzero reads and missing values for each CpG site
  for(i in 2:length(data)){ #loop through the CpG sites
    reads = length(which(data[i]>0)) 
    nonzero.reads = c(nonzero.reads, reads)
    NAs = sum(is.na(data[i]))
    NA.values = c(NA.values, NAs)
  }
  # percentage of nonzero reads and NA values
  NA.perc = NA.values/nrow(data)
  nonzero.perc = nonzero.reads/nrow(data)
  
  # sites with <15% nonzero reads
  zeros = names[which(nonzero.perc<.15)]
  
  # sites with >20% NAs 
  missing = names[which(NA.perc>.2)]
  
  # combine nonzero reads or NAs
  joint = union(zeros, missing)
  
  return(joint)
}

# define function for imputation
# data: chromosome data

impute <- function(data) {
  
  for(i in 2:length(data)){ #loop through the CpG sites
    
    # replace NA values with the mean methylation value at each site
    data[,i][is.na(data[,i])] <- mean(data[,i], na.rm=T)
  }
  return(data)
}

chrom = c(1:22, "X")

removed = c()

# loop through chromosomes
for (i in chrom){
  
  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/Megan_August2019/GenomeWide/Data_files/RDA/")
  
  # load the chromosome data
  file.name = paste0("Chr", i, "_data4analysis.rda")
  load(file.name)
  
  # filter samples and remove variables (leave infant_id and CpG sites)
  data = df %>% filter(infant_id %in% samples2) %>%
    select(-c(id, wf_id, bord, sitename, multb, sex, primary_outcome, cluster, arm,
                          BMI_calc2, zlen_wlms_24h, seq.batch, Library.Prep.Round))
  
  # extract CpG site names
  sites = names(data)[2:length(data)]
  
  # perform quality control
  CpGs = quality(sites, data)
  data.qc = data[,!(colnames(data) %in% CpGs)]
  
  # perform imputation
  data.qc.imp = impute(data.qc)
  
  # convert data into a matrix and transpose it
  data.mat = as.matrix(data.qc.imp)
  data.mat = t(data.mat)
  
  # remove first row of sample IDs (to make the datasets easier to combine later on)
  data.mat = data.mat[-1,]
  
  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/CellTypeAdjustment/data/")
  out.name = paste0("Chr", i, "_filtered.txt")
  
  # save transposed methylation data
  write.table(data.mat, out.name, row.names=T, col.names=F, quote=F, sep='\t')
  
  # store number of CpG sites removed per chromosome
  num = length(CpGs)
  removed = c(removed, num)
  
  # print how many sites were removed from each chromosome
  print(paste0(num, " sites were removed from Chromosome ", i))
}

# print the total number of sites removed from all chromosomes
print(paste0(sum(removed), " total sites were removed from all chromosomes"))

# save infant IDs
samples = c("infant_id", samples)
write.table(t(samples), "sample_IDs.txt", row.names=F, col.names=F, quote=F, sep='\t')
