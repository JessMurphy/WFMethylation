############################################################
##  Women First Trial: Guatemala Methylation Analysis (EWAS)         
##  Written by Jessica Murphy 
##  Last edited on August 24, 2020
##  This script performs CpG site quality control for each chromosome and
##  combines the necessary covariates with the methylation data. An EWAS
##  will then be performed for each chromosome (methylation_analysis.R).
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/Megan_August2019/GenomeWide/Data_files/RDA/")

#################### FUNCTIONS ####################

# define function for quality control
# names: CpG site names
# data: chromosome data

quality <- function(names, data) {
  
  nonzero.reads = NA.values = c()
  
  # number of nonzero reads and missing values for each CpG site
  for(i in 2:length(data)){ #loop through the CpG sites
    reads = length(which(data[,i]>0)) 
    nonzero.reads = c(nonzero.reads, reads)
    NAs = sum(is.na(data[,i]))
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


#################### DATA ####################

# read in covariates
covs = read.delim("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/CellTypeAdjustment/metadata.txt")

# read in refactor components
refactor = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/CellTypeAdjustment/refactor.out.components.txt")
colnames(refactor) = c("refactor1", "refactor2")

# read in SNPs
snps = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/CellTypeAdjustment/removed_snps.txt")

# extract sample IDs
samples = as.character(covs[,1]) 

removed = c()

# loop through chromosomes
for (i in 1:22){
  
  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/Megan_August2019/GenomeWide/Data_files/RDA/")
  
  # load the chromosome data (samples by sites)
  file.name = paste0("Chr", i, "_data4analysis.rda")
  load(file.name)
  
  # filter samples and remove variables (leave infant_id and CpG sites)
  data = df %>% filter(infant_id %in% samples) %>%
    select(-c(id, wf_id, bord, sitename, multb, sex, primary_outcome, cluster, arm,
              BMI_calc2, zlen_wlms_24h, seq.batch, Library.Prep.Round))
  
  # remove SNPs
  data2 = data[,!(colnames(data) %in% snps[,1])]
  
  # extract CpG site names
  sites = names(data2)[2:length(data2)]
  
  # perform quality control
  CpGs = quality(sites, data2)
  data.qc = data2[,!(colnames(data2) %in% CpGs)]
  
  # combine covariates and refactor components
  covs2 = bind_cols(covs, refactor)
  
  # add additional covariates to methylation data
  combined = inner_join(covs2, data.qc, by="infant_id")
  
  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/data/")
  out.name = paste0("Chr", i, "_data4analysis.rda")
  
  # save methylation data
  save(combined, file=out.name)
  
  # store number of CpG sites removed per chromosome
  num = length(CpGs)
  removed = c(removed, num)
  
  # print how many sites were removed from each chromosome
  print(paste0(num, " sites were removed from Chromosome ", i))
}

# print the total number of sites removed from all chromosomes
print(paste0(sum(removed), " total sites were removed from all chromosomes"))
