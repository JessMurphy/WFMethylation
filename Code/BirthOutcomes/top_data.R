############################################################
##  Women First Trial: Guatemala Methylation Analysis (Birth Outcomes)         
##  Written by Jessica Murphy 
##  Last edited on January 24, 2021.
##  This script saves the methylation data for the top FDR results
##  (used for emmeans/emtrends).
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")

# read in FDR and ME results
top.LGAZ = read.table("top_LGAZ_1000.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.05) 
top.WGAZ = read.table("top_WGAZ_1000.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.05)
top.ME = read.table("ME_results_LGAZ.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.1)

# define top_data function
# inputs: results - a dataframe of results 
#         name - a character string to name the output file
# output: a tab-separated text file of the data for the top results
top_data <- function(results, name) {
  
  # extract the chromosome numbers from the results file
  results$num = as.numeric(sapply(strsplit(results$Chr, "r", fixed=T), tail, 1))
  results$num = as.factor(results$num)
  nums = levels(results$num)
  
  # create data frame to store data in
  out = data.frame(rep(1, 105))
  
  # loop through the chromosomes
  for (j in 1:length(nums)){
    
    setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/data/")
       
    # load the chromosome data (samples by sites)
    file.name = paste0("Chr", nums[j], "_data4analysis.rda")
    load(file.name)
    
    # subset the results based on the specified chromosome
    temp = results %>% filter(num==nums[j])
    
    # subset the methylation data for the top CpG sites
    data = combined %>% select(temp$methyl)
    
    # store the data per chromosome
    out = cbind(out, data)
  }
  
  # remove the initial column of ones
  out = out[,-1]
  
  # extract the covariates
  covs = combined[,1:29]
  
  # merge the covariates with the combined data
  out2 = cbind(covs, out)

  # save the combined data
  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")
  write.table(out2, paste0(name, ".txt"), row.names=F, quote=F, sep='\t')
}

# call function
top_data(top.LGAZ, "top_LGAZ_data")
top_data(top.WGAZ, "top_WGAZ_data")
top_data(top.ME, "top_ME_data")
