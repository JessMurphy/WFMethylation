############################################################
##  Women First Trial: Guatemala Methylation Analysis (Methylation Outcome)         
##  Written by Jessica Murphy 
##  Last edited on May 3, 2021.
##  This script saves the methylation data for the top FDR and
##  comb-p results (used for emmeans/emtrends).
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/")

# read in top FDR and ME results
ME.results = read.table("ME_results.txt", header=T, sep='\t')
top.int.fdr = read.table("top_int_fdr.txt", header=T, sep='\t') 
top.arm2.fdr = read.table("top_arm2_fdr.txt", header=T, sep='\t') %>% filter(arm.fdr2 <= 0.05)
top.BMI2.fdr = read.table("top_BMI2_fdr.txt", header=T, sep='\t') %>% filter(BMI.fdr2 <= 0.05)

# read in combp results
arm.sites = read.table("combp_arm_sites.txt", header=T, sep='\t')
bmi.sites = read.table("combp_bmi_sites.txt", header=T, sep='\t')
int.sites = read.table("combp_int_sites.txt", header=T, sep='\t')

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
    load(file.name) #combined
    
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
  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/")
  write.table(out2, paste0(name, ".txt"), row.names=F, quote=F, sep='\t')
}

# call function for the FDR and ME results
top_data(ME.results, "ME_data")
top_data(top.int.fdr, "top_int_data")
top_data(top.arm2.fdr, "top_arm2_data")
top_data(top.BMI2.fdr, "top_BMI2_data")

# call function for the comb-p results
top_data(arm.sites, "arm_combp_data")
top_data(bmi.sites, "bmi_combp_data")
top_data(int.sites, "int_combp_data")