
library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/")

# read in top results
ME.results = read.table("ME_results.txt", header=T, sep='\t')
top.int.fdr = read.table("top_int_fdr.txt", header=T, sep='\t') 
top.arm2.fdr = read.table("top_arm2_fdr.txt", header=T, sep='\t') %>% filter(arm.fdr2 <= 0.05)
top.BMI2.fdr = read.table("top_BMI2_fdr.txt", header=T, sep='\t') %>% filter(BMI.fdr2 <= 0.05)

top_data <- function(results, name) {
  
  results$num = as.numeric(sapply(strsplit(results$Chr, "r", fixed=T), tail, 1))
  results$num = as.factor(results$num)
  
  nums = levels(results$num)
  
  out = data.frame(rep(1, 105))
  
  for (j in 1:length(nums)){
    
    setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/data/")
       
    # load the chromosome data (samples by sites)
    file.name = paste0("Chr", nums[j], "_data4analysis.rda")
    load(file.name)
    
    temp = results %>% filter(num==nums[j])
    data = combined %>% select(temp$methyl)
    
    out = cbind(out, data)
  }
  
  out = out[,-1]
  covs = combined[,1:29]
  out2 = cbind(covs, out)

  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/")
  write.table(out2, paste0(name, ".txt"), row.names=F, quote=F, sep='\t')
}

# call function
top_data(ME.results, "ME_data")
top_data(top.int.fdr, "top_int_data")
top_data(top.arm2.fdr, "top_arm2_data")
top_data(top.BMI2.fdr, "top_BMI2_data")