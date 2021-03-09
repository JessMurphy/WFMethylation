
library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")

# read in top results
top.LGAZ = read.table("top_LGAZ_1000.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.05) 
top.WGAZ = read.table("top_WGAZ_1000.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.05)
top.ME = read.table("ME_results_LGAZ.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.1)

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

  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")
  write.table(out2, paste0(name, ".txt"), row.names=F, quote=F, sep='\t')
}

# call function
top_data(top.LGAZ, "top_LGAZ_data")
top_data(top.WGAZ, "top_WGAZ_data")
top_data(top.ME, "top_ME_data")
