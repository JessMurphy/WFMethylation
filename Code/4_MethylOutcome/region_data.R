############################################################
##  Women First Trial: Guatemala Methylation Analysis (Methylation Outcome)         
##  Written by Jessica Murphy 
##  Last edited on January 19, 2021.
##  This script subsets the EWAS results to the sites within 1,000
##  bp of the top FDR results (used to make the region plots).
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/")

# read in top FDR and ME results
MEs = read.table("ME_results.txt", header=T, sep='\t') 
top.int.fdr = read.table("top_int_fdr.txt", header=T, sep='\t') 
top.arm2.fdr = read.table("top_arm2_fdr.txt", header=T, sep='\t') 
top.BMI2.fdr = read.table("top_BMI2_fdr.txt", header=T, sep='\t')  

# filter the FDR and ME results (0.1 FDR for ME, 0.05 FDR otherwise)
top.MEs = MEs %>% filter(int.fdr < 0.1 | arm.fdr < 0.1 | BMI.fdr < 0.1 | arm.fdr2 < 0.1 | BMI.fdr2 < 0.1) %>% select(methyl, Chr, pos)
top.int = top.int.fdr %>% filter(int.fdr < 0.05 | arm.fdr < 0.05 | BMI.fdr < 0.05) %>% select(methyl, Chr, pos)
top.arm2 = top.arm2.fdr %>% arrange(arm.p2) %>% filter(arm.fdr2 < 0.05) %>% top_n(-10, arm.fdr2) %>% select(methyl, Chr, pos)
top.BMI2 = top.BMI2.fdr %>% arrange(BMI.p2) %>% filter(BMI.fdr2 < 0.05) %>% select(methyl, Chr, pos)

#top = rbind(top.arm2, top.BMI2, top.int)

# define region_data function
# inputs: top - a dataframe of top results 
#         name - a character string to name the output file
# output: a tab-separated text file of the results within 1,000 bp of the top results

region_data <- function(top, name) {
  
  # extract the chromosome numbers from the results file
  top$num = as.numeric(sapply(strsplit(top$Chr, "r", fixed=T), tail, 1))
  
  # define the region around each top results
  top$pos = as.numeric(top$pos)
  top$start = top$pos-1000
  top$end = top$pos+1000
  
  # create data frame to store data in & vector to store results in
  out = data.frame(rep(1, 105))
  results.keep = c()
  
  # loop through the top sites
  for (i in 1:nrow(top)){
    
    # chromosome results filename
    file.name2 = paste0("Chr", top$num[i], "_filtered.txt")
    
    # every chromosome except Chr1 does not have a header
    if (top$num[i] != 1){
      
      # load the chromosome results 
      results = read.table(file.name2, sep='\t')
      
      # add header back in
      colnames(results) = c('methyl', 'BMI.p', 'arm.p', 'int.p', 'BMI.p2', 'arm.p2', 'arm2.beta1', 'arm3.beta1',
                            'BMI.beta1', 'arm2xbmi.beta1', 'arm3xbmi.beta1', 'arm2.beta2', 'arm3.beta2', 'BMI.beta2',
                            "NAs", "zeros", "non.zeros", "ones", "non.ones", "outliers", "diff.miss", "ME",
                            "Chr", "pos", "genes", "genes_flanking", "genes_type")
    } else {
      
      # load the chromosome results 
      results = read.table(file.name2, header=T, sep='\t')
    }
    
    # determine which sites lie within the region
    keep = between(results$pos, top$start[i], top$end[i])
    sites = results$methyl[keep]
    
    # subset the results to the region sites
    results2 = results %>% filter(methyl %in% sites)
    
    # keep track of the top result associated with each region
    bin = rep(top$methyl[i], length(sites))
    results2$bin = bin
    
    # store the results per top site
    results.keep = rbind(results.keep, results2)
    
    print(i)
  }
  
  # save the top region results
  write.table(results.keep, paste0(name, ".txt"), row.names=F, quote=F, sep='\t')
}

# call function
region_data(top.MEs, "top_ME_regions")
region_data(top.int, "top_int_regions")
region_data(top.arm2, "top_arm2_regions")
region_data(top.BMI2, "top_BMI2_regions")
