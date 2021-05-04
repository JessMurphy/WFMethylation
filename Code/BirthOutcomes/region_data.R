############################################################
##  Women First Trial: Guatemala Methylation Analysis (Birth Outcomes)         
##  Written by Jessica Murphy 
##  Last edited on January 24, 2021.
##  This script subsets the EWAS results to the sites within 1,000
##  bp of the top FDR results (used to make the region plots).
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")

#top = c("chr2.229974407", "chr10.121593892", "chr13.75326388")

# read in top FDR and ME results
top.LGAZ = read.table("top_LGAZ_1000.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.05) %>% select(methyl, Chr, pos)
top.WGAZ = read.table("top_WGAZ_1000.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.05) %>% select(methyl, Chr, pos)
top.ME = read.table("ME_results_LGAZ.txt", header=T, sep='\t') %>% filter(methyl.fdr < 0.1) %>% select(methyl, Chr, pos)

# define region_data function
# inputs: top - a dataframe of top results 
#         name - a character string to name the output file
#         outcome - a character string of the specific birth outcome (LGAS, WGAZ, HCGAZ, or WLGAZ)
# output: a tab-separated text file of the results within 1,000 bp of the top results
region_data <- function(top, outcome, name) {
  
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
    file.name2 = paste0("Chr", top$num[i], "_filtered_", outcome, ".txt")
    
    # every chromosome except Chr1 does not have a header
    if (top$num[i] != 1){
      
      # load the chromosome results 
      results = read.table(file.name2, sep='\t')
      
      # add header back in
      colnames(results) = c('methyl', 'methyl.p', 'arm2.p', 'arm3.p', 'BMI.p', 'age.p', 'expose.p', 'PC1.p',
                            'methyl.beta', 'arm2.beta', 'arm3.beta', 'BMI.beta', 'age.beta', 'expose.beta', 'PC1.beta',
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

# call function (only LGAZ and WGAZ had sites pass the FDR threshold)
region_data(top.LGAZ, "LGAZ", "top_LGAZ_regions")
region_data(top.WGAZ, "WGAZ", "top_WGAZ_regions")
region_data(top.ME, "LGAZ", "top_ME_regions")
