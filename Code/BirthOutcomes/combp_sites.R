############################################################
##  Women First Trial: Guatemala Methylation Analysis (Birth Outcomes)         
##  Written by Jessica Murphy 
##  Last edited on March 29, 2021.
##  This script subsets the epigenome-wide results to the sites within 
##  within the top comb-p regions. It also combines the epigenome-wide 
##  results with the comb-p fdr site results.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/combp/results/")

# read in the top comb-p region results
LGAZ = read.table("combp_LGAZ_results.regions-p.bed")
WGAZ = read.table("combp_WGAZ_results.regions-p.bed")
HCGAZ = read.table("combp_HCGAZ_results.regions-p.bed")
WLGAZ = read.table("combp_WLGAZ_results.regions-p.bed")

# read in the comb-p fdr site results
LGAZ2 = read.table("combp_LGAZ_results.fdr.bed")
WGAZ2 = read.table("combp_WGAZ_results.fdr.bed")
HCGAZ2 = read.table("combp_HCGAZ_results.fdr.bed")
WLGAZ2 = read.table("combp_WLGAZ_results.fdr.bed")

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")

# read in the epigenome-wide results
LGAZ_results = read.table("LGAZ_results_FDR.txt", header=T, sep='\t')
WGAZ_results = read.table("WGAZ_results_FDR.txt", header=T, sep='\t')
HCGAZ_results = read.table("HCGAZ_results_FDR.txt", header=T, sep='\t')
WLGAZ_results = read.table("WLGAZ_results_FDR.txt", header=T, sep='\t')

# define region_sites function
# inputs: combp - a dataframe of comb-p region results 
#         fdr - a dataframe of comb-p fdr site results
#         results - a dataframe of epigenome-wide results
#         name - a character string to name the output file
# output: a tab-separated text file of the results within the comb-p regions
region_sites <- function(combp, fdr, results, name) {
  
  # add the header for the region results
  colnames(combp) = c("chrom", "start", "end", "min_p", "n_probes", "z_p", "z_sidak_p")
  
  # extract the chromosome numbers from the results file
  combp$num = as.numeric(sapply(strsplit(combp$chrom, "r", fixed=T), tail, 1))
  
  # define the region based on the start and end positions
  combp$chrom = as.character(combp$chrom)
  combp$region = paste0(combp$chrom, ":", combp$start, "-", combp$end)

  # add the header for the fdr site results
  colnames(fdr) = c("chrom", "start", "end", "p", "region-p", "region-q")
  fdr$chrom = as.character(fdr$chrom)
  
  # create a vector to store the overall results in
  results.keep = c()
  
  # loop through the chromosomes
  for (i in 1:22){
        
    # subset all of the results to the specified chromosome
    regions = combp %>% filter(num==i)
    sites = results %>% filter(Chr==paste0("chr",i))
    fdr2 = fdr %>% filter(chrom==paste0("chr",i))
    
    # create a vector to store the chromosome-specific results in
    results2 = c()
    
    # loop through the regions in the specified chromosome
    for (j in 1:nrow(regions)){

      # filter the comb-p fdr results & EWAS results to the sites that lie within the region
      keep.sites = sites %>% filter(between(pos, regions$start[j], regions$end[j])) %>% mutate(region=regions$region[j])
      keep.fdr = fdr2 %>% filter(between(start, regions$start[j], regions$end[j]))

      # merge the comb-p fdr results & EWAS results by position
      keep = merge(keep.sites, keep.fdr, by.x="pos", by.y="start")
            
      # store the results per region
      results2 = rbind(results2, keep)
    }
    
    # store the results per chromosome
    results.keep = rbind(results.keep, results2)
     
    print(i)
  }  
  
  # save the top comb-p single site results
  write.table(results.keep, paste0(name, ".txt"), row.names=F, quote=F, sep='\t')
}

# call function
region_sites(LGAZ, LGAZ2, LGAZ_results, "combp_LGAZ_sites")
region_sites(WGAZ, WGAZ2, WGAZ_results, "combp_WGAZ_sites")
region_sites(HCGAZ, HCGAZ2, HCGAZ_results, "combp_HCGAZ_sites")
region_sites(WLGAZ, WLGAZ2, WLGAZ_results, "combp_WLGAZ_sites")