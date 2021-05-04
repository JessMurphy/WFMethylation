############################################################
##  Women First Trial: Guatemala Methylation Analysis (Birth Outcomes)         
##  Written by Jessica Murphy 
##  Last edited on August 30, 2020
##  This script performs the analysis for chromosome mychrom.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

# specify chromosome
j = mychrom

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/data/")

# link to analysis function
source("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/analysis/analysis_function.R")

# read in chromosome data (samples by sites)
file.name = paste0("Chr", j, "_data4analysis.rda")
load(file.name)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")

# perform analysis (call analysis function)
pheno(combined, "LGAZ", j)
pheno(combined, "WGAZ", j)
pheno(combined, "HCGAZ", j)
pheno(combined, "WLGAZ", j)