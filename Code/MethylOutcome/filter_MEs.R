
library(dplyr)

# read in Megan's ME results
MEs = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/data/ME_results_final.txt", header=T)

# read in EWAS results
results = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/EWAS_results.txt", header=T, sep='\t')

# subset ME results from EWAS
ME_results = results %>% filter(ME == "ME")

# calculate total NA values (NAs + outliers)
ME_results$total.NAs = ME_results$NAs + ME_results$outliers

# filter ME results
ME_results2 = ME_results %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15)

# calculate how many sites were removed from filtering
print(paste0(nrow(ME_results)-nrow(ME_results2), " sites were removed from filtering"))

# add fdr adjusted p values
ME_results2$BMI.fdr = p.adjust(ME_results2$BMI.p, "fdr")
ME_results2$arm.fdr = p.adjust(ME_results2$arm.p, "fdr")
ME_results2$int.fdr = p.adjust(ME_results2$int.p, "fdr")
ME_results2$BMI.fdr2 = p.adjust(ME_results2$BMI.p2, "fdr")
ME_results2$arm.fdr2 = p.adjust(ME_results2$arm.p2, "fdr")

# merge Megan's results with new results
merged_results = merge(MEs, ME_results2, by.y = "methyl") # 1-37 Megan results

# save results
write.table(ME_results2, "ME_results.txt", row.names=F, quote=F, sep='\t')
write.table(merged_results, "MEs_merged.txt", row.names=F, quote=F, sep='\t')