
library(dplyr)

# read in Megan's ME results
MEs = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/data/ME_results_final.txt", header=T)

# read in birth outcomes results
LGAZ = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/LGAZ_results.txt", header=T, sep='\t')
WGAZ = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/WGAZ_results.txt", header=T, sep='\t')
HCGAZ = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/HCGAZ_results.txt", header=T, sep='\t')
WLGAZ = read.table("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/WLGAZ_results.txt", header=T, sep='\t')

# subset ME results 
LGAZ_ME = LGAZ %>% filter(ME == "ME")
WGAZ_ME = WGAZ %>% filter(ME == "ME")
HCGAZ_ME = HCGAZ %>% filter(ME == "ME")
WLGAZ_ME = WLGAZ %>% filter(ME == "ME")

# calculate total NA values (NAs + outliers)
LGAZ_ME$total.NAs = LGAZ_ME$NAs + LGAZ_ME$outliers
WGAZ_ME$total.NAs = WGAZ_ME$NAs + WGAZ_ME$outliers
HCGAZ_ME$total.NAs = HCGAZ_ME$NAs + HCGAZ_ME$outliers
WLGAZ_ME$total.NAs = WLGAZ_ME$NAs + WLGAZ_ME$outliers

# filter the results
LGAZ_ME2 = LGAZ_ME %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 
WGAZ_ME2 = WGAZ_ME %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 
HCGAZ_ME2 = HCGAZ_ME %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 
WLGAZ_ME2 = WLGAZ_ME %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 

# calculate how many sites were removed from filtering
print(paste0(nrow(LGAZ_ME)-nrow(LGAZ_ME2), " LGAZ sites were removed from filtering"))
print(paste0(nrow(WGAZ_ME)-nrow(WGAZ_ME2), " WGAZ sites were removed from filtering"))
print(paste0(nrow(HCGAZ_ME)-nrow(HCGAZ_ME2), " HCGAZ sites were removed from filtering"))
print(paste0(nrow(WLGAZ_ME)-nrow(WLGAZ_ME2), " WLGAZ sites were removed from filtering"))

# add fdr adjusted p values
LGAZ_ME2$methyl.fdr = p.adjust(LGAZ_ME2$methyl.p, "fdr")
WGAZ_ME2$methyl.fdr = p.adjust(WGAZ_ME2$methyl.p, "fdr")
HCGAZ_ME2$methyl.fdr = p.adjust(HCGAZ_ME2$methyl.p, "fdr")
WLGAZ_ME2$methyl.fdr = p.adjust(WLGAZ_ME2$methyl.p, "fdr")

# merge Megan's results with new results
LGAZ_merged = merge(MEs, LGAZ_ME2, by.y = "methyl") # 1-37 Megan results
WGAZ_merged = merge(MEs, WGAZ_ME2, by.y = "methyl") # 1-37 Megan results
HCGAZ_merged = merge(MEs, HCGAZ_ME2, by.y = "methyl") # 1-37 Megan results
WLGAZ_merged = merge(MEs, WLGAZ_ME2, by.y = "methyl") # 1-37 Megan results

# save results
write.table(LGAZ_ME2, "ME_results_LGAZ.txt", row.names=F, quote=F, sep='\t')
write.table(WGAZ_ME2, "ME_results_WGAZ.txt", row.names=F, quote=F, sep='\t')
write.table(HCGAZ_ME2, "ME_results_HCGAZ.txt", row.names=F, quote=F, sep='\t')
write.table(WLGAZ_ME2, "ME_results_WLGAZ.txt", row.names=F, quote=F, sep='\t')

write.table(LGAZ_merged, "MEs_merged_LGAZ.txt", row.names=F, quote=F, sep='\t')
write.table(WGAZ_merged, "MEs_merged_WGAZ.txt", row.names=F, quote=F, sep='\t')
write.table(HCGAZ_merged, "MEs_merged_HCGAZ.txt", row.names=F, quote=F, sep='\t')
write.table(WLGAZ_merged, "MEs_merged_WLGAZ.txt", row.names=F, quote=F, sep='\t')