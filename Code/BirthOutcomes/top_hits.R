library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/")

# read in results
LGAZ = read.table("./LGAZ_results.txt", header=T, sep='\t')
WGAZ = read.table("./WGAZ_results.txt", header=T, sep='\t')
HCGAZ = read.table("./HCGAZ_results.txt", header=T, sep='\t' )
WLGAZ = read.table("./WLGAZ_results.txt", header=T, sep='\t')

# calculate total NA values (NAs + outliers)
LGAZ$total.NAs = LGAZ$NAs + LGAZ$outliers
WGAZ$total.NAs = WGAZ$NAs + WGAZ$outliers
HCGAZ$total.NAs = HCGAZ$NAs + HCGAZ$outliers
WLGAZ$total.NAs = WLGAZ$NAs + WLGAZ$outliers

# filter the results
LGAZ2 = LGAZ %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 
WGAZ2 = WGAZ %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 
HCGAZ2 = HCGAZ %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 
WLGAZ2 = WLGAZ %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 

print(paste0(length(LGAZ2$diff_miss < 0.001), " LGAZ sites with diff.miss p-value < 0.001 after filtering"))
print(paste0(length(WGAZ2$diff_miss < 0.001), " WGAZ sites with diff.miss p-value < 0.001 after filtering"))
print(paste0(length(HCGAZ2$diff_miss < 0.001), " HCGAZ sites with diff.miss p-value < 0.001 after filtering"))
print(paste0(length(WLGAZ2$diff_miss < 0.001), " WLGAZ sites with diff.miss p-value < 0.001 after filtering"))

# calculate how many sites were removed from filtering
print(paste0(nrow(LGAZ)-nrow(LGAZ2), " LGAZ sites were removed from filtering"))
print(paste0(nrow(WGAZ)-nrow(WGAZ2), " WGAZ sites were removed from filtering"))
print(paste0(nrow(HCGAZ)-nrow(HCGAZ2), " HCGAZ sites were removed from filtering"))
print(paste0(nrow(WLGAZ)-nrow(WLGAZ2), " WLGAZ sites were removed from filtering"))

# add FDR adjusted p-values
LGAZ2$methyl.fdr = p.adjust(LGAZ2$methyl.p, "fdr")
WGAZ2$methyl.fdr = p.adjust(WGAZ2$methyl.p, "fdr")
HCGAZ2$methyl.fdr = p.adjust(HCGAZ2$methyl.p, "fdr")
WLGAZ2$methyl.fdr = p.adjust(WLGAZ2$methyl.p, "fdr")

# save FDR results
write.table(LGAZ2, "LGAZ_results_FDR.txt", row.names=F, quote=F, sep='\t')
write.table(WGAZ2, "WGAZ_results_FDR.txt", row.names=F, quote=F, sep='\t')
write.table(HCGAZ2, "HCGAZ_results_FDR.txt", row.names=F, quote=F, sep='\t')
write.table(WLGAZ2, "WLGAZ_results_FDR.txt", row.names=F, quote=F, sep='\t')

# top 1,000 results
top.LGAZ = LGAZ2 %>% top_n(-1000, methyl.p)
top.WGAZ = WGAZ2 %>% top_n(-1000, methyl.p)
top.HCGAZ = HCGAZ2 %>% top_n(-1000, methyl.p)
top.WLGAZ = WLGAZ2 %>% top_n(-1000, methyl.p)

# save top 1,000 results
write.table(top.LGAZ, "top_LGAZ_1000.txt", row.names=F, quote=F, sep='\t')
write.table(top.WGAZ, "top_WGAZ_1000.txt", row.names=F, quote=F, sep='\t')
write.table(top.HCGAZ, "top_HCGAZ_1000.txt", row.names=F, quote=F, sep='\t')
write.table(top.WLGAZ, "top_WLGAZ_1000.txt", row.names=F, quote=F, sep='\t')