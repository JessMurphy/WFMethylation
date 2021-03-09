
library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/")

results = read.table("./EWAS_results.txt", header = T, sep='\t')

# calculate total NA values (NAs + outliers)
results$total.NAs = results$NAs + results$outliers

# filter the results
results2 = results %>% filter(non.zeros>=15, total.NAs<=20, non.ones>=15) 

print(paste0(length(results2$diff_miss < 0.001), " sites with diff.miss p-value < 0.001 after filtering"))

# calculate how many sites were removed from filtering
print(paste0(nrow(results)-nrow(results2), " sites were removed from filtering"))

# add FDR adjusted p-values
results2$BMI.fdr = p.adjust(results2$BMI.p, "fdr")
results2$arm.fdr = p.adjust(results2$arm.p, "fdr")
results2$int.fdr = p.adjust(results2$int.p, "fdr")
results2$BMI.fdr2 = p.adjust(results2$BMI.p2, "fdr")
results2$arm.fdr2 = p.adjust(results2$arm.p2, "fdr")

write.table(results2, "EWAS_results_FDR.txt", row.names=F, quote=F, sep='\t')

# filter top interaction results
top.int.fdr = results2 %>% filter(int.fdr <= 0.1) 
top.int = results2 %>% top_n(-1000, int.p)
top.arm = results2 %>% top_n(-1000, arm.p)
top.BMI = results2 %>% top_n(-1000, BMI.p)

# save top interaction results
write.table(top.int.fdr, "top_int_fdr.txt", row.names=F, quote=F, sep='\t')
write.table(top.int, "top_int_1000.txt", row.names=F, quote=F, sep='\t')
write.table(top.arm, "top_arm_1000.txt", row.names=F, quote=F, sep='\t')
write.table(top.BMI, "top_BMI_1000.txt", row.names=F, quote=F, sep='\t')

# filter top non-interaction results
top.BMI2.fdr = results2 %>% filter(int.fdr > 0.1) %>% filter(BMI.p2 < 1E-4)
top.arm2.fdr = results2 %>% filter(int.fdr > 0.1) %>% filter(arm.p2 < 1E-4)
top.arm2 = results2 %>% top_n(-1000, arm.p2)
top.BMI2 = results2 %>% top_n(-1000, BMI.p2)

# save top non-interaction results
write.table(top.BMI2.fdr, "top_BMI2_fdr.txt", row.names=F, quote=F, sep='\t')
write.table(top.arm2.fdr, "top_arm2_fdr.txt", row.names=F, quote=F, sep='\t')
write.table(top.arm2, "top_arm2_1000.txt", row.names=F, quote=F, sep='\t')
write.table(top.BMI2, "top_BMI2_1000.txt", row.names=F, quote=F, sep='\t')

