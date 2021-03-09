##################
### Setup File for comb-p
### January 21, 2021
##################

setwd('/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes')
library(data.table)

for(phen in c( 'LGAZ', 'WGAZ', 'WLGAZ')){ #'HCGAZ
  setwd('/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes')
  
re <- fread(paste0('results/', phen, '_results_FDR.txt'), header = TRUE, sep = '\t')

names(re)
re <- re[,c(24,25,25,2)]
head(re)
colnames(re) <- c('chrom', 'start', 'end', 'rawp')

### end should be one more than the start
re$end <- as.numeric(as.character(re$end)) + 1

### Chrom should only be the number, remove the 'chr'
re$chrom <- as.character(re$chrom)
re$chrom <- substr(re$chrom, 4, nchar(re$chrom)) 
re$chrom <- as.numeric(as.character(re$chrom))

table(re$chrom)

re <-re[order(re$chrom,re$start),]

### turn off scientific notation
options(scipen=999)

setwd('combp/FilePrep')
write.table(re, paste0('ToSort_', phen, '4combp.bed'), 
            sep = '\t', quote = FALSE, row.names = FALSE)

}

#### these now need to be sorted in bedtools!
#### See Comb_p_script.sh for instructions

#### the column names need to be chrom start end rawp ()


setwd('/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/combp/FilePrep')

library(data.table)

for(phen in c('HCGAZ', 'LGAZ', 'WGAZ', 'WLGAZ')){
df <- fread(paste0('Sorted_', phen, '4combp.bed'))
colnames(df) <- c('chrom', 'start', 'end', 'rawp')
summary(df$rawp)

options(scipen=999)

head(df)
df$chrom <- paste0('chr', df$chrom)

write.table(df, paste0('Sorted_', phen, '4combp.bed'),
            row.names = FALSE, quote = FALSE,
            sep = '\t')

}




###################
## ME's
###################


setwd('/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/')
library(data.table)
re <- fread('results/EWAS_results_FDR.txt', header = TRUE, sep = '\t')

#### Subset to ME's
re <- re[which(re$ME == 'ME'),]

### Need a separate file for interaction, arm, BMI
### Interaction: from the interaction model
### Arm and BMI - from the no interaction model
names(re)
re <- re[,c(23,24,24,4,5,6)]
head(re)
colnames(re)[1:3] <- c('chrom', 'start', 'end')

### end should be one more than the start
re$end <- as.numeric(as.character(re$end)) + 1

### Chrom should only be the number, remove the 'chr'
re$chrom <- as.character(re$chrom)
re$chrom <- substr(re$chrom, 4, nchar(re$chrom)) 
re$chrom <- as.numeric(as.character(re$chrom))

table(re$chrom)

re <-re[order(re$chrom,re$start),]

re_int <- re[,c(1:4)]
re_arm <- re[,c(1:3,6)]
re_BMI <- re[,c(1:3,5)]

colnames(re_int)[4] <- 'rawp'
colnames(re_arm)[4] <- 'rawp'
colnames(re_BMI)[4] <- 'rawp'

### turn off scientific notation
options(scipen=999)

setwd('combp/EWAS/FilePrep')
write.table(re_int, 'ME_ToSort_interaction4combp.bed', 
            sep = '\t', quote = FALSE, row.names = FALSE)
write.table(re_arm, 'ME_ToSort_arm4combp.bed', 
            sep = '\t', quote = FALSE, row.names = FALSE)
write.table(re_BMI, 'ME_ToSort_bmi4combp.bed', 
            sep = '\t', quote = FALSE, row.names = FALSE)

#### these now need to be sorted in bedtools!
#### See Comb_p_script.sh for instructions

#### the column names need to be chrom start end rawp ()


setwd('/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/combp/EWAS/FilePrep')

library(data.table)
df <- fread('ME_Sorted_arm4combp.bed')
colnames(df) <- c('chrom', 'start', 'end', 'rawp')
summary(df$rawp)

options(scipen=999)

head(df)
df$chrom <- paste0('chr', df$chrom)

write.table(df, 'ME_Sorted_arm4combp.bed', row.names = FALSE, quote = FALSE,
            sep = '\t')


### interaction

df <- fread('ME_Sorted_interaction4combp.bed')
colnames(df) <- c('chrom', 'start', 'end', 'rawp')
summary(df$rawp)

options(scipen=999)

head(df)
df$chrom <- paste0('chr', df$chrom)

write.table(df, 'ME_Sorted_interaction4combp.bed', row.names = FALSE, quote = FALSE,
            sep = '\t')


#### bmi

df <- fread('ME_Sorted_bmi4combp.bed')
colnames(df) <- c('chrom', 'start', 'end', 'rawp')
summary(df$rawp)

options(scipen=999)

head(df)
df$chrom <- paste0('chr', df$chrom)

write.table(df, 'ME_Sorted_bmi4combp.bed', row.names = FALSE, quote = FALSE,
            sep = '\t')






