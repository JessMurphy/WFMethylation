
# This script produces the required phenotype files needed to
# run QTLTools.

# The file produced is TAB delimited and each line corresponds to a single 
# molecular phenotype. 

# The first 6 columns are:
# - Chr: Chromosome ID [string]
# - start: Start genomic position of the phenotype [integer, 0-based]
# - end: End genomic position of the phenotype [integer, 1-based]
# - pid: Phenotype ID [string]
# - gid: Phenotype group ID [string]
# - strand: Strand orientation [+/-]

# Then each additional column gives the quantification for a sample. 
# Quantifications are encoded with floating numbers. The file will
# have P lines and N+6 columns where P and N are the numbers of 
# phenotypes and samples, respectively (sites by samples).

# Sample IDs are specified in the header line. This line needs to 
# start with a hash key (i.e. #).

# source: https://qtltools.github.io/qtltools/ 
# (Preparing input files -> Phenotype data [BED])

library(dplyr)

setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/")

# read in vcf sample names
ids = read.table("separated_vcf_ids.txt", header=T)

# loop through chromosomes
for (i in 1:20){
  
  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/data/")
  
  # load the chromosome data (samples by sites)
  file.name = paste0("Chr", i, "_data4analysis.rda")
  load(file.name)
  
  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/Megan_August2019/GenomeWide/Data_files/")
  
  # read in the map file for the chromosome (Chromosome, position, strand, methyl)
  file.name2 = paste0("Chr", i, "_map_file.txt")
  map = read.table(file.name2, header=T)
  
  # remove covariates (leave just wf_id and methylation sites)
  data = combined[,-c(1, 3:29)]
  
  # combine vcf ids with the data (the ids need to be in the same format)
  data = merge(ids, data, by="wf_id")
  data = data %>% select(-c(prefix, wf_id))
  
  # define required variables for bed format
  start = sapply(map$Position, function(x) x-1)
  end = map$Position
  strand = sapply(strsplit(as.character(map$Strand), "G"), tail, 1)
  Chr = sapply(strsplit(as.character(map$Chromosome), "r"), tail, 1) #need just the # to correspond with the vcf files
  
  map2 = map %>% mutate(start, end, strand, Chr) %>% 
    select(Chr, start, end, pid=methyl, gid=methyl, strand)
    
  # convert data into a matrix and transpose it
  data.mat = as.matrix(data)
  data.mat = t(data.mat)
  
  # add sample IDs as column names (remove from first row) & add CpG site names as first column
  colnames(data.mat) = data.mat[1,]
  data.mat = data.mat[-1,]
  data.mat2 = cbind(row.names(data.mat), data.mat)
  data.mat2 = as.data.frame(data.mat2)
  colnames(data.mat2)[1] = "pid"
  
  # combine the map data with the methylation matrix
  data.out = merge(map2, data.mat2, by="pid", all.y=T)
  data.out = data.out %>% select(Chr, start, end, pid, everything())
  
  # add required header to the data
  colnames(data.out)[1] = "#chr"
  
  setwd("/nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/bed_files/")
  out.name = paste0("Chr", i, ".bed")
  
  # save methylation data
  write.table(data.out, out.name, sep="\t", quote=F, row.names=F, col.names=T)
  
  # print chromosome #
  print(i)
}
