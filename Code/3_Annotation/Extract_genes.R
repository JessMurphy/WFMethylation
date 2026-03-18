############################################################
##  Women First Trial: Guatemala Methylation Analysis       
##  Written by Jessica Murphy 
##  Last edited on February 5, 2021.
##  This script extracts and formats the gene information from
##  the Gencode annotation file.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

#setwd('/nfs/storage/math/gross-s2/projects/Gencode/')

# read in the Gencode annotation file
df <- read.table('gencode.v31.annotation.gtf', sep = '\t')

# subset the annotation file to just the gene annotations
dff <- df[which(df$V3 == 'gene'),]
dim(dff)

# make sure the 9th column (a semi colon list) is a character string
dff$V9 <- as.character(dff$V9)

# extract the first 8 columns into a new dataframe
all <- dff[,c(1:8)]

# add columns to the dataframe to store the 9th column once separated
all$gene_id <- '.'
all$gene_type <- '.'
all$gene_name <- '.'
all$level <- '.'
all$hgnc_id <- '.'
all$havana_gene <- '.'

# loop through the gene annotations
for(i in 1:nrow(dff)){
  
  # separate the 9th column (semi colon list)
  temp <- strsplit(dff$V9[i], ';')
  
  # some annotations don't contain the same amount of information in the list
  # so determine how much information is stored and save it to the right columns
  if(length(as.vector(temp[[1]])) == 8){
    all[i,9:14] <- as.vector(temp[[1]][-c(6,7)])
  }
  if(length(as.vector(temp[[1]])) == 7){
    all[i,9:14] <- as.vector(temp[[1]][-6])
  }
  if(length(as.vector(temp[[1]])) == 6){
    all[i,9:14] <- as.vector(temp[[1]])
  }
  if(length(as.vector(temp[[1]])) == 5){
    all[i,9:13] <- as.vector(temp[[1]])
  }
  if(length(as.vector(temp[[1]])) == 4){
    all[i,9:12] <- as.vector(temp[[1]])
  }
}

# remove the empty space in the new columns
all$gene_id = sapply(strsplit(all$gene_id, " ", fixed=T), tail, 1)
all$gene_type = sapply(strsplit(all$gene_type, " ", fixed=T), tail, 1)
all$gene_name = sapply(strsplit(all$gene_name, " ", fixed=T), tail, 1)
all$level = sapply(strsplit(all$level, " ", fixed=T), tail, 1)
all$hgnc_id = sapply(strsplit(all$hgnc_id, " ", fixed=T), tail, 1)
all$havana_gene = sapply(strsplit(all$havana_gene, " ", fixed=T), tail, 1)

# determine how many entries are missing a gene id or name
length(which(all$gene_id == '.'))
length(which(all$gene_name == '.'))

# save all of the updated gene annotation information
write.table(all, 'Gencode.v31.annotation.genes_extra_info.txt',
            row.names = FALSE, col.names = FALSE, sep = '\t',
            quote = FALSE)

# subset the data to just the first 8 columns and the gene name and type
df <- all[,c(1:9,11)]

# save the reduced gene annotation information
write.table(df, 'Gencode.v31.annotation.genes.txt',
            row.names = FALSE, col.names = FALSE, sep = '\t',
            quote = FALSE)

###### Now we have the info to annotate the genes!
