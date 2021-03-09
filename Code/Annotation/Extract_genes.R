

setwd('/nfs/storage/math/gross-s2/projects/Gencode/')

df <- read.table('gencode.v31.annotation.gtf', sep = '\t')

dff <- df[which(df$V3 == 'gene'),]
dim(dff)

#### separate the semi colon list

dff$V9 <- as.character(dff$V9)

all <- dff[,c(1:8)]
all$gene_id <- '.'
all$gene_type <- '.'
all$gene_name <- '.'
all$level <- '.'
all$hgnc_id <- '.'
all$havana_gene <- '.'

for(i in 1:nrow(dff)){
  
  temp <- strsplit(dff$V9[i], ';')
  
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

all$gene_id = sapply(strsplit(all$gene_id, " ", fixed=T), tail, 1)
all$gene_type = sapply(strsplit(all$gene_type, " ", fixed=T), tail, 1)
all$gene_name = sapply(strsplit(all$gene_name, " ", fixed=T), tail, 1)
all$level = sapply(strsplit(all$level, " ", fixed=T), tail, 1)
all$hgnc_id = sapply(strsplit(all$hgnc_id, " ", fixed=T), tail, 1)
all$havana_gene = sapply(strsplit(all$havana_gene, " ", fixed=T), tail, 1)

length(which(all$gene_id == '.'))
length(which(all$gene_name == '.'))

write.table(all, 'Gencode.v31.annotation.genes_extra_info.txt',
            row.names = FALSE, col.names = FALSE, sep = '\t',
            quote = FALSE)

df <- all[,c(1:9,11)]

write.table(df, 'Gencode.v31.annotation.genes.txt',
            row.names = FALSE, col.names = FALSE, sep = '\t',
            quote = FALSE)

###### Now we have the info to annotate the genes!

