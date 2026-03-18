############################################################
##  Women First Trial: Guatemala Methylation Analysis (Birth Outcomes)       
##  Written by Jessica Murphy 
##  Last edited on May 4, 2021.
##  This script creates region plots of the beta estimates
##  for the top combp regions annotated to genes.
##  Please send any questions to jessica.murphy@ucdenver.edu
############################################################

# load libraries
library(dplyr)
library(ggplot2)

#setwd("~/RESEARCH/Methylation/results/combp/")
setwd("C:/Users/jimur/OneDrive/Documents/RESEARCH/Methylation/results/combp/")

# read in combp results (filter and add region)
LGAZ.regions = read.table("combp_LGAZ_results.txt", header=T, sep='\t') %>%
  mutate(region=paste0(chrom, ":", start, "-", end)) %>% filter(z_sidak_p<=0.05, n_probes>3)
WGAZ.regions = read.table("combp_WGAZ_results.txt", header=T, sep='\t') %>%
  mutate(region=paste0(chrom, ":", start, "-", end)) %>% filter(z_sidak_p<=0.05, n_probes>3)
HCGAZ.regions = read.table("combp_HCGAZ_results.txt", header=T, sep='\t') %>%
  mutate(region=paste0(chrom, ":", start, "-", end)) %>% filter(z_sidak_p<=0.05, n_probes>3)
WLGAZ.regions = read.table("combp_WLGAZ_results.txt", header=T, sep='\t') %>%
  mutate(region=paste0(chrom, ":", start, "-", end)) %>% filter(z_sidak_p<=0.05, n_probes>3)

# just look at the top five results per outcome annotated to a gene
top.LGAZ.regions = LGAZ.regions %>% filter(genes!='.') %>% arrange(z_sidak_p) %>% slice(1:5)
top.WGAZ.regions = WGAZ.regions %>% filter(genes!='.') %>% arrange(z_sidak_p) %>% slice(1:5)
top.HCGAZ.regions = HCGAZ.regions %>% filter(genes!='.') %>% arrange(z_sidak_p) %>% slice(1:5)
top.WLGAZ.regions = WLGAZ.regions %>% filter(genes!='.') %>% arrange(z_sidak_p) %>% slice(1:5)

# read in the single site results corresponding to the top five regions for each outcome
LGAZ.sites = read.table("combp_LGAZ_sites.txt", header=T, sep='\t') %>% filter(region %in% top.LGAZ.regions$region) %>% 
  select(-chrom, end) %>% mutate(outcome="LGAZ")
WGAZ.sites = read.table("combp_WGAZ_sites.txt", header=T, sep='\t') %>% filter(region %in% top.WGAZ.regions$region) %>% 
  select(-chrom, end) %>% mutate(outcome="WGAZ")
HCGAZ.sites = read.table("combp_HCGAZ_sites.txt", header=T, sep='\t') %>% filter(region %in% top.HCGAZ.regions$region) %>% 
  select(-chrom, end) %>% mutate(outcome="HCGAZ")
WLGAZ.sites = read.table("combp_WLGAZ_sites.txt", header=T, sep='\t') %>% filter(region %in% top.WLGAZ.regions$region) %>% 
  select(-chrom, end) %>% mutate(outcome="WLGAZ")

# define the regions that are present in more than one outcome
overlap = c("chr7:229225-229778", "chr16:88497527-88497986", "chr17:604210-604464",
            "chr10:121593753-121594682", "chr8:6627063-6627284", "chr4:32716247-32716414")

# create an empty list to store the output in
out = med = list()

# loop through the overlapped regions
for (i in overlap){
  
  # separate the region into its start and end positions
  chrom = strsplit(i, split=":")[[1]][1]
  region = strsplit(i, split=":")[[1]][2]
  start = strsplit(region, split="-")[[1]][1]
  end = strsplit(region, split="-")[[1]][2]
  
  # filter the results to the sites within the region
  temp.LGAZ = LGAZ.sites %>% filter(Chr==chrom, pos>=start, pos<=end) 
  temp.WGAZ = WGAZ.sites %>% filter(Chr==chrom, pos>=start, pos<=end) 
  temp.HCGAZ = HCGAZ.sites %>% filter(Chr==chrom, pos>=start, pos<=end) 
  temp.WLGAZ = WLGAZ.sites %>% filter(Chr==chrom, pos>=start, pos<=end) 
  
  # combine the results from all of the outcomes 
  temp = bind_rows(temp.LGAZ, temp.WGAZ, temp.HCGAZ, temp.WLGAZ)
  temp$outcome = factor(temp$outcome, levels=c("LGAZ", "WGAZ", "HCGAZ", "WLGAZ"))
  
  # add a variable to differentiate which sites were included in the combp n_probes
  temp$combp = ifelse(temp$region.q>0.05, "no", "yes")
  temp$combp = factor(temp$combp, levels=c("yes", "no"))
  
  # store the combined results per region
  out[[i]] = temp
  
  # calculate the median beta estimate per outcome for the specified region
  med.temp = temp %>% group_by(outcome) %>% summarize(median=median(methyl.beta))
  
  # store the median values per region
  med[[i]] = med.temp
}

# colorblind friendly palette
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

########## OVERLAPPING OUTCOMES ##########

# plot each region with overlapping outcomes

# separate the data into yes/no for the combp variable (chr7:229225-229778)
no.out1 = out[[1]] %>% filter(combp=="no")
yes.out1 = out[[1]] %>% filter(combp=="yes")

# chr7:229225-229778 (combp excluded sites in grey)
ggplot() +
  geom_point(data=yes.out1, aes(x=pos, y=methyl.beta, col=outcome, fill=outcome, size=methyl.p), alpha=0.5) +
  geom_point(data=no.out1, aes(x=pos, y=methyl.beta, size=methyl.p), col=cbPalette[1], fill=cbPalette[1], alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) + 
  geom_hline(data=med[[1]], aes(yintercept=median, col=outcome)) +
  labs(x="Position", y="Beta Estimate", title=names(out)[[1]]) +
  theme_bw() + scale_fill_manual(values=cbPalette[c(4,6,2)]) +
  scale_color_manual(values=cbPalette[c(4,6,2)])

# chr16:88497527-88497986
ggplot(out[[2]]) +
  geom_point(aes(x=pos, y=methyl.beta, color=outcome, size=methyl.p), alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) +
  geom_hline(data=med[[2]], aes(yintercept=median, col=outcome)) +
  labs(x="Position", y="Beta Estimate", title=names(out)[[2]]) +
  theme_bw() + scale_color_manual(values=cbPalette[c(4,6,7,2)])

# chr17:604210-604464
ggplot(out[[3]]) +
  geom_point(aes(x=pos, y=methyl.beta, color=outcome, size=methyl.p), alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) +
  geom_hline(data=med[[3]], aes(yintercept=median, col=outcome)) +
  labs(x="Position", y="Beta Estimate", title=names(out)[[3]]) +
  theme_bw() + scale_color_manual(values=cbPalette[c(4,6)])

# separate the data into yes/no for the combp variable (chr10:121593753-121594682)
no.out4 = out[[4]] %>% filter(combp=="no")
yes.out4 = out[[4]] %>% filter(combp=="yes")

# chr10:121593753-121594682 (combp excluded sites in grey)
ggplot() +
  geom_point(data=yes.out4, aes(x=pos, y=methyl.beta, fill=outcome, color=outcome, size=methyl.p), alpha=0.5) +
  geom_point(data=no.out4, aes(x=pos, y=methyl.beta, size=methyl.p), col=cbPalette[1], fill=cbPalette[1], alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) + scale_shape_manual(values=c(21, 22)) +
  geom_hline(data=med[[4]], aes(yintercept=median, col=outcome)) +
  labs(x="Position", y="Beta Estimate", title=names(out)[[4]]) +
  theme_bw() + scale_fill_manual(values=cbPalette[c(6,2)]) +
  scale_color_manual(values=cbPalette[c(6,2)])

# chr8:6627063-6627284
ggplot(out[[5]]) +
  geom_point(aes(x=pos, y=methyl.beta, color=outcome, size=methyl.p), alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) +
  geom_hline(data=med[[5]], aes(yintercept=median, col=outcome)) +
  labs(x="Position", y="Beta Estimate", title=names(out)[[5]]) +
  theme_bw() + scale_color_manual(values=cbPalette[c(6,2)])

# chr4:32716247-32716414
ggplot(out[[6]]) +
  geom_point(aes(x=pos, y=methyl.beta, color=outcome, size=methyl.p), alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) +
  geom_hline(data=med[[6]], aes(yintercept=median, col=outcome)) +
  labs(x="Position", y="Beta Estimate", title=names(out)[[6]]) +
  theme_bw() + scale_color_manual(values=cbPalette[c(7,2)])

########## INDIVIDUAL OUTCOME ##########

# chr2:237669618-237669957
ggplot(LGAZ.sites %>% filter(region=="chr2:237669618-237669957") %>% mutate(median=median(methyl.beta))) +
  geom_point(aes(x=pos, y=methyl.beta, size=methyl.p), col=cbPalette[4], alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) + theme_bw() +
  geom_hline(aes(yintercept=median), col=cbPalette[4]) +
  labs(x="Position", y="Beta Estimate", title="chr2:237669618-237669957")

# chr8:48428414-48428476
ggplot(LGAZ.sites %>% filter(region=="chr8:48428414-48428476") %>% mutate(median=median(methyl.beta))) +
  geom_point(aes(x=pos, y=methyl.beta, size=methyl.p), col=cbPalette[4], alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) + theme_bw() +
  geom_hline(aes(yintercept=median), col=cbPalette[4]) +
  labs(x="Position", y="Beta Estimate", title="chr8:48428414-48428476")
  
# chr12:114699849-114700285
ggplot(HCGAZ.sites %>% filter(region=="chr12:114699849-114700285") %>% mutate(median=median(methyl.beta))) +
  geom_point(aes(x=pos, y=methyl.beta, size=methyl.p), col=cbPalette[7], alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) + theme_bw() +
  geom_hline(aes(yintercept=median), col=cbPalette[7]) +
  labs(x="Position", y="Beta Estimate", title="chr12:114699849-114700285")
  
# chr7:102933918-102934080
ggplot(HCGAZ.sites %>% filter(region=="chr7:102933918-102934080") %>% mutate(median=median(methyl.beta))) +
  geom_point(aes(x=pos, y=methyl.beta, size=methyl.p), col=cbPalette[7], alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) + theme_bw() +
  geom_hline(aes(yintercept=median), col=cbPalette[7]) +
  labs(x="Position", y="Beta Estimate", title="chr7:102933918-102934080")

# chr11:68374730-68374870
ggplot(HCGAZ.sites %>% filter(region=="chr11:68374730-68374870") %>% mutate(median=median(methyl.beta))) +
  geom_point(aes(x=pos, y=methyl.beta, size=methyl.p), col=cbPalette[7], alpha=0.5) +
  scale_size("p-value", trans="log10", range=c(14, 0.5), breaks=c(1e-7, 1e-5, 1e-3, 1e-1)) +
  geom_hline(aes(yintercept=0), show.legend=F) + theme_bw() +
  geom_hline(aes(yintercept=median), col=cbPalette[7]) +
  labs(x="Position", y="Beta Estimate", title="chr11:68374730-68374870")
