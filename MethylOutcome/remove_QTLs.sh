#!/bin/bash

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/EWAS/results/

# remove the QTL sites for all chromosomes
for i in {1..22..1}
do

# double check the number of lines before filtering
wc -l Chr${i}_annotated.txt

# remove the significant QTL sites
fgrep -wvf /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/Chr${i}_QTL.txt Chr${i}_annotated.txt > Chr${i}_filtered.txt

# double check the number of lines after filtering
wc -l Chr${i}_filtered.txt

done 