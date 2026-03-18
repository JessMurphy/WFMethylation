#!/bin/bash

cd /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/BirthOutcomes/results/

# remove the QTL sites for all chromosomes
for i in {1..22..1}
do

# double check the number of lines before filtering
wc -l Chr${i}_annotated_LGAZ.txt
wc -l Chr${i}_annotated_WGAZ.txt
wc -l Chr${i}_annotated_HCGAZ.txt
wc -l Chr${i}_annotated_WLGAZ.txt

# remove the significant QTL sites
fgrep -wvf /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/Chr${i}_QTL.txt Chr${i}_annotated_LGAZ.txt > Chr${i}_filtered_LGAZ.txt
fgrep -wvf /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/Chr${i}_QTL.txt Chr${i}_annotated_WGAZ.txt > Chr${i}_filtered_WGAZ.txt
fgrep -wvf /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/Chr${i}_QTL.txt Chr${i}_annotated_HCGAZ.txt > Chr${i}_filtered_HCGAZ.txt
fgrep -wvf /nfs/storage/math/gross-s2/projects/guatemala/Guatemala-omics/Methylation/QTL/results/Chr${i}_QTL.txt Chr${i}_annotated_WLGAZ.txt > Chr${i}_filtered_WLGAZ.txt

# double check the number of lines after filtering
wc -l Chr${i}_filtered_LGAZ.txt
wc -l Chr${i}_filtered_WGAZ.txt
wc -l Chr${i}_filtered_HCGAZ.txt
wc -l Chr${i}_filtered_WLGAZ.txt

done 